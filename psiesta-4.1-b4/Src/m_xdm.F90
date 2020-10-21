      MODULE m_xdm
! ^LML 2019-06-11
! Please cite L. M. LeBlanc, J. A. Weatherby, A. Otero-de-la-Roza and E. R. Johnson
! J. Chem. Theory Comput. (2018) 14, 5715-5724. DOI: https://pubs.acs.org/doi/full/10.1021/acs.jctc.8b00797
! when using the XDM dispersion model in SIESTA
      use precision, only:dp, grid_p
      use fdf
      use sys, only: die
      use parallel, only: IOnode, Node, Nodes
      use alloc, only: re_alloc, de_alloc

      IMPLICIT NONE

      PUBLIC :: initxdm, coefxdm, enfrxdm
      PRIVATE :: inversebr

      logical,save,PUBLIC :: xdm_bool, xdm_boolc, coef_once, coef_savd
      integer,save,PUBLIC :: xdm_mesh
      real(grid_p),allocatable,save,PUBLIC :: rho_r(:,:,:) ! nsp, np, nspin
      real(grid_p),allocatable,save,PUBLIC :: grrho_r(:,:,:) ! nsp, np, nspin
      real(grid_p),allocatable,save,PUBLIC :: laprho_r(:,:,:) ! nsp, np, npsin
      real(grid_p),allocatable,save,PUBLIC :: tau_r(:,:,:) ! nsp, np, nspin

      real(dp),save,PRIVATE :: xdm_a1
      real(dp),save,PRIVATE :: xdm_a2
      real(dp),save,PRIVATE :: xdm_sc
      real(dp),save,PRIVATE :: xdm_cutoff
      real(dp),allocatable,save,PRIVATE :: cx(:,:,:) ! na_u, na_u, 2:4
      real(dp),allocatable,save,PRIVATE :: rvdw(:,:) ! na_u, na_u

      CONTAINS

      SUBROUTINE initxdm(nspin)

      integer,intent(in) :: nspin

      real(dp) :: default = 1.0000_dp ! for XDM a1 and a2 parameters

!--------------------------------------------------------------------------- BEGIN

      ! Set relevant control switches for XDM initialization
      if (fdf_boolean('XDM',.false.)) then

        if (nspin .gt. 2) then 
          call die('XDM not implemented for non-collinear spins')

        else
          xdm_bool = .true.
          xdm_boolc = .false.
          xdm_mesh  = fdf_integer('XDM.Mesh',1)
          ! Default =     : 1, Siesta integration mesh
          ! Future options(?) : 2, Becke "Fuzzy cell" integration mesh over atoms using Franchini weights,
          ! see M. Franchini et al. J. Comput. Chem. (2013) 34, 1819, and
          ! A. D. Becke J. Chem. Phys. (1988) 88, 2547 **Not currently implemented**
          if (xdm_mesh .gt. 1 .or. xdm_mesh .le. 0) then 
            call die('XDM.Mesh has an unknown value. Check input.')
          end if
          coef_once = fdf_boolean('XDM.Once',.true.)
          coef_savd = .false.

          ! Read in XDM parameters from fdf input and write to output
          xdm_a1 = fdf_double('XDM.a1',default)
          xdm_a2 = fdf_double('XDM.a2',default)
          xdm_sc = fdf_double('XDM.sc',default)

          if (IOnode) then
            write(6,*)
            write(6,"(a,f8.4,a,f8.4,a)") "initxdm: a1=", xdm_a1,&
                  ", a2=", xdm_a2, " ang,",&
                  ", sc=", xdm_sc, "."
            write(6,*)
          end if

        end if

      else
        xdm_bool = .false.
        xdm_boolc = .false.
        coef_once = .true.
        coef_savd = .false.

      endif

!--------------------------------------------------------------------------- END
      END SUBROUTINE initxdm

      SUBROUTINE coefxdm( np, nspin, rhoatm, rhoscf, iaorb, iphorb,&
                          isa, dvol, xa, cell )

      use atmfuncs, only: chcore_sub, izofis, phiatm, rcut
      use atomlist, only: indxua, indxuo, no_s, Datm
      use mesh, only: dxa, nmuc, nsp, xdop, xdsp
      use meshphi
      use siesta_geom, only: na_u
#ifdef MPI
      use mpi
#endif

      IMPLICIT NONE

      integer,intent(in) :: np
      integer,intent(in) :: nspin
      integer,intent(in) :: iaorb(*)
      integer,intent(in) :: iphorb(*)
      integer,intent(in) :: isa(*)

      real(grid_p),intent(in) :: rhoatm(nsp,np)
      real(grid_p),intent(in) :: rhoscf(nsp,np,nspin)

      real(dp),intent(in) :: dvol
      real(dp),intent(in) :: xa(3,*)
      real(dp),intent(in) :: cell(3,3)

      ! Local variables
      integer :: atomic_number
      integer :: i, ia, iep, iix,&
                 io, iop, ip, iphi, is, isp, ispin,&
                 iu, iua, iuo, iui, iuj, ix
      integer :: j, kn
      integer :: nmucl(3)
      integer :: x1, y2, z3

      real(dp) :: apol(na_u), atvol(na_u)
      real(dp) :: Ci, ch
      real(dp) :: denom, dxsp(3)
      real(dp) :: freevol(na_u)
      real(dp) :: gradCi(3), grch(3), grrho
      real(dp) :: lap, ml(2:4,na_u)
      real(dp) :: phip
      real(dp) :: Qi, qh(na_u)
      real(dp) :: r, r2o, r2sp, rc(na_u,na_u), rho
      real(dp) :: tau, x(3), xroot

      real(grid_p),allocatable :: rhoatc(:,:) ! nsp, np

      real(dp),allocatable :: bi(:,:,:),& ! nsp, np, nspin
                              wei(:,:) ! nsp, np

#ifdef MPI
      real(dp) :: atvolS(na_u), apolS(na_u), freevolS(na_u),&
                  mlS(2:4,na_u), qhS(na_u)
      integer :: MPIerror
#endif

      ! Constants
      real(dp),parameter :: pi = 3.14159265358979323846_dp
      real(dp),parameter :: epsilon = 1.0D-19 ! Avoid division by zero
      real(grid_p),parameter :: rhoatm_tolerance = 1.0e-12_grid_p ! Avoid division by zero

      ! Free atomic polarizabilities from CRC handbook, 88th ed. (ang^3)
      real(dp),parameter :: afree(1:102) = (/ 0.6668_DP, 0.2051_DP,&
      24.3300_DP, 5.6000_DP, 3.0300_DP, 1.7600_DP, 1.1000_DP,&
      0.8020_DP, 0.5570_DP, 0.3956_DP, 24.1100_DP, 10.6000_DP,&
      6.8000_DP, 5.3800_DP, 3.6300_DP, 2.9000_DP, 2.1800_DP, 1.6411_DP,&
      43.4000_DP, 22.8000_DP, 17.8000_DP, 14.6000_DP, 12.4000_DP,&
      11.6000_DP, 9.4000_DP, 8.4000_DP, 7.5000_DP, 6.8000_DP,&
      6.2000_DP, 5.7500_DP, 8.1200_DP, 6.0700_DP, 4.3100_DP, 3.7700_DP,&
      3.0500_DP, 2.4844_DP, 47.3000_DP, 27.6000_DP, 22.7000_DP,&
      17.9000_DP, 15.7000_DP, 12.8000_DP, 11.4000_DP, 9.6000_DP,&
      8.6000_DP, 4.8000_DP, 7.2000_DP, 7.3600_DP, 10.2000_DP,&
      7.7000_DP, 6.6000_DP, 5.5000_DP, 5.3500_DP, 4.0440_DP,&
      59.4200_DP,39.7000_DP, 31.1000_DP, 29.6000_DP, 28.2000_DP,&
      31.4000_DP, 30.1000_DP, 28.8000_DP, 27.7000_DP, 23.5000_DP,&
      25.5000_DP,24.5000_DP, 23.6000_DP, 22.7000_DP, 21.8000_DP,&
      21.0000_DP,21.9000_DP, 16.2000_DP, 13.1000_DP, 11.1000_DP,&
      9.7000_DP,8.5000_DP, 7.6000_DP, 6.5000_DP, 5.8000_DP, 5.0200_DP,&
      7.6000_DP,6.8000_DP, 7.4000_DP, 6.8000_DP, 6.0000_DP, 5.3000_DP,&
      48.6000_DP,38.3000_DP, 32.1000_DP, 32.1000_DP, 25.4000_DP,&
      24.9000_DP,24.8000_DP, 24.5000_DP, 23.3000_DP, 23.0000_DP,&
      22.7000_DP,20.5000_DP, 19.7000_DP, 23.8000_DP, 18.2000_DP,&
      17.5000_DP /)

      ! Conversion factors
      character(len=80) :: scale
      real(dp) :: Dscale, Escale
        scale = fdf_string( 'MM.UnitsEnergy','eV' )
        Escale = fdf_convfac(scale,'Ry')
        scale = fdf_string( 'MM.UnitsDistance','Ang' )
        Dscale = fdf_convfac(scale,'Bohr')

!--------------------------------------------------------------------------- BEGIN

      ! Start timer
      call timer('coefXDM', 1)

      if (IOnode) then
        write(6,*)
        write(6,"(a)") "* XDM dispersion"
        write(6,"(a,f8.4)") "a1 =", xdm_a1
        write(6,"(a,f8.4)") "a2 (ang) =", xdm_a2
        write(6,"(a,f8.4)") "a2 (bohr) =", xdm_a2*Dscale
        write(6,"(a,f8.4)") "sc =", xdm_sc
        write(6,*)
      end if

      ! Compute pairwise dispersion coefficients of the XDM dispersion model

      ! Which mesh will be used to integrate the free and atomic
      ! volumes, as well as the  multipole moments?

!!!!!!
      if (xdm_mesh == 1) then
      !xdm_mesh = 1: Use Siesta integration mesh
!!!!!!
      ! Allocations and Initializations
      allocate ( bi(nsp,np,nspin), wei(nsp,np),&
                 rhoatc(nsp,np) )

      rhoatc(:,:) = 0._grid_p

      bi(:,:,:) = 0._dp
      freevol(:) = 0._dp
      atvol(:) = 0._dp
      qh(:) = 0._dp
      wei(:,:) = 0._dp
      ml(:,:) = 0._dp
#ifdef MPI
      mlS(:,:) = 0._dp
      freevolS(:) = 0._dp
      atvolS(:) = 0._dp
#endif

      ! Compute Becke-Roussel hole at each mesh point
      do ip = 1,np ! Loop over mesh points
        do isp = 1,nsp ! Loop over mesh sub-points
          do ispin = 1, nspin ! Loop over spins

            ! Assign the density, its gradient and laplacian, as well as the kinetic-
            ! energy density from the quantities computed in a previous call to
            ! rhoofd in dhscf
            if (rho_r(isp,ip,ispin).lt.1E-14_dp) cycle
            rho = rho_r(isp,ip,ispin)
            grrho = grrho_r(isp,ip,ispin)
            lap = laprho_r(isp,ip,ispin)
            tau = tau_r(isp,ip,ispin)

            if (nspin.eq.1) then
              rho = max(rho_r(isp,ip,ispin),1E-14_dp)/2._dp
              grrho = grrho/2._dp
              lap = lap/2._dp
              tau = abs(tau)/2._dp ! tau must be positive definite
            end if

            ! Solve for the b parameter
            call inversebr(rho,grrho,lap,tau,xroot)

            ! Compute Becke-Roussel hole at each mesh point
            bi(isp,ip,ispin) = xroot * (exp(-xroot)/&
                                  (8._dp*pi*rho + epsilon))**(1._dp/3._dp)

          end do ! Loop over spins
        end do ! Loop over mesh sub-points
      end do ! Loop over mesh points

      ! Perform Hirshfeld partitioning and compute atomic polarizabilities and volumes
      ! ^adapted from m_partial_charges.F

      ! Compute promolecular 'Harris' density
      call rhooda( no_s, np, Datm, rhoatc, iaorb, iphorb, isa,&
                   .true. ) ! and include core-corrections

      do ip = 1, np ! Loop over mesh points
        do kn = 1+endpht(ip-1), endpht(ip) ! Loop over non-zero orbitals
                                           ! crossing mesh point
          i = lstpht(kn)
          iu = indxuo(i)

          ! Generate phi value
          iphi = iphorb(i)
          ia = iaorb(i)
          iua = indxua(ia)
          is = isa(ia)
          r2o = rcut(is,iphi)**2
          iop = listp2(kn)

          do isp = 1, nsp ! Loop on mesh sub-points

            ! Ignore vanishing densities
            if (abs(rhoatm(isp,ip)).lt.rhoatm_tolerance) cycle

            dxsp(:) = xdop(:,iop) + xdsp(:,isp) - dxa(:,ia)
            r2sp = sum(dxsp(:)**2)
            if (r2sp.lt.r2o) then
              call phiatm(is,iphi,dxsp,phip,gradCi)
              call chcore_sub( is, dxsp, ch, grch )

              Ci = phip
              Qi = Datm(iu)*Ci*Ci
              freevol(iua) = freevol(iua) + r2sp**(3._dp/2._dp)&
                             * (Qi+ch) * dvol
              wei(isp,ip) = (Qi+ch) / rhoatc(isp,ip)
              do ispin = 1, nspin
                qh(iua) = qh(iua) + (Qi+ch)&
                          *(rho_r(isp,ip,ispin)+ch) / rhoatc(isp,ip)

                ! Ignore vanishing Hirshfeld weights
                if (wei(isp,ip).lt.1E-14_dp) cycle

                if (nspin.gt.1) then
                  if (rho_r(isp,ip,ispin).le.1E-14_dp) cycle
                  atvol(iua) = atvol(iua) + r2sp**(3._dp/2._dp)&
                              *0.5*(rho_r(isp,ip,1)+rho_r(isp,ip,2)+ch)&
                              *wei(isp,ip) *dvol
                else
                  if (rho_r(isp,ip,ispin).le.1E-14_dp) cycle
                  atvol(iua) = atvol(iua) + r2sp**(3._dp/2._dp)&
                               *(rho_r(isp,ip,ispin)+ch)&
                               *wei(isp,ip) *dvol
                end if

                ! Compute multipole moments for each atom
                do iix = 2,4 ! nth-order multipole moments
                  ml(iix,iua) = ml(iix,iua) + rho_r(isp,ip,ispin)&
                               *( (sqrt(r2sp))**(iix-1)&
                                  - (max( 0._dp, ((sqrt(r2sp)&
                                  - bi(isp,ip,ispin))**(iix-1)) )) )**2&
                               *wei(isp,ip) *dvol
                end do ! nth-order multipole moments
              end do ! Loop over spins
            end if
          end do ! Loop over mesh sub-points
        end do ! Loop over non-zero orbitals crossing mesh point
      end do ! Loop over mesh points

#ifdef MPI
      call MPI_Barrier(MPI_Comm_World,MPIerror)
      call MPI_AllReduce(freevol,freevolS,na_u,MPI_double_precision,&
                         MPI_sum,MPI_Comm_World,MPIerror)
      call MPI_AllReduce(atvol,atvolS,na_u,MPI_double_precision,&
                         MPI_sum,MPI_Comm_World,MPIerror)
      call MPI_AllReduce(ml,mlS,3*na_u,MPI_double_precision,&
                         MPI_sum,MPI_Comm_World,MPIerror)
      freevol(1:na_u) = freevolS(1:na_u)
      atvol(1:na_u) = atvolS(1:na_u)
      ml(2:4,1:na_u) = mlS(2:4,1:na_u)
#endif

      ! Deallocations
      deallocate( rho_r, grrho_r, laprho_r, tau_r, wei, bi )

!!!!!!
      end if ! xdm_mesh
!!!!!!

      ! Compute atomic polarizabilities
      do iu = 1, na_u ! Loop over atoms
        is = isa(iu)
        atomic_number = izofis(is)
        apol(iu) = min(atvol(iu) / freevol(iu),1._dp)&
                   * afree(atomic_number)*Dscale**3
      end do ! Loop over atoms
 
      ! Proceed with dispersion coefficients
      if (IOnode) then
        write(6,*)
        write(6,"(a)") "+ Volumes and moments"
        write(6,"(a)") "# All results in atomic units (Hartree,bohr)"
        write(6,"(a)") "# i        V             Vfree           M1&
      &             M2             M3"
        do iui = 1, na_u
          write(6,"(i3,1p,5(1x,e14.6))") iui, atvol(iui), freevol(iui),&
                ml(2:4,iui)
        end do
        write(6,*)
      end if

      ! Compute dispersion coefficients, critical and vdW radii
      do iui = 1, na_u ! Loop over first atom
        do iuj = 1, iui ! Loop over second atom

          denom = ml(2,iui)*apol(iuj) + ml(2,iuj)*apol(iui)

          cx(iui,iuj,2) = apol(iui)*apol(iuj)*ml(2,iui)*ml(2,iuj)&
                          /(denom + epsilon)
          cx(iuj,iui,2) = cx(iui,iuj,2)

          cx(iui,iuj,3) = 3._dp/2._dp*(apol(iui)*apol(iuj)&
                          *(ml(2,iui)*ml(3,iuj) + ml(3,iui)*ml(2,iuj))&
                          /(denom  + epsilon))
          cx(iuj,iui,3) = cx(iui,iuj,3)

          cx(iui,iuj,4) = 2._dp*(apol(iui)*apol(iuj)*(ml(2,iui)&
                          *ml(4,iuj) + ml(4,iui)* ml(2,iuj))/(denom + epsilon))&
                          + 21._dp/5._dp*(apol(iui)*apol(iuj)&
                          *ml(3,iui)*ml(3,iuj)/(denom + epsilon))
          cx(iuj,iui,4) = cx(iui,iuj,4)

          rc(iui,iuj) = 1._dp/3._dp * ((cx(iui,iuj,3)/(cx(iui,iuj,2) + epsilon))**0.5_dp&
                        + (cx(iui,iuj,4)/(cx(iui,iuj,2) + epsilon))**0.25_dp&
                        + (cx(iui,iuj,4)/(cx(iui,iuj,3) + epsilon))**0.5_dp)
          rvdw(iui,iuj) = xdm_a1*rc(iui,iuj) + xdm_a2*Dscale

        end do ! Loop over second atom
      end do ! Loop over first atom

      if (IOnode) then
        write(6,"(a)") "+ Dispersion coefficients"
        write(6,"(a)") "# All results in atomic units (Hartree,bohr)"
        write(6,"(a)") "# i   j      C6             C8             C10 &
      &           Rc            Rvdw"

        do iui = 1, na_u ! Loop over first atom
          do iuj = 1, iui ! Loop over second atom
            write(6,"(i3,1x,i3,1p,5(1x,e14.6))") iui, iuj,&
                  cx(iui,iuj,2), cx(iui,iuj,3), cx(iui,iuj,4),&
                  rc(iui,iuj), rvdw(iui,iuj)
          end do ! Loop over second atom
        end do ! Loop over first atom
      end if

      coef_savd = .true.

      ! Stop timer
      call timer('coefXDM', 2)

!--------------------------------------------------------------------------- END
      END SUBROUTINE coefxdm

      SUBROUTINE enfrxdm(na, xa, isa, cell, exdmt, ifa, fa, istr,&
                         stress)

      use units, only: kBar

      IMPLICIT NONE

      integer,intent(in) :: na
      real(dp),intent(in) :: xa(3,*)
      integer,intent(in) :: isa(*)
      real(dp), intent(in) :: cell(3,3)
      integer,intent(in) :: ifa ! compute forces if > 0
      integer,intent(in) :: istr ! compute stess if > 0

      real(dp),intent(out) :: exdmt

      real(dp),intent(inout) :: fa(3,*)
      real(dp),intent(inout) :: stress(3,3)

      ! Local variables
      integer :: i, ii, idir, imid, ix, iix
      integer :: j, jj, jdir, jmid, jx
      integer :: k, kk, kdir, kmid
      integer :: nvec, nvecj0, nveck0

      logical :: lallfound1, lallfound2, lallfound3

      real(dp) :: a, alpha
      real(dp) :: b, beta
      real(dp) :: c
      real(dp) :: dbcount, denom
      real(dp) :: exdm(2:4), exdmc(2:4)
      real(dp) :: fdamp(na,na,2:4), faxdm(1:3,na), frxdm(2:4), frxdmt
      real(dp) :: gamma
      real(dp) :: proj1, proj2, proj3
      real(dp) :: rcx1, rcx2, rcx3, rx, rxi, rxj
      real(dp) :: rcy1, rcy2, rcy3, ry, ryi, ryj 
      real(dp) :: rcz1, rcz2, rcz3, rz, rzi, rzj
      real(dp) :: recipa, recipb, recipc
      real(dp) :: r2, r2i, r2j, r2k, rnorm, rvol
      real(dp) :: stxdm(1:3,1:3)
      real(dp) :: vol
      real(dp),external :: volcel
      real(dp) :: x, xdm_cutoff2, y, z

      ! Constants
      real(dp),parameter :: epsilon = 1.0D-19 ! Avoid division by zero

      ! Conversion factors
      real(dp),parameter :: ha2ry = 2._dp
      real(dp),parameter :: ry2gpa = 14710.508_dp

      character(len=80) :: scale
      real(dp) :: Dscale, Escale
        scale = fdf_string ( 'MM.UnitsEnergy','eV' )
        Escale = fdf_convfac(scale,'Ry')
        scale = fdf_string( 'MM.UnitsDistance','Ang' )
        Dscale = fdf_convfac(scale,'Bohr')

      ! Get potential cutoff
      xdm_cutoff = fdf_physical('XDM.Cutoff',40.0d0,'Bohr')

!--------------------------------------------------------------------------- BEGIN

      ! Start timer
      call timer('enfrXDM', 1)

      ! Quick return if XDM is not required
      if (.not.xdm_bool) then
        call timer('enfrXDM',2)
        return
      end if

      ! Compute XDM dispersion energy

      ! Nullify and allocate dynamic arrays.
      if (.not.coef_savd) then
        allocate ( cx(na,na,2:4), rvdw(na,na) )
        cx(:,:,:) = 0._dp
        rvdw(:,:) = 0._dp
      end if

      ! The following has been adapted from molecularmechanics.F90
      ! Initialize relevant quantities
      fdamp(:,:,:) = 0.0_dp

      exdm(:) = 0.0_dp
      exdmc(:) = 0.0_dp
      exdmt = 0.0_dp

      frxdm(:) = 0.0_dp
      frxdmt = 0.0_dp
      faxdm(:,:) = 0.0_dp

      stxdm(:,:) = 0.0_dp

      ! Find number of cell images required
      xdm_cutoff2 = xdm_cutoff**2
      call uncell(cell,a,b,c,alpha,beta,gamma,1.0_dp)
      recipa = 1.0_dp/a
      recipb = 1.0_dp/b
      recipc = 1.0_dp/c

      ! Find volume if required
      if (istr.ne.0) then
        vol = volcel(cell)
        rvol = 1.0_dp/vol
      end if

      do i = 1, na ! Loop over first atom
        do j = 1, i ! Loop over second atom

          ! Avoid double-counting
          if (i.eq.j) then
            dbcount = 0.5_dp
          else
            dbcount = 1.0_dp
          endif

          ! Find image of j nearest to i
          x = xa(1,j) - xa(1,i)
          y = xa(2,j) - xa(2,i)
          z = xa(3,j) - xa(3,i)

          ! Find projection of cell vector 3 on to i - j vector
          rnorm = x**2 + y**2 + z**2
          if (rnorm .gt. 1.0d-12) rnorm = 1.0_dp/sqrt(rnorm)
          proj3 = rnorm*recipc*(x*cell(1,3) + y*cell(2,3) + z*cell(3,3))
          kmid = nint(proj3)
          x = x - kmid*cell(1,3)
          y = y - kmid*cell(2,3)
          z = z - kmid*cell(3,3)

          ! Find projection of cell vector 2 on to i - j vector
          rnorm = x**2 + y**2 + z**2
          if (rnorm .gt. 1.0d-12) rnorm = 1.0_dp/sqrt(rnorm)
          proj2 = rnorm*recipb*(x*cell(1,2) + y*cell(2,2) + z*cell(3,2))
          jmid = nint(proj2)
          x = x - jmid*cell(1,2)
          y = y - jmid*cell(2,2)
          z = z - jmid*cell(3,2)

          ! Find projection of cell vector 1 on to i - j vector
          rnorm = x**2 + y**2 + z**2
          if (rnorm .gt. 1.0d-12) rnorm = 1.0_dp/sqrt(rnorm)
          proj1 = rnorm*recipa*(x*cell(1,1) + y*cell(2,1) + z*cell(3,1))
          imid = nint(proj1)
          x = x - imid*cell(1,1)
          y = y - imid*cell(2,1)
          z = z - imid*cell(3,1)

          ! Initialize counter for number of valid vectors
          nvec = 0

          ! Outer loop over first cell vector direction
          do idir = 1,-1,-2

            ! Reinitialize distance squared
            r2i = 10000.0_dp*xdm_cutoff2

            ! Loop over first cell vector
            lallfound1 = .false.
            if (idir.eq.1) then
              ii = 0
            else
              ii = -1
            endif

            ! Set initial coordinate vector
            rxi = x + dble(ii)*cell(1,1)
            ryi = y + dble(ii)*cell(2,1)
            rzi = z + dble(ii)*cell(3,1)

            ! Set increment vector
            rcx1 = dble(idir)*cell(1,1)
            rcy1 = dble(idir)*cell(2,1)
            rcz1 = dble(idir)*cell(3,1)

            do while (.not.lallfound1)

              ! Save number of vectors before search over second direction
              nvecj0 = nvec

              ! Outer loop over second cell vector direction
              do jdir = 1,-1,-2

                ! Reinitialize saved distance squared
                r2j = 10000.0_dp*xdm_cutoff2

                ! Loop over first cell vector
                lallfound2 = .false.
                if (jdir.eq.1) then
                  jj = 0
                else
                  jj = -1
                endif

                ! Set initial coordinate vector
                rxj = rxi + dble(jj)*cell(1,2)
                ryj = ryi + dble(jj)*cell(2,2)
                rzj = rzi + dble(jj)*cell(3,2)

                ! Set increment vector
                rcx2 = dble(jdir)*cell(1,2)
                rcy2 = dble(jdir)*cell(2,2)
                rcz2 = dble(jdir)*cell(3,2)

                do while (.not.lallfound2)

                  ! Save number of vectors before search over second direction
                  nveck0 = nvec

                  ! Outer loop over second cell vector direction
                    do kdir = 1,-1,-2

                    ! Reinitialize saved distance squared
                    r2k = 10000.0_dp*xdm_cutoff2

                    ! Loop over third cell vector
                    lallfound3 = .false.
                    if (kdir.eq.1) then
                      kk = 0
                    else
                      kk = -1
                    endif

                    ! Set initial coordinate vector
                    rx = rxj + dble(kk)*cell(1,3)
                    ry = ryj + dble(kk)*cell(2,3)
                    rz = rzj + dble(kk)*cell(3,3)

                    ! Set increment vector
                    rcx3 = dble(kdir)*cell(1,3)
                    rcy3 = dble(kdir)*cell(2,3)
                    rcz3 = dble(kdir)*cell(3,3)

                    do while (.not.lallfound3)

                      ! Calculate square of distance
                      r2 = rx**2 + ry**2 + rz**2

                      ! Check distance squared against cutoff squared
                      if (r2 .le. xdm_cutoff2) then
                        if (r2 .gt. 1.0d-10) then

                        ! Valid distance, so increment counter
                        nvec = nvec + 1

                        ! Evaluate dispersion energy contributions,
                        ! forces, and stresses for this valid distance
                        do ix = 2, 4
                          iix = 2*ix + 2
                          fdamp(i,j,ix) = 1/(rvdw(i,j)**iix&
                                          + r2**(iix/2) + epsilon)

                          exdm(ix) = -xdm_sc*cx(i,j,ix)*fdamp(i,j,ix)&
                                     *dbcount *ha2ry
                          exdmc(ix) = exdmc(ix) + exdm(ix)
                          exdmt = exdmt + exdm(ix)

                          frxdm(ix) = iix*r2**ix&
                                      *fdamp(i,j,ix)*exdm(ix)
                          frxdmt = frxdmt + frxdm(ix)
                        end do

                        if (ifa .ne. 0) then
                          ! Gradients
                          faxdm(1,i) = faxdm(1,i) + frxdmt*rx
                          faxdm(2,i) = faxdm(2,i) + frxdmt*ry
                          faxdm(3,i) = faxdm(3,i) + frxdmt*rz
                          faxdm(1,j) = faxdm(1,j) - frxdmt*rx
                          faxdm(2,j) = faxdm(2,j) - frxdmt*ry
                          faxdm(3,j) = faxdm(3,j) - frxdmt*rz

                          fa(1,i) = fa(1,i) + frxdmt*rx
                          fa(2,i) = fa(2,i) + frxdmt*ry
                          fa(3,i) = fa(3,i) + frxdmt*rz
                          fa(1,j) = fa(1,j) - frxdmt*rx
                          fa(2,j) = fa(2,j) - frxdmt*ry
                          fa(3,j) = fa(3,j) - frxdmt*rz
                        end if ! ifa.ne.0

                        if (istr .ne. 0) then
                          ! Stress
                          frxdmt = frxdmt*rvol
                          stxdm(1,1) = stxdm(1,1) - frxdmt*rx*rx
                          stxdm(2,1) = stxdm(2,1) - frxdmt*ry*rx
                          stxdm(3,1) = stxdm(3,1) - frxdmt*rz*rx
                          stxdm(1,2) = stxdm(1,2) - frxdmt*rx*ry
                          stxdm(2,2) = stxdm(2,2) - frxdmt*ry*ry
                          stxdm(3,2) = stxdm(3,2) - frxdmt*rz*ry
                          stxdm(1,3) = stxdm(1,3) - frxdmt*rx*rz
                          stxdm(2,3) = stxdm(2,3) - frxdmt*ry*rz
                          stxdm(3,3) = stxdm(3,3) - frxdmt*rz*rz
                        endif ! istr.ne.0

                        frxdmt = 0.0_dp

                      end if ! .gt. 1.0d-10
                    endif ! .le. xdm_cutoff2

                    ! Increment by third vector
                    kk = kk + kdir
                    rx = rx + rcx3
                    ry = ry + rcy3
                    rz = rz + rcz3

                    ! Check to see if this direction is complete
                    lallfound3 = (r2.gt.r2k .and. r2.gt.xdm_cutoff2)
                    r2k = r2
                    end do ! End loop lallfound3
                  end do ! End loop kdir

                ! Increment by second vector
                jj = jj + jdir
                rxj = rxj + rcx2
                ryj = ryj + rcy2
                rzj = rzj + rcz2

              ! Check to see if this distance is complete
              lallfound2 = (r2.gt.r2j .and. r2.gt.xdm_cutoff2 .and.&
                                 nvec.eq.nveck0)
              r2j = r2
              end do ! End loop lallfound2
            end do ! End loop jdir

            ! Increment by first vector
            ii = ii + idir
            rxi = rxi + rcx1
            ryi = ryi + rcy1
            rzi = rzi + rcz1

            ! Check to see if this direction is complete
            lallfound1 = (r2.gt.r2i .and. r2.gt.xdm_cutoff2 .and.&
                         nvec.eq.nvecj0)
            r2i = r2
            end do ! End loop allfound1
          end do ! End loop idir
        end do ! Loop over second atom
      end do ! Loop over first atom

      if (istr.ne.0) then
        stress = stress + stxdm
      end if

      if (IOnode) then
        write(6,*)
        write(6,"(a)") "+ van der Waals energies, forces and stresses&
                        &(Ry/bohr)"
        write(6,"(a, e20.12)") "Evdw(total,Ry) = ", exdmt
        write(6,"(a, e20.12)") "Evdw(C6,Ry) = ", exdmc(2)
        write(6,"(a, e20.12)") "Evdw(C8,Ry) = ", exdmc(3)
        write(6,"(a, e20.12)") "Evdw(C10,Ry) = ", exdmc(4)
        write(6,*)
        if (ifa.ne.0) then
          do i = 1, na
          write (6,'("  Fvdw (",i3.3,",Ry/bohr) = ",1p,3(e20.12,1x))')&
                 i, faxdm(:,i)
          end do
          write(6,*)
        end if
        if (istr.ne.0) then
        write (6,'("  sigma_vdw (Ry/bohr**3) = ",1p,3(e20.12,1x)," ")')&
               stxdm(1,:)
        write (6,'("                           ",1p,3(e20.12,1x)," ")')&
               stxdm(2,:)
        write (6,'("                           ",1p,3(e20.12,1x)," ")')&
               stxdm(3,:)
        write (6,'("  sigma_vdw (GPa) = ",1p,3(e20.12,1x)," ")')&
               stxdm(1,:)*ry2gpa
        write (6,'("                    ",1p,3(e20.12,1x)," ")')&
               stxdm(2,:)*ry2gpa
        write (6,'("                    ",1p,3(e20.12,1x)," ")')&
               stxdm(3,:)*ry2gpa
        write(6,'(/,a,6f12.2)')  '  stress-tensor-voigt_vdw (kbar):',  &
             (stxdm(jx,jx)/kbar,jx=1,3),            &
              stxdm(1,2)/kbar,                      &
              stxdm(2,3)/kbar,                      &
              stxdm(1,3)/kbar
        end if
      end if

      ! Stop timer
      call timer('enfrXDM', 2)

!--------------------------------------------------------------------------- END
      END SUBROUTINE enfrxdm

      ! Solve for the b parameter in the Becke-Roussel exchange hole model
      SUBROUTINE inversebr(rho,grrho,lap,tau,xroot)

      IMPLICIT NONE

      real(dp),intent(in) :: rho, grrho, lap, tau
      real(dp),intent(out) :: xroot

      ! Local variables
      real(dp) :: wz, ds, qs, rhs, xshift, xold, expx, ffx, fx, gx

      ! Constants
      real(dp),parameter :: pi = 3.14159265358979323846_dp
      real(dp),parameter :: epsilon = 1.0D-19 ! Avoid division by zero
      
!--------------------------------------------------------------------------- BEGIN

      ! Compute relevant quantities from the density,its  gradient and
      ! laplacian, and from the kinetic-energy density
      wz = grrho**2 / (rho + epsilon) ! Compute Weizsaecker term
      ds = tau - 0.25_dp * wz ! Dsigma
      qs = 1._dp/6._dp * (lap - 2._dp * ds) ! Qsigma

      ! Solve for x in non-linear equation
      rhs = 2._dp/3._dp * pi**(2._dp/3._dp)&
            * (rho)**(5._dp/3._dp) / (qs + epsilon)

      ! Newton seed
      if (rhs .gt. 0._dp) then
        xroot = 3._dp
        xshift = 1._dp
        do while ((xroot * exp(-2._dp*xroot/3._dp))/&
                  (xroot-2._dp) .lt. rhs)
          xshift = xshift * 0.1_dp
          xroot = 2._dp + xshift
        end do
      else
        xroot = 1._dp
        xshift = 1._dp
        do while ((xroot * exp(-2._dp*xroot/3._dp))/&
                  (xroot-2._dp) .gt. rhs)
          xshift = xshift * 0.1_dp
          xroot = 2._dp - xshift
        end do
      end if

      ! Newton's method
      xold = 2._dp
      do while (abs(xroot - xold) .gt. 1e-10_dp)
        xold = xroot
        expx = exp(-2._dp * xroot / 3._dp)
        gx = (xroot * expx) / (xroot - 2._dp)
        fx = gx - rhs
        ffx = gx * (1._dp / xroot - 2._dp/3._dp - 1._dp / (xroot&
              - 2._dp))
        xroot = xroot - fx / ffx
      end do

!--------------------------------------------------------------------------- END
      END SUBROUTINE inversebr

      END MODULE m_xdm

