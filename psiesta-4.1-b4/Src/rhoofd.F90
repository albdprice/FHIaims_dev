! 
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!
module m_rhoofd

  implicit none

  private
  public :: rhoofd

contains


subroutine rhoofd( no, np, maxnd, numd, listdptr, listd, nspin, &
     Dscf, rhoscf, nuo, nuotot, iaorb, iphorb, isa )
! ********************************************************************
! Finds the SCF density at the mesh points from the density matrix.
! Written by P.Ordejon and J.M.Soler. May'95.
! Re-ordered so that mesh is the outer loop and the orbitals are
! handled as lower-half triangular. J.D.Gale and J.M.Soler, Feb'99
! Version of rhoofd that optionally uses a direct algorithm to save 
! memory. Modified by J.D.Gale, November'99
! *********************** InpUT **************************************
! integer no              : Number of basis orbitals
! integer np              : Number of mesh points
! integer maxnd           : First dimension of listD and Dscf, and
!                           maximum number of nonzero elements in
!                           any row of Dscf
! integer numd(nuo)       : Number of nonzero elemts in each row of Dscf
! integer listdptr(nuo)   : Pointer to start of rows in listd
! integer listd(maxnd)    : List of nonzero elements in each row of Dscf
! integer nspin           : Number of spin components
! real*8  Dscf(maxnd)     : Rows of Dscf that are non-zero 
! integer nuo             : Number of orbitals in unit cell locally
! integer nuotot          : Number of orbitals in unit cell in total
! integer iaorb(*)        : Pointer to atom to which orbital belongs
! integer iphorb(*)       : Orbital index within each atom
! integer isa(*)          : Species index of all atoms
! *********************** OUTPUT **************************************
! real    rhoscf(nsp,np)  : SCF density at mesh points
! *********************************************************************

!  Modules

  use precision,     only: dp, grid_p
  use atmfuncs,      only: rcut, all_phi
  use atm_types,     only: nsmax=>nspecies
  use atomlist,      only: indxuo
  use m_spin,        only: SpOrb
  use listsc_module, only: LISTSC
  use mesh,          only: nsp, dxa, xdop, xdsp, meshLim
  use meshdscf,      only: matrixOtoM
  use meshdscf,      only: nrowsDscfL, listdl, listdlptr, NeedDscfL, &
       numdl, DscfL
  use meshphi,       only: DirectPhi, endpht, lstpht, listp2, phi
  use parallel,      only: Nodes, node
  use sys,           only: die
  use alloc,         only: re_alloc, de_alloc
  use parallelsubs,  only: GlobalToLocalOrb
  use m_xdm,         only: xdm_bool, xdm_boolc, xdm_mesh, &
                           rho_r, grrho_r, laprho_r, tau_r                    
! ^LML                                                                            
#ifdef MPI
  use mpi_siesta
#endif
  implicit none

!     Argument types and dimensions
  integer, intent(in) :: no, np, nuo, maxnd, numd(nuo), nspin, &
       nuotot, iaorb(*), iphorb(*), &
       isa(*), listdptr(nuo), listd(maxnd)
  real(dp), intent(in) :: Dscf(:,:)
  real(grid_p), intent(out) :: rhoscf(nsp,np,nspin)
  external                  :: memory, timer
!     Internal variables and arrays
  integer, parameter :: minloc = 1000 ! Min buffer size for local copy of Dscf
  integer, parameter :: maxoa  = 100  ! Max # of orbitals per atom
  logical :: ParallelLocal
  integer :: i, ia, ic, ii, ik, ijl, il, imp, ind
! ^LML
  integer :: ispin, io, iop, ip, iphi, is
  integer :: isp, iu, iul, j, jc, jk, last, lasta
! ^LML
  integer :: lastop, maxloc, maxloc2, triang, nc
  integer :: maxndl, nphiloc, lenx, leny, lenxy, lenz

  ! Total hamiltonian size
  integer :: h_spin_dim
  
  real(dp) :: r2sp, dxsp(3)
  real(dp) :: Dij, Cij(nsp), &
              E1ij(3,nsp), E2ij(3,nsp), Gij(3,nsp), &
              F1ij(3,3,nsp), F2ij(3,3,nsp), &
              Kij(3,3,nsp), grrho(3,nsp,np,nspin), &
              Ggrrho(3,3,nsp,np,nspin), &
              tau(3,nsp,np,nspin)
! ^LML
  integer, pointer :: ilc(:), ilocal(:), iorb(:)
  real(dp), pointer :: r2cut(:), Clocal(:,:), Dlocal(:,:), phia(:,:), &
                       Glocal(:,:,:), grphia(:,:,:), &
                       Flocal(:,:,:,:), Ggrphia(:,:,:,:)
! ^LML
#ifdef _TRACE_
  integer :: MPIerror
#endif

#ifdef DEBUG
  call write_debug( '    PRE rhoofd' )
#endif
#ifdef _TRACE_
  call MPI_Barrier( MPI_Comm_World, MPIerror )
  call MPItrace_event( 1000, 1 )
#endif
!     Start time counter
  call timer('rhoofd',1)

  ! Get spin-size
  h_spin_dim = size(Dscf, 2)

!     Set algorithm logical
  ParallelLocal = (Nodes > 1)

  if (ParallelLocal) then
     if (nrowsDscfL > 0) then
        maxndl = listdlptr(nrowsDscfL) + numdl(nrowsDscfL)
     else
        maxndl = 1
     end if
     nullify(DscfL)
     call re_alloc( DscfL, 1, maxndl, 1, h_spin_dim, 'DscfL', 'rhoofd' )
!       Redistribute Dscf to DscfL form
     call matrixOtoM( maxnd, numd, listdptr, maxndl, nuo, &
          h_spin_dim, Dscf, DscfL )
  end if

!     Find atomic cutoff radii
  nullify(r2cut)
  call re_alloc( r2cut, 1, nsmax, 'r2cut', 'rhoofd' )
  r2cut = 0.0_dp
  do i = 1,nuotot
     ia = iaorb(i)
     is = isa(ia)
     io = iphorb(i)
     r2cut(is) = max( r2cut(is), rcut(is,io)**2 )
  end do

! Find size of buffers to store partial copies of Dscf and C
  maxloc2 = maxval(endpht(1:np)-endpht(0:np-1))
  maxloc = maxloc2 + minloc
  maxloc = min( maxloc, no )
  triang = (maxloc+1)*(maxloc+2)/2

  lenx  = meshLim(2,1) - meshLim(1,1) + 1
  leny  = meshLim(2,2) - meshLim(1,2) + 1
  lenz  = meshLim(2,3) - meshLim(1,3) + 1
  lenxy = lenx*leny

!$OMP parallel default(shared), &
!$OMP&shared(rhoscf), &
!$OMP&private(ilocal,ilc,iorb,Dlocal,Clocal,phia), &
!$OMP&private(ip,nc,ic,imp,i,il,last,j,iu,iul,ii,ind,io), &
!$OMP&private(ijl,ispin,lasta,lastop,ia,is,iop,isp,iphi), &
!$OMP&private(jc,nphiloc,dxsp,r2sp)

! Allocate local memory
  nullify ( ilocal, ilc, iorb, Dlocal, Clocal, phia )
  if (xdm_boolc .and. xdm_mesh == 1) then
    nullify ( Glocal, Flocal, grphia, Ggrphia )
  end if
! ^LML
!$OMP critical
  allocate( ilocal(no), ilc(maxloc2), iorb(maxloc) )
  allocate( Dlocal(triang,nspin), Clocal(nsp,maxloc2) )
  if (xdm_boolc .and. xdm_mesh == 1) then
    allocate( Glocal(3,nsp,maxloc2), Flocal(3,3,nsp,maxloc2), &
              rho_r(nsp,np,nspin), grrho_r(nsp,np,nspin), &
              laprho_r(nsp,np,nspin), tau_r(nsp,np,nspin) )
    allocate( phia(maxoa,nsp), grphia(3,maxoa,nsp), Ggrphia(3,3,maxoa,nsp) )
  else
    if ( DirectPhi ) allocate( phia(maxoa,nsp) )
  end if
! ^LML
!$OMP end critical

! Initializations
  Dlocal(:,:) = 0.0_dp
  ilocal(:)   = 0
  iorb(:)     = 0

  last = 0

!$OMP do
  do ip = 1,np

!    Initializations
     rhoscf(:,ip,:) = 0.0_grid_p
     if (xdm_boolc .and. xdm_mesh == 1) then
       rho_r(:,ip,:) = 0.0_grid_p
       grrho(:,:,ip,:) = 0.0_grid_p
       grrho_r(:,ip,:) = 0.0_grid_p
       Ggrrho(:,:,:,ip,:) = 0.0_grid_p
       laprho_r(:,ip,:) = 0.0_grid_p
       tau(:,:,ip,:) = 0.0_grid_p
       tau_r(:,ip,:) = 0.0_grid_p
     end if
! ^LML

!    Find number of nonzero orbitals at this point
     nc = endpht(ip) - endpht(ip-1)
!       iorb(il)>0 means that row il of Dlocal must not be overwritten
!       iorb(il)=0 means that row il of Dlocal is empty
!       iorb(il)<0 means that row il of Dlocal contains a valid row of 
!             Dscf, but which is not required at this point
     do ic = 1,nc
        imp = endpht(ip-1) + ic
        i = lstpht(imp)
        il = ilocal(i)
        if (il > 0) iorb(il) = i
     end do

!    Look for required rows of Dscf not yet stored in Dlocal
     do ic = 1,nc
        imp = endpht(ip-1) + ic
        i = lstpht(imp)
        if (ilocal(i) == 0) then
!          Look for an available row in Dlocal
           do il = 1,maxloc
!             last runs circularly over rows of Dlocal
              last = last + 1
              if (last > maxloc) last = 1
              if (iorb(last) <= 0) goto 10
           end do
           call die('rhoofd: no slot available in Dlocal')
10         continue
!          Copy row i of Dscf into row last of Dlocal
           j = abs(iorb(last))
           if (j /= 0) ilocal(j) = 0
           ilocal(i)  = last
           iorb(last) = i
           il = last
           iu = indxuo(i)
           if ( ParallelLocal ) then
              iul = NeedDscfL(iu)
              if ( i == iu ) then
                 do ii = 1, numdl(iul)
                    ind = listdlptr(iul) + ii
                    j   = listdl(ind)
                    ijl = idx_ijl(il,ilocal(j))
                    if ( SpOrb ) then
                       Dlocal(ijl,1) = DscfL(ind,1)
                       Dlocal(ijl,2) = DscfL(ind,2)
                       Dlocal(ijl,3) = 0.5*(DscfL(ind,3)+DscfL(ind,7))
                       Dlocal(ijl,4) = 0.5*(DscfL(ind,4)+DscfL(ind,8))
                    else
                       Dlocal(ijl,:) = DscfL(ind,:)
                    end if
                 end do
              else
                 do ii = 1, numdl(iul)
                    ind = listdlptr(iul)+ii
                    j   = LISTSC( i, iu, listdl(ind) )
                    ijl = idx_ijl(il,ilocal(j))
                    if ( SpOrb ) then
                       Dlocal(ijl,1) = DscfL(ind,1)
                       Dlocal(ijl,2) = DscfL(ind,2)
                       Dlocal(ijl,3) = 0.5*(DscfL(ind,3)+DscfL(ind,7))
                       Dlocal(ijl,4) = 0.5*(DscfL(ind,4)+DscfL(ind,8))
                    else
                       Dlocal(ijl,:) = DscfL(ind,:)
                    end if
                 end do
              end if
           else
              call GlobalToLocalOrb( iu, Node, Nodes, iul )
              if ( i == iu ) then
                 do ii = 1, numd(iul)
                    ind = listdptr(iul)+ii
                    j   = listd(ind)
                    ijl = idx_ijl(il,ilocal(j))
                    if ( SpOrb ) then
                       Dlocal(ijl,1) = Dscf(ind,1)
                       Dlocal(ijl,2) = Dscf(ind,2)
                       Dlocal(ijl,3) = 0.5*(Dscf(ind,3)+Dscf(ind,7))
                       Dlocal(ijl,4) = 0.5*(Dscf(ind,4)+Dscf(ind,8))
                    else
                       Dlocal(ijl,:) = Dscf(ind,:)
                    end if
                 end do
              else
                 do ii = 1, numd(iul)
                    ind = listdptr(iul)+ii
                    j   = LISTSC( i, iu, listd(ind) )
                    ijl = idx_ijl(il,ilocal(j))
                    if ( SpOrb ) then
                       Dlocal(ijl,1) = Dscf(ind,1)
                       Dlocal(ijl,2) = Dscf(ind,2)
                       Dlocal(ijl,3) = 0.5*(Dscf(ind,3)+Dscf(ind,7))
                       Dlocal(ijl,4) = 0.5*(Dscf(ind,4)+Dscf(ind,8))
                    else
                       Dlocal(ijl,:) = Dscf(ind,:)
                    end if
                 end do
              end if
           end if
        end if
     end do

!    Check algorithm
     if ( DirectPhi ) then
        lasta = 0
        lastop = 0
        do ic = 1,nc
           imp = endpht(ip-1) + ic
           i   = lstpht(imp)
           il  = ilocal(i)
           ia  = iaorb(i)
           iop = listp2(imp)
           ilc(ic) = il
           
!          Generate or retrieve phi values
           if ( ia /= lasta .or. iop /= lastop ) then
              lasta  = ia
              lastop = iop
              is = isa(ia)
              do isp = 1 , nsp
                 dxsp(:) = xdsp(:,isp) + xdop(:,iop) - dxa(:,ia)
                 r2sp = sum(dxsp**2)
                 if ( r2sp < r2cut(is) ) then
!$OMP critical
                   if (xdm_boolc .and. xdm_mesh == 1) then
                     call all_phi( is, +1, dxsp, nphiloc, phia(:,isp), &
                                   grphia(:,:,isp), Ggrphia(:,:,:,isp) )
                   else
                     call all_phi( is, +1, dxsp, nphiloc, phia(:,isp) )
                   end if
! ^LML
!$OMP end critical
                 else
                    phia(:,isp) = 0.0_dp
                  if (xdm_boolc .and. xdm_mesh == 1) then
                    grphia(:,:,isp) = 0.0_dp
                    Ggrphia(:,:,:,isp) = 0.0_dp
                  end if
! ^LML
                 end if
              end do
           end if
           iphi = iphorb(i)

!          Retrieve phi values
           Clocal(:,ic) = dsqrt(2._dp) * phia(iphi,:)
           if (xdm_boolc .and. xdm_mesh == 1) then
             Glocal(:,:,ic) = grphia(:,iphi,:)
             Flocal(:,:,:,ic) = Ggrphia(:,:,iphi,:)
           end if
! ^LML

!          Loop on second orbital of mesh point
           do jc = 1, ic - 1
              ijl = idx_ijl(il,ilc(jc))
              
!             Loop over sub-points
              do isp = 1,nsp
                 if (xdm_boolc .and. xdm_mesh == 1) then
                   Gij(:,isp) = Glocal(:,isp,ic) * Glocal(:,isp,jc)
                   E1ij(:,isp) = Clocal(isp,ic) / dsqrt(2._dp) * Glocal(:,isp,jc)
                   E2ij(:,isp) = Glocal(:,isp,ic) * Clocal(isp,jc) / dsqrt(2._dp)
                   F1ij(:,:,isp) = Clocal(isp,ic) / dsqrt(2._dp) * Flocal(:,:,isp,jc)
                   F2ij(:,:,isp) = Flocal(:,:,isp,ic) * Clocal(isp,jc) / dsqrt(2._dp)
                   do ik = 1,3
                     Kij(ik,ik,isp) = Glocal(ik,isp,ic) * Glocal(ik,isp,jc)
                   end do
                 end if
! ^LML
                 do ispin = 1,nspin
                    rhoscf(isp,ip,ispin) = rhoscf(isp,ip,ispin) + &
                         Dlocal(ijl,ispin) * Clocal(isp,ic) * Clocal(isp,jc)
                    if (xdm_boolc .and. xdm_mesh == 1) then
                      Dij = 2*Dlocal(ijl,ispin)
                      grrho(:,isp,ip,ispin) = grrho(:,isp,ip,ispin) &
                                            + Dij*E1ij(:,isp) + Dij*E2ij(:,isp)
                      Ggrrho(:,:,isp,ip,ispin) = Ggrrho(:,:,isp,ip,ispin) &
                                               + 2*Dij*Kij(:,:,isp) &
                                               + Dij*F1ij(:,:,isp) &
                                               + Dij*F2ij(:,:,isp)
                      tau(:,isp,ip,ispin) = tau(:,isp,ip,ispin) &
                                          + Dij*Gij(:,isp)
                    end if
! ^LML
                 end do
              end do

           end do

!          ilc(ic) == il
           ijl = idx_ijl(il,ilc(ic))
           
!          Loop over sub-points
           do isp = 1,nsp
              if (xdm_boolc .and. xdm_mesh == 1) then
                Gij(:,isp) = Glocal(:,isp,ic) * Glocal(:,isp,ic)
                E1ij(:,isp) = Clocal(isp,ic) / dsqrt(2._dp) * Glocal(:,isp,ic)
                E2ij(:,isp) = Glocal(:,isp,ic) * Clocal(isp,ic) / dsqrt(2._dp)
                F1ij(:,:,isp) = Clocal(isp,ic) / dsqrt(2._dp) * Flocal(:,:,isp,ic)
                F2ij(:,:,isp) = Flocal(:,:,isp,ic) * Clocal(isp,ic) / dsqrt(2._dp)
                do ik = 1,3
                  Kij(ik,ik,isp) = Glocal(ik,isp,ic) * Glocal(ik,isp,ic)
                end do
              end if
! ^LML
              do ispin = 1,nspin
                 rhoscf(isp,ip,ispin) = rhoscf(isp,ip,ispin) + &
                      Dlocal(ijl,ispin) * 0.5_dp * Clocal(isp,ic) ** 2
                 if (xdm_boolc .and. xdm_mesh == 1) then
                   Dij = Dlocal(ijl,ispin)
                   grrho(:,isp,ip,ispin) = grrho(:,isp,ip,ispin) &
                                         + Dij*E1ij(:,isp) + Dij*E2ij(:,isp)
                   Ggrrho(:,:,isp,ip,ispin) = Ggrrho(:,:,isp,ip,ispin) &
                                            + 2*Dij*Kij(:,:,isp) &
                                            + Dij*F1ij(:,:,isp) &
                                            + Dij*F2ij(:,:,isp)
                   tau(:,isp,ip,ispin) = tau(:,isp,ip,ispin) &
                                       + Dij*Gij(:,isp)
                  end if
! ^LML
              end do
              
           end do
           
        end do

     else

        lasta = 0
        lastop = 0
! ^LML
!       Store values
        do ic = 1 , nc
           imp = endpht(ip-1) + ic
           i   = lstpht(imp)
           il  = ilocal(lstpht(imp))
           ia  = iaorb(lstpht(imp))
           is = isa(ia)
           iop = listp2(imp)
           ilc(ic) = il
! ^LML

!          Retrieve phi values
           Clocal(:,ic) = dsqrt(2._dp) * phi(:,imp)
           if (xdm_boolc .and. xdm_mesh == 1) then
             if (ia.ne.lasta .or. iop.ne.lastop) then
               lasta = ia
               lastop = iop
               do isp = 1,nsp
                 dxsp(:) = xdsp(:,isp) + xdop(:,iop) - dxa(:,ia)
                 r2sp = sum(dxsp**2)
                 if (r2sp.lt.r2cut(is)) then
                   call all_phi( is, +1, dxsp, nphiloc, phia(:,isp), &
                                 grphia(:,:,isp), Ggrphia(:,:,:,isp) )
                 else
                   grphia(:,:,isp) = 0.0_dp
                   Ggrphia(:,:,:,isp) = 0.0_dp
                 endif
               enddo
             endif
             iphi = iphorb(i)
             Glocal(:,:,ic) = grphia(:,iphi,:)
             Flocal(:,:,:,ic) = Ggrphia(:,:,iphi,:)
           endif
! ^LML
 
!          Loop on second orbital of mesh point
           do jc = 1, ic - 1
              ijl = idx_ijl(il,ilc(jc))

!             Loop over sub-points
              do isp = 1,nsp
                 if (xdm_boolc .and. xdm_mesh == 1) then
                   Gij(:,isp) = Glocal(:,isp,ic) * Glocal(:,isp,jc)
                   E1ij(:,isp) = Clocal(isp,ic) / dsqrt(2._dp) * Glocal(:,isp,jc)
                   E2ij(:,isp) = Glocal(:,isp,ic) * Clocal(isp,jc) / dsqrt(2._dp)
                   F1ij(:,:,isp) = Clocal(isp,ic) / dsqrt(2._dp) * Flocal(:,:,isp,jc)
                   F2ij(:,:,isp) = Flocal(:,:,isp,ic) * Clocal(isp,jc) / dsqrt(2._dp)
                   do ik = 1,3
                     Kij(ik,ik,isp) = Glocal(ik,isp,ic) * Glocal(ik,isp,jc)
                   end do
                 end if
! ^LML
                 do ispin = 1,nspin
                    rhoscf(isp,ip,ispin) = rhoscf(isp,ip,ispin) + &
                         Dlocal(ijl,ispin) * Clocal(isp,ic) * Clocal(isp,jc)
                    if (xdm_boolc .and. xdm_mesh == 1) then
                      Dij = 2*Dlocal(ijl,ispin)
                      grrho(:,isp,ip,ispin) = grrho(:,isp,ip,ispin) &
                                            + Dij*E1ij(:,isp) + Dij*E2ij(:,isp)
                      Ggrrho(:,:,isp,ip,ispin) = Ggrrho(:,:,isp,ip,ispin) &
                                               + 2*Dij*Kij(:,:,isp) &
                                               + Dij*F1ij(:,:,isp) &
                                               + Dij*F2ij(:,:,isp)
                      tau(:,isp,ip,ispin) = tau(:,isp,ip,ispin) &
                                          + Dij*Gij(:,isp)
                    end if
! ^LML
                 end do
              end do

           end do

!          ilc(ic) == il
           ijl = idx_ijl(il,ilc(ic))
           
!          Loop over sub-points
           do isp = 1,nsp
              if (xdm_boolc .and. xdm_mesh == 1) then
                  Gij(:,isp) = Glocal(:,isp,ic) * Glocal(:,isp,ic)
                  E1ij(:,isp) = Clocal(isp,ic) / dsqrt(2._dp) * Glocal(:,isp,ic)
                  E2ij(:,isp) = Glocal(:,isp,ic) * Clocal(isp,ic) / dsqrt(2._dp)
                  F1ij(:,:,isp) = Clocal(isp,ic) / dsqrt(2._dp) * Flocal(:,:,isp,ic)
                  F2ij(:,:,isp) = Flocal(:,:,isp,ic) * Clocal(isp,ic) / dsqrt(2._dp)
                  do ik = 1,3
                    Kij(ik,ik,isp) = Glocal(ik,isp,ic) * Glocal(ik,isp,ic)
                  end do
              end if
! ^LML
              do ispin = 1,nspin
                 rhoscf(isp,ip,ispin) = rhoscf(isp,ip,ispin) + &
                       Dlocal(ijl,ispin) * 0.5_dp * Clocal(isp,ic) ** 2
                 if (xdm_boolc .and. xdm_mesh == 1) then
                   Dij = Dlocal(ijl,ispin)
                   grrho(:,isp,ip,ispin) = grrho(:,isp,ip,ispin) &
                                         + Dij*E1ij(:,isp) + Dij*E2ij(:,isp)
                   Ggrrho(:,:,isp,ip,ispin) = Ggrrho(:,:,isp,ip,ispin) &
                                            + 2*Dij*Kij(:,:,isp) &
                                            + Dij*F1ij(:,:,isp) &
                                            + Dij*F2ij(:,:,isp)
                   tau(:,isp,ip,ispin) = tau(:,isp,ip,ispin) &
                                       + Dij*Gij(:,isp)
                 end if
! ^LML
              end do
              
           end do
           
        end do

     end if
     
!    Restore iorb for next point
     do imp = 1+endpht(ip-1), endpht(ip)
        i  = lstpht(imp)
        il = ilocal(i)
        iorb(il) = -i
     end do

  end do
!$OMP end do

  if (xdm_boolc .and. xdm_mesh == 1) then
    do ip = 1,np ! Loop over mesh points
      do isp = 1,nsp ! Loop over sub-points
        do ispin = 1,nspin ! Loop over spins
          rho_r(isp,ip,ispin) = rhoscf(isp,ip,ispin)
          grrho_r(isp,ip,ispin) = sqrt(grrho(1,isp,ip,ispin)**2 &
                                     + grrho(2,isp,ip,ispin)**2 &
                                     + grrho(3,isp,ip,ispin)**2)
          laprho_r(isp,ip,ispin) = Ggrrho(1,1,isp,ip,ispin) &
                                 + Ggrrho(2,2,isp,ip,ispin) &
                                 + Ggrrho(3,3,isp,ip,ispin)
          tau_r(isp,ip,ispin) = sqrt(tau(1,isp,ip,ispin)**2 &
                                   + tau(2,isp,ip,ispin)**2 &
                                   + tau(3,isp,ip,ispin)**2)
          end do ! Loop over spins
        end do ! Loop over sub-points
      end do ! Loop over mesh points

    end if
! ^LML

! Free local memory
!$OMP critical
  deallocate( ilocal, ilc, iorb, Dlocal, Clocal )
  if (xdm_boolc .and. xdm_mesh == 1) then
    deallocate( Glocal, Flocal )
    deallocate( phia, grphia, Ggrphia )
  else
    if ( DirectPhi ) deallocate( phia )
  end if
! ^LML
!$OMP end critical

!$OMP end parallel

  call de_alloc( r2cut, 'r2cut', 'rhoofd' )
  if (ParallelLocal) then
     call de_alloc( DscfL, 'DscfL', 'rhoofd' )
  end if

#ifdef _TRACE_
  call MPI_Barrier( MPI_Comm_World, MPIerror )
  call MPItrace_event( 1000, 0 )
#endif
  call timer('rhoofd',2)

#ifdef DEBUG
  call write_debug( '    POS rhoofd' )
#endif

contains

  ! In any case will the compiler most likely inline this
  ! small routine. So it should not pose any problem.
  pure function idx_ijl(i,j) result(ij)
    integer, intent(in) :: i,j
    integer :: ij
    if ( i > j ) then
       ij = i * (i + 1)/2 + j + 1
    else
       ij = j * (j + 1)/2 + i + 1
    end if
  end function idx_ijl
  
end subroutine rhoofd

end module m_rhoofd
