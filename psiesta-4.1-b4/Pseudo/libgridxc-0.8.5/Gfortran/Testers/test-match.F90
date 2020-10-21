PROGRAM test_match

  ! Program to compare the Soler-Balbas and the libxc
  ! numerical results.
  ! It uses a radial charge density, and atomxc as driver

  ! nspin and the relativistic flag are hard-wired
  ! Note the exchange differences for PBE in relativistic mode.
  ! Relativistic exchange is not implemented for PBE in libxc.
  ! In the Soler-Balbas package, it is inherited from LDA_X.
  
  ! Used module procedures
  USE gridXC, only: atomXC => gridxc_atomXC
  USE gridXC, only: setXC  => gridxc_setXC, setXC_libxc => gridxc_setXC_libxc
  USE gridXC, only: gridxc_init

  ! Used module parameters
  USE gridXC, only: dp

  implicit none

  ! Tester parameters
  integer, parameter:: irel  =  1 ! Relativistic? 0=>no, 1=>yes
  integer, parameter:: nSpin =  2 ! Number of spin components
  integer, parameter:: nr    = 101 ! Number of radial points
  integer, parameter:: n1cut =  8 ! Cutoff parameter
  integer, parameter:: n2cut =  2 ! Cutoff parameter:
                                  !    fCut(r)=(1-(r/rMax)**n1cut)**n2cut
  real(dp),parameter:: dWidth = 2._dp ! Width of density distribution, in Bohr
  real(dp),parameter:: Qtot = 10._dp  ! Integral of density distribution
  real(dp),parameter:: spinPol= 2._dp ! Integral of densUp - densDown
  real(dp),parameter:: rMax = 20._dp  ! Cutoff radius, in Bohr
  real(dp),parameter:: deltaDens = 1.e-8_dp  ! Finite diff. change
  real(dp),parameter:: densMin  = 1.e-9_dp  ! Min. density to proceed


  ! Tester variables and arrays
  integer :: iDelta, ir, irmax, ismax, iSpin, one, two
  real(dp):: avgDiffVxc, dDensdr, dens(nr,nSpin), dens0(nr,nSpin), &
             d0tot, d0(nSpin), dEdDens, dDens, diffVxc, &
             Dc, Dc0, dr, dVol, Dx, Dx0, Ec, Ec0, Ex, Ex0, &
             kf, kg, maxDiffVxc, pi, r, rMesh(nr), &
             Vxc(nr,nSpin), Vxc0(nr,nSpin),  wr

  
  ! Find radial mesh points and gaussian density
  pi = acos(-1._dp)
  d0tot = Qtot / (2*pi*dWidth**2)**1.5_dp    ! Total density at origin
  if (nSpin==1) then
    d0(1) = d0tot
  else
    one = 1   ! A silly thing to satisfy the compiler when nSpin=1
    two = 2
    d0(one) = d0tot * (Qtot + spinPol) / Qtot / 2 ! Spin up density at origin
    d0(two) = d0tot * (Qtot - spinPol) / Qtot / 2 ! Spin down density at origin
  end if
  dr = rmax / (nr-1)                      ! Interval between radial points
  do ir = 1,nr
    rMesh(ir) = dr * (ir-1)               ! Radial point values
    dens0(ir,:) = DensOfR( d0(:), rMesh(ir) )
  end do

  ! Find exchange and correlation energy and potential from radial density
  call setXC(1,['GGA'],['PBE'], [1.0_dp], [1.0_dp])
!  call setXC(1,['LDA'],['PZ'], [1.0_dp], [1.0_dp])

  call atomXC( irel, nr, nr, rMesh, nSpin, dens0, Ex0, Ec0, Dx0, Dc0, Vxc0 )
  print *, "Ex, Ec:", Ex0, Ec0
  print "(6f10.5)",  Vxc0

!  call setXC_libxc(2,[1,9])
  call setXC_libxc(2,[101,130])
  
  call atomXC( irel, nr, nr, rMesh, nSpin, dens0, Ex0, Ec0, Dx0, Dc0, Vxc0 )
  print *, "Ex, Ec:", Ex0, Ec0
  print "(6f10.5)",  Vxc0

  
  
CONTAINS

FUNCTION DensOfR( d0, r )

  ! Returns a radial density distribution

  implicit none
  real(dp),intent(in):: d0(nSpin)  ! Density at center of charge distribution
  real(dp),intent(in):: r          ! Distance to center of charge distribution
  real(dp)           :: DensOfR(nSpin)  ! Electron density

  ! Use a simple gaussian distribution
  DensOfR = d0 * exp(-r**2/2/dWidth**2)

  ! Impose a smooth radial cutoff
  DensOfR = DensOfR * ( 1 - (r/rMax)**n1cut )**n2cut

END FUNCTION DensOfR

END PROGRAM test_match

