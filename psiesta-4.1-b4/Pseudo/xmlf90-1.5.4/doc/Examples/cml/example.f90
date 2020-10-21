 program example

   use xmlf90_wxml
   use xmlf90_cml

  integer,  parameter ::  sp = selected_real_kind(6,30)
  integer,  parameter ::  dp = selected_real_kind(14,100)
!
! NB normally you will be writting to the xml file
! from mulitple fortran files/subroutines, therefore
! type(xmlf_t) :: myfile   (below)
! would normally need to be treated as a global
! variable, either in a module or a common block.
! 
  type(xmlf_t) :: myfile
!

  integer           :: num, na
  character(len=10) :: jon
  character(len=2)  :: elements(3)
  real(kind=dp)      :: coords(3,3)
  real(kind=dp)      :: adp

  data coords(1:3,1)/0.0, 0.0, 0.0/
  data coords(1:3,2)/0.5, 0.5, 0.5/
  data coords(1:3,3)/0.4, 0.4, 0.4/

  adp=1.234567890
  na=3
  elements(1) = 'Ca'
  elements(2) = 'Si'
  elements(3) = 'O'
  num = 20
  jon = '   jon'
  
  call xml_OpenFile('output.xml', myfile, indent=.true.)

  ! Start element
  call xml_NewElement(myfile, 'foo')

  ! Add molecule
  call cmlAddMolecule(xf=myfile, natoms=na,elements=elements,coords=coords)                               

  ! Add molecule output in short style
  call cmlAddMolecule(xf=myfile, natoms=na,elements=elements,coords=coords, style='xyz3')

  ! Add molecule output in short style in user supplied format
  call cmlAddMolecule(xf=myfile, natoms=na,elements=elements,coords=coords, style='xyz3', fmt='(f12.6)')

  ! End and Close
  call xml_EndElement(myfile, 'foo')
  call xml_Close(myfile)
 
end program example
