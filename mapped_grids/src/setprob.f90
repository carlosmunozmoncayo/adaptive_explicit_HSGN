subroutine setprob
  use cparam_mod       
  implicit none
  integer :: iunit
  character(len=25) :: fname

  fname = 'setprob.data'
  iunit = 7

  print *, 'In setprob:'
  call opendatafile(iunit, fname)

  ! Read and assign common variables (now handled with a module)
  read(iunit,*) grav
  print *, ' grav = ', grav
  read(iunit,*) sea_level
  print *, ' sea_level = ', sea_level
  read(iunit,*) dry_tolerance
  print *, ' dry_tolerance = ', dry_tolerance
  read(iunit,*) boussEquations
  print *, ' boussEquations = ', boussEquations
  read(iunit,*) boussEDCcsq
  print *, ' boussEDCcsq = ', boussEDCcsq
  read(iunit,*) boussEDCgamma
  print *, ' boussEDCgamma = ', boussEDCgamma
  read(iunit,*) inclusion_idx
  print *, ' inclusion_idx = ', inclusion_idx
  read(iunit,*) inclusion_type
  print *, ' inclusion_type = ', inclusion_type
  read(iunit,*) srcintegrator
  print *, ' srcintegrator = ', trim(srcintegrator)
  

  close(iunit)
  

end subroutine setprob