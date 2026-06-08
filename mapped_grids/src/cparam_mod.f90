module cparam_mod
    !Module to hold common parameters that will be written in setprob
    !The idea is to drop common blocks
    implicit none
    double precision :: grav, sea_level, dry_tolerance
    integer :: boussEquations
    double precision :: boussEDCcsq
    double precision :: boussEDCgamma

    integer :: inclusion_idx
    character(len=25) :: srcintegrator, inclusion_type

end module cparam_mod