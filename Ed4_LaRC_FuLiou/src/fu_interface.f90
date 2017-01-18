
program fu_interface
  use fulioumulti 
  use extras
  implicit none
  
  real(kind=8), allocatable, dimension(:) :: fuir
  real(kind=8), allocatable, dimension(:) :: fdir
  real(kind=8), allocatable, dimension(:) :: fusw
  real(kind=8), allocatable, dimension(:) :: fdsw
  integer np ,nl
 
 print*, 'begin Fu function test'
  call getatmosphere('../testatms/jmls.lay ',	  &  !!! Fill atmosphere Structure from an input file
 fi%nv,  &    ! # model layers
 fi%pp,  &    ! Pressure (hPa)	   @nv+1 levels
 fi%pt,  &    ! Temperature (K)	     "
 fi%ph,  &    ! H20 Mix Ratio (g/g)  "
 fi%po,  &    ! O3 Mix Ratio (g/g)   "
 fi%pts)      !Skin Temperature (K)
 print*, 'read jmls' 
 nl = fi%nv+1
 fi%pt(1:nl) = fi%pt(nl:1:-1)
 
 
 np = size(fi%pp)
 allocate(fuir(np), fdir(np), fusw(np), fdsw(np))
 
 call init
 print*, 'init Fu radiation model'
 
 call print_in_fu
 call rad(fi%nv+1, &
               fi%pp, fi%pt, fi%pts, &
               fi%ph, fi%po, &
               356.0, 0.0, 0.0, 0.0, &
               0.0,0.0,0.0, 0.0, &
               1.0, 0.3, 0.5, 0.5, 1361., &
               fuir, fdir, fusw, fdsw) 
call print_out_fu
print*, 'pass radiation routine'
print*, ' fu_interface.f90 normal end'

end program fu_interface

subroutine init
  use fuinput
  
  call set_default_options_fu ! Sets some of the more obsure inputs to reasonable values.
  !fi%lscm = (/.false., .false., .true., .false./)
end subroutine init

subroutine rad(nlev, &
               plev, tlev, tsfc, &
               qlev, o3lev, &
               co2ppmv, ch4vmr, n2ovmr, o2vmr, &
               cfc11vmr,cfc12vmr,cfc22vmr, ccl4vmr, &
               emis, albedo, coszen, fday, scon, &
               fuir, fdir, fusw, fdsw)
  use fulioumulti
  implicit none
  
  integer, intent(in) :: nlev
  real(kind=8), dimension(nlev), intent(in) :: plev
  real(kind=8), dimension(nlev), intent(in) :: tlev
  real(kind=8), intent(in) :: tsfc
  real(kind=8), dimension(nlev), intent(in) :: qlev
  real(kind=8), dimension(nlev), intent(in) :: o3lev
  real(kind=8), intent(in) :: co2ppmv
  real(kind=8), intent(in) :: ch4vmr
  real(kind=8), intent(in) :: n2ovmr
  real(kind=8), intent(in) :: o2vmr
  real(kind=8), intent(in) :: cfc11vmr
  real(kind=8), intent(in) :: cfc12vmr
  real(kind=8), intent(in) :: cfc22vmr
  real(kind=8), intent(in) :: ccl4vmr
  real(kind=8), intent(in) :: emis !greybody emission
  real(kind=8), intent(in) :: albedo
  real(kind=8), intent(in) :: coszen
  real(kind=8), intent(in) :: fday
  real(kind=8), intent(in) :: scon
  
  real(kind=8), dimension(nlev), intent(out) :: fuir
  real(kind=8), dimension(nlev), intent(out) :: fdir
  real(kind=8), dimension(nlev), intent(out) :: fusw
  real(kind=8), dimension(nlev), intent(out) :: fdsw
  
  
  
  fi%nv=nlev-1 !nv is the number of model layers 
  fi%pp = plev
  fi%pt = tlev
  fi%pts = tsfc
  fi%ph = qlev
  fi%po = o3lev
  fi%umco2 = co2ppmv
  fi%umch4 = ch4vmr*dble(1.0e06)
  fi%umn2o = n2ovmr*dble(1.0e06)
  fi%cfc_conc = (/ cfc11vmr, cfc12vmr, cfc22vmr /)
  fi%ee(:) = emis
  fi%sfcalb(:,:,:) = albedo
  fi%ss = scon*fday
  fi%u0 = coszen 
  fi%ur = coszen
  
  ! RADIATVE TRANSFER --------------------------------------------------
!   call print_in_fu		   ! PRINTS INPUTS  AS ASCII 
  call rad_multi_fu  ! CALL THE CODE !!!
!   call print_out_fu
  
  fuir = fo(3)%fuir
  fdir = fo(3)%fdir
  fusw = fo(3)%fus
  fdsw = fo(3)%fds
end
