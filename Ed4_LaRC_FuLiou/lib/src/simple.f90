USE FULIOUMULTI
USE GENERATE_FULIOU_LEVELS ,only : gflq, generate_level_scheme
USE EXTRAS       ,only : getatmosphere, aer_scale_hgt
USE CALIPSO_OUTPUT, only : pack_sky,print_pack_sky,skyp,SKYP_TYPE

USE ICEDIRSFC,only: tau_uc !! Debug Diagnostic

implicit none 
TYPE (SKYP_TYPE) ut,tu
integer k,icase,irp,itxt,iasp,ib,ie
real psfc
integer kk

do kk = 1,100
 call set_default_options_fu ! Sets some of the more obsure inputs to reasonable values.
 fi%txt=1
 
do itxt = 0,0 !0,2 
 fi%txt=itxt
if (itxt .eq. 0 )  then
 ib=0
 ie=0
else
 ib=-20
 ie=20
endif

do iasp =  ib,ie,2

do icase = 2,2,2
fi%curvedearth= .false.
if ( icase == 2) fi%curvedearth= .true.

!

do irp=60,60
!InPut Profile Assignment
 call getatmosphere('../../testatms/jmls.lay ', &
! call getatmosphere('./testatms/jmls.lay ',&
 FI%VI%nlev,&
 FI%VI%pp,&
 FI%VI%pt,&
 FI%VI%ph,&
 FI%VI%po,&
 FI%pts)
 FI%VI%nlev = FI%VI%nlev+1  ! LAYER(getatm) to LEVEL

 FI%VI%hsfc = 0.00 !! SURFACE GEOPOTENTIAL OF FI%VI profile

 gflq%hsfc = 1563. !Meters Surface elev. of ACTUAL FOV... to nearest 120m Multiple
 gflq%hsfc = -4.0! 
 gflq%hsfc = 1500.0! 



 gflq%mode = 'CALIP'
 gflq%mode = 'CERES'

fi%isksolve= 1  ! Solver Method (0=fu 1=gwtsa) 
fi%ss	   = 1365 ! Solar Constant wm-2
fi%u0      =  1.0 ! Cosine Solar Zenith Angle
fi%ur      =  0.8 ! Cosine View Zenith Angle (for IR Radiance)

!-------Cnd 2
fi%fc(1)%dpi%ldpi = .false.
fi%fc(1)%cldfrac   = 1.00000    ! Cloud Fraction (0-1) 
fi%fc(1)%novl      =   1 

FI%VD%cldpres(1:2, 1,1) = (/704,925/)
FI%VD%cldpres(1:2, 1,1) = (/204,215/)
!FI%VD%cldpres(1:2, 1,1) = (/400,800/)

fi%fc(1)%rphase(1)    =  2.0    ! Cloud Phase 1=Water 2=Ice
fi%fc(1)%de(1) = irp
fi%fc(1)%re(1) = 15.

fi%fc(1)%asp(1) = exp(iasp*0.1) !! Fu 20006 Ice AspectRatio !!!!! NEW FOR 20010130

fi%fc(1)%tau_vis(1)       = 4	    ! Cloud Visible Optical Depth ( Minnis)
fi%fc(1)%sc(1)%mn_lin_tau =  fi%fc(1)%tau_vis(1) 



!Surface Properties --------------------------------------------------

!Allow different albedos for Aerosol Vs. NO Aerosol cases , And for each Clear/Cloud Conditions
fi%sfcalb(1:18,1,0)  = 0.0 ! Clear sky -Spectral Surface Albedo SW
fi%sfcalb(1:18,2,0)  = 0.0 ! Pristine sky -Spectral Surface Albedo SW
fi%sfcalb(1:18,1,1:)  = 0.0  ! CLOUDY w/AOT  sky -Spectral Surface Albedo SW
fi%sfcalb(1:18,2,1:)  = 0.0  ! CLOUDY w/o AOT sky -Spectral Surface Albedo SW

fi%ee(1:12)  = 0.99 ! Spectral Surface Emissivity LW

!Aerosols ------------------------------------------------------------
fi%nac		     = 1	   ! 2 aerosol types 
fi%itps(1)	     = 2	   ! Continental see types (1-18)
!fi%itps(2)	     = 11	   ! Soot	  see types (1-18)

fi%n_atau	      = 1	   ! 1 wavelength input for aerosols
fi%a_wli(1)	      = 0.63	   ! AOT wavelength(microns) of a_taus
fi%a_taus(1,1)	      =  0.0001	   ! AOT for constituent 1
!fi%a_taus(1,2)	      =  0.0001	   ! AOT for constituent 2

!----------------------------------------------------------------------

 call generate_level_scheme !! Define model Fixed layer structure pre-cloud by fixed DZ intervals...
!call print_vla_in 
 call prepare_model_profile_fu !! CALL After all FI%VD and FI%VI structures are defined.
 call vla_interface_fu     ! uses FI%VO !! Assign Model ATM Profile and CLD Levels
! call print_vla_out

!Aerosol Profile (after fi%pp is created )-----------------------------

 call aer_scale_hgt(fi%nv,fi%pp,3.0,fi%aprofs(1:fi%nv,1) )
 call aer_scale_hgt(fi%nv,fi%pp,3.0,fi%aprofs(1:fi%nv,2) )
! RADIATVE TRANSFER --------------------------------------------------

if (kk == 100) call print_in_fu		   ! PRINTS INPUTS  AS ASCII 
!fi%swonlycomp=.true.
 call rad_multi_fu  ! CALL THE CODE !!!

if (kk == 100 ) call print_out_fu		   ! PRINTS OUTPUTS  AS ASCII


!stop ' normal end'

! call print_out_hr	!PRINTS ASCII HEATING RATE PROFILES

!--------------------------------------------------------------------
!  call pack_sky
!  call print_pack_sky
!  ut= skyp
!  write(11) ut
if ( fi%ierror /= 0 ) print*,' fi%ierror=',fi%ierror
!print'(f5.1,3f10.2,f10.3,3f10.2,10f8.2)',&
!fi%fc(1)%de(1),&
!ftoa(2)%olr,&
!ftoa(2)%swdn,&
!ftoa(2)%swup, &
!ftoa(2)%swup/ftoa(2)%swdn,&
!fsfc(2)%swdn,&
!fsfc(2)%swdir,&
!fsfc(2)%swdif,&
!tau_uc(:,1)

print'(I4,f8.2,f8.0,f8.2,2f8.2,3f8.3)',fi%txt, &
fi%fc(1)%tau_vis(1),&
fi%fc(1)%de(1),&
fi%fc(1)%asp(1), &
ftoa(2)%swdn,&
ftoa(2)%swup, &
ftoa(2)%swalb, &
fos(2)%swalbtoa(10),&
fos(2)%swalbtoa(18)

enddo
enddo
enddo
enddo


enddo !kk
end
