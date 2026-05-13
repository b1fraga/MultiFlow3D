!##########################################################################
      module module_LSM
use, intrinsic :: iso_fortran_env, only: dp => real64
!##########################################################################
          SAVE
!============================== LSM VARIABLES =============================
          real(dp) :: reldif_LSM,length,densl,densg,nul,nug,mul,mug
          real(dp) :: cfl_lsm
          real(dp) :: uprev
          real(dp) :: slope
          real(dp) :: densip12,densjp12,denskp12
          real(dp) :: densim12,densjm12,denskm12
          real(dp) :: muip12,mujp12,mukp12,muim12,mujm12,mukm12
          real(dp) :: Mdef_w
          integer :: ntime_reinit,accuracy
          integer :: ngrid
          integer :: numfile3
          logical :: REINIT,L_LSMinit,LENDS
          logical :: L_anim_phi,L_anim_grd

!=========================================================================
      end module module_LSM
!##########################################################################
