!##########################################################################
      module module_LSM
!##########################################################################
          SAVE
!============================== LSM VARIABLES =============================
          double precision :: reldif_LSM,length,densl,densg,nul,nug,mul,mug
          double precision :: cfl_lsm
          double precision :: uprev
          double precision :: slope
          double precision :: densip12,densjp12,denskp12
          double precision :: densim12,densjm12,denskm12
          double precision :: muip12,mujp12,mukp12,muim12,mujm12,mukm12
          double precision :: Mdef_w
          integer :: ntime_reinit,accuracy
          integer :: ngrid
          integer :: numfile3
          logical :: REINIT,L_LSMinit,LENDS
          logical :: L_anim_phi,L_anim_grd

!=========================================================================
      end module
!##########################################################################
