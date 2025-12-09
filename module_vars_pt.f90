      !##########################################################################
      module vars_pt
use, intrinsic :: iso_fortran_env, only: dp => real64
!##########################################################################
          SAVE
          integer :: np_loc,npg_loc
          logical :: DF,PSIcell,Lcol,Lcolwall
          real(dp) :: k_n
          integer,allocatable,dimension(:)::    ptsinproc
          integer,allocatable,dimension(:)::  ptsinproc_g                   ! ghost particle 11/2020 boyang
          integer,allocatable,dimension(:)::  id
          real(dp), allocatable, dimension(:):: xp_pt,yp_pt,zp_pt
          real(dp), allocatable, dimension(:):: up_pt,vp_pt,wp_pt
          real(dp), allocatable, dimension(:):: xp_loc,yp_loc,zp_loc
          real(dp), allocatable, dimension(:):: uop_pt,vop_pt,wop_pt
          real(dp), allocatable, dimension(:):: uop_loc,vop_loc
          real(dp), allocatable, dimension(:):: wop_loc,dp_loc
          real(dp), allocatable, dimension(:):: rhop_loc
          real(dp), allocatable, dimension(:):: Fu,Fv,Fw
          real(dp), allocatable, dimension(:):: Fpu,Fpv,Fpw
          real(dp), allocatable, dimension(:):: xpold,ypold,zpold
          real(dp), allocatable, dimension(:):: uopold,vopold,wopold
          real(dp), allocatable, dimension(:):: dp_pt,dp_old
          real(dp), allocatable, dimension(:):: rho_pt,rhop_old
!     --------------------------------
!     ghost particles
!     --------------------------------
          integer,allocatable,dimension(:)::  idg
          real(dp), allocatable, dimension(:):: xpg_loc,ypg_loc
          real(dp), allocatable, dimension(:):: zpg_loc
          real(dp), allocatable, dimension(:):: uopg_loc,vopg_loc
          real(dp), allocatable, dimension(:):: wopg_loc,dpg_loc
          real(dp), allocatable, dimension(:):: rhopg_loc

      end module vars_pt

