  !##########################################################################
      module vars_pt
!##########################################################################
	SAVE
      integer np_loc,npg_loc
      logical DF,PSIcell,Lcol,Lcolwall
      integer,allocatable,dimension(:)::	ptsinproc
      integer,allocatable,dimension(:)::  ptsinproc_g                   ! ghost particle 11/2020 boyang 
	integer,allocatable,dimension(:)::  id
      double precision, allocatable, dimension(:):: xp_pt,yp_pt,zp_pt
      double precision, allocatable, dimension(:):: up_pt,vp_pt,wp_pt
      double precision, allocatable, dimension(:):: xp_loc,yp_loc,zp_loc
      double precision, allocatable, dimension(:):: uop_pt,vop_pt,wop_pt
      double precision, allocatable, dimension(:):: uop_loc,vop_loc
      double precision, allocatable, dimension(:):: wop_loc,dp_loc
      double precision, allocatable, dimension(:):: rhop_loc
      double precision, allocatable, dimension(:):: Fu,Fv,Fw
      double precision, allocatable, dimension(:):: Fpu,Fpv,Fpw
      double precision, allocatable, dimension(:):: xpold,ypold,zpold
      double precision, allocatable, dimension(:):: uopold,vopold,wopold
      double precision, allocatable, dimension(:):: dp_pt,dp_old
      double precision, allocatable, dimension(:):: rho_pt,rhop_old
!     --------------------------------
!     ghost particles
!     -------------------------------- 
      integer,allocatable,dimension(:)::  idg
      double precision, allocatable, dimension(:):: xpg_loc,ypg_loc
      double precision, allocatable, dimension(:):: zpg_loc
      double precision, allocatable, dimension(:):: uopg_loc,vopg_loc
      double precision, allocatable, dimension(:):: wopg_loc,dpg_loc
      double precision, allocatable, dimension(:):: rhopg_loc

      end module

