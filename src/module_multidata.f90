!##########################################################################
      module multidata
        use, intrinsic :: iso_fortran_env, only: dp => real64
        use mpi_f08
!##########################################################################
          SAVE
          integer :: nbp,nbpmax,num_domains
          integer :: rdivmax,idom,jdom,kdom
          integer,allocatable,dimension(:) :: rdiv
          integer,allocatable,dimension(:,:) :: rdv
          integer,allocatable,dimension(:) :: dom_id,dom_ad,dom_indid
          integer,allocatable,dimension(:) :: i_unst,j_unst,k_unst
          integer,allocatable,dimension(:) :: id_unst
          integer,allocatable,dimension(:) :: imbinblk  !Pablo
          real(dp), allocatable,dimension(:,:)::xcor,ycor,zcor  !Pablo

          type multidom
              integer :: inext,iprev,jnext,jprev,knext,kprev
              integer :: corprev1,corprev2,corprev3,corprev4
              integer :: cornext1,cornext2,cornext3,cornext4
              integer :: edgprev1,edgprev2,edgprev3
              integer :: edgprev4,edgprev5,edgprev6
              integer :: edgnext1,edgnext2,edgnext3
              integer :: edgnext4,edgnext5,edgnext6
              integer :: per_in,per_ip,per_jn,per_jp,per_kn,per_kp
              integer :: ttc_i,ttc_j,ttc_k,ttc_ijk,ngrid
              integer :: isp,iep,jsp,jep,ksp,kep
              integer :: isu,ieu,jsu,jeu,ksu,keu
              integer :: isv,iev,jsv,jev,ksv,kev
              integer :: isw,iew,jsw,jew,ksw,kew
              integer :: nwork,nvars
              integer :: mximb
              type(MPI_Request) :: rq_m1, rq_p1, rq_m2, rq_p2, rq_m3, rq_p3
              type(MPI_Request) :: rq_c1m, rq_c2m, rq_c3m, rq_c4m
              type(MPI_Request) :: rq_c1p, rq_c2p, rq_c3p, rq_c4p
              type(MPI_Request) :: rq_e1m, rq_e2m, rq_e3m, rq_e4m, rq_e5m, rq_e6m
              type(MPI_Request) :: rq_e1p, rq_e2p, rq_e3p, rq_e4p, rq_e5p, rq_e6p
              real(dp)    :: xsl,ysl,zsl,xel,yel,zel,dx,dy,dz
              real(dp),pointer,dimension(:) :: tauw
              real(dp), pointer, dimension(:,:,:) :: S,So,Sm,Stm
              real(dp), pointer, dimension(:,:,:) :: sfactor
              real(dp), pointer, dimension(:) :: x,y,z,xc,yc,zc
              real(dp), pointer, dimension(:,:,:) :: u,v,w,p,pp
              real(dp), pointer, dimension(:,:,:) :: ksgs,ksgso
              real(dp), pointer, dimension(:,:,:) :: eps,epso
              real(dp), pointer, dimension(:,:,:) :: T,To,Tm,Ttm
              real(dp), pointer, dimension(:,:,:) :: su,ap,sup
              real(dp), pointer, dimension(:,:,:) :: ae,aw,as,an
              real(dp), pointer, dimension(:,:,:) :: at,ab
              real(dp), pointer, dimension(:,:,:) :: um,vm,wm
              real(dp), pointer, dimension(:,:,:) :: pm,ppm,vis
              real(dp), pointer, dimension(:,:,:) :: uum,vvm,wwm
              real(dp), pointer, dimension(:,:,:) :: uvm,uwm,vwm
              real(dp), pointer, dimension(:,:,:) :: ustar,vstar
              real(dp), pointer, dimension(:,:,:) :: wstar
              real(dp), pointer, dimension(:,:,:) :: uo,uoo,vo
              real(dp), pointer, dimension(:,:,:) :: voo,wo,woo
              real(dp), pointer, dimension(:,:,:) :: stfcinf
              real(dp), pointer, dimension(:,:,:) :: facp1,facp2
              real(dp), pointer, dimension(:,:,:) :: facm1,facm2
              real(dp), pointer, dimension (:) :: dh1,dh2,dh3
              integer, pointer, dimension (:,:,:,:) :: ndimb,imbbodynum
              integer, pointer, dimension (:) :: faz,cntp
              integer, pointer, dimension (:,:,:) :: ntav1,ntav2
              integer, dimension (26) :: tg
              integer, pointer, dimension (:,:,:) :: ibfactor
              real(dp), pointer, dimension(:)     :: cof
              real(dp), pointer, dimension(:,:)   :: tauwe,tauww
              real(dp), pointer, dimension(:,:)   :: tauws,tauwn
              real(dp), pointer, dimension(:,:)   :: tauwt,tauwb
              real(dp), pointer, dimension(:,:)   :: tauwe2,tauww2
              real(dp), pointer, dimension(:,:)   :: tauws2,tauwn2
              real(dp), pointer, dimension(:,:)   :: tauwt2,tauwb2
              real(dp), pointer, dimension(:) :: sendb_m1,sendb_p1
              real(dp), pointer, dimension(:) :: recvb_m1,recvb_p1
              real(dp), pointer, dimension(:) :: sendb_m2,sendb_p2
              real(dp), pointer, dimension(:) :: recvb_m2,recvb_p2
              real(dp), pointer, dimension(:) :: sendb_m3,sendb_p3
              real(dp), pointer, dimension(:) :: recvb_m3,recvb_p3
              real(dp), pointer, dimension(:)::sc1m,sc1p,rc1m,rc1p
              real(dp), pointer, dimension(:)::sc2m,sc2p,rc2m,rc2p
              real(dp), pointer, dimension(:)::sc3m,sc3p,rc3m,rc3p
              real(dp), pointer, dimension(:)::sc4m,sc4p,rc4m,rc4p
              real(dp), pointer, dimension(:)::se1m,se1p,re1m,re1p
              real(dp), pointer, dimension(:)::se2m,se2p,re2m,re2p
              real(dp), pointer, dimension(:)::se3m,se3p,re3m,re3p
              real(dp), pointer, dimension(:)::se4m,se4p,re4m,re4p
              real(dp), pointer, dimension(:)::se5m,se5p,re5m,re5p
              real(dp), pointer, dimension(:)::se6m,se6p,re6m,re6p
              real(dp), pointer, dimension(:,:,:) ::d1,dphi_dxplus
              real(dp), pointer, dimension(:,:,:) ::dphi_dyplus, &
        dphi_dzplus,dphi_dxminus,dphi_dyminus,dphi_dzminus

              real(dp), pointer, dimension (:,:) :: u_unst,v_unst
              real(dp), pointer, dimension (:,:) :: w_unst
              real(dp), pointer, dimension (:,:) :: um_unst,vm_unst
              real(dp), pointer, dimension (:,:) :: wm_unst
              real(dp), pointer, dimension (:,:) :: p_unst,pm_unst
              real(dp), pointer, dimension (:,:) :: ksgs_unst
              real(dp), pointer, dimension (:,:) :: eps_unst
              real(dp), pointer, dimension (:,:) :: T_unst,Tm_unst  !Aleks 04/24
!============================== LSM VARIABLES ============================
              real(dp), pointer, dimension(:,:,:) :: phi_init, &
        phi_new,phi_reinit,phi,dphi_dx,dphi_dy,dphi_dz,s_phi0,h_phi, &
        dens,mu,phim,abs_dphi_check
              real(dp), pointer, dimension(:)     :: dens_mg
              real(dp), pointer, dimension(:) :: sendb_m,sendb_p
              real(dp), pointer, dimension(:) :: recvb_m,recvb_p
              integer, pointer, dimension (:) :: ijkp_lsm
              integer :: tot
              integer :: niul,njul,nkul,nivl,njvl,nkvl,niwl,njwl,nkwl
              integer :: nipl,njpl,nkpl,nigl,njgl,nkgl
              integer :: nipl2,njpl2,nkpl2
              real(dp), pointer, dimension(:,:,:) :: resmax
              real(dp), pointer, dimension(:,:,:) :: resfact
              real(dp), pointer, dimension(:,:,:) :: resfact1
!==========================================================================
              integer :: bc_west,bc_east,bc_south,bc_north,bc_bottom,bc_top
              integer :: Tbc_west,Tbc_east,Tbc_south,Tbc_north
              integer :: Tbc_bottom,Tbc_top
              logical :: coarse_ng,fine_ng
          end type multidom

          type (multidom), pointer, dimension(:) :: dom

      end module multidata
!##########################################################################
