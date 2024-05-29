!######################################################################
      module imb
!######################################################################
	  SAVE
      double precision  :: xt(5),yt(5),xdt(5),ydt(5),xddt(5),yddt(5),yto
	double precision  :: lambda,sigma,nxl

	INTEGER :: imp_proc_master,imb_block_master,forcefilej,a,b
        INTEGER :: bodynum,maxnode,master,maxnodeIBS,mdfsteps,yangcase

	LOGICAL,allocatable,dimension (:):: rotating,LSELFST,intflow

      double precision,allocatable ::Cx(:),Cxor(:),Cy(:),Cyor(:)
      double precision,allocatable ::Cz(:),Czor(:),Czz(:),Cxx(:)
      double precision,allocatable ::l2norm(:)
	double precision,allocatable,dimension(:) :: xaero,yaero,zaero
	double precision,allocatable,dimension(:) :: iniT_selfST
	double precision,allocatable,dimension(:) :: SUMtorque_ST
	double precision,allocatable,dimension(:) :: acc_selfST,massI
	double precision,allocatable,dimension(:) :: radsin,rads,pitch
	double precision,allocatable,dimension(:) :: R,reddelta
	double precision,allocatable,dimension(:) :: FX1_MASTER,FX2_MASTER
        double precision,allocatable,dimension(:) :: FXSp_MASTER
        double precision,allocatable,dimension(:) :: FXT_MASTER
	double precision,allocatable,dimension(:) :: FX3_MASTER,FX3_loc
        double precision,allocatable,dimension(:) :: nodex_loc,nodey_loc
        double precision,allocatable,dimension(:) :: nodez_loc
        double precision,allocatable,dimension(:) :: FXSp_loc,FXT_loc
        double precision,allocatable,dimension(:) :: FX1_loc,FX2_loc
        double precision,allocatable,dimension(:) :: R0_loc,alpha0_loc
        double precision,allocatable,dimension(:) :: U_Beta1_loc
        double precision,allocatable,dimension(:) :: Sp_Beta_loc
        double precision,allocatable,dimension(:) :: T_Beta_loc
        double precision,allocatable,dimension(:) :: U_Beta2_loc,zend
        double precision,allocatable,dimension(:) :: U_Beta3_loc,zini
        double precision,allocatable,dimension(:,:) :: dh1_loc,dh2_loc
        double precision,allocatable,dimension(:,:) :: dh3_loc,delvol
        double precision,allocatable,dimension(:,:) :: dh4_loc!,dh5_loc
        double precision,allocatable,dimension(:,:) :: nodexlocal
        double precision,allocatable,dimension(:,:) :: nodeylocal
        double precision,allocatable,dimension(:,:) :: nodezlocal
        double precision,allocatable,dimension(:,:) :: FX1NF,FX2NF,FX3NF
        double precision,allocatable,dimension(:,:) :: FXSpNF,FXTNF
        double precision,allocatable,dimension(:,:) :: nodex,nodey,nodez
        double precision,allocatable,dimension(:,:) :: torque,U_Beta3
        double precision,allocatable,dimension(:,:) :: U_Beta1,U_Beta2
        double precision,allocatable,dimension(:,:) :: FX1,FX2,FX3
        double precision,allocatable,dimension(:,:) :: FXSp,FXT
        double precision,allocatable,dimension(:,:) :: alpha0,R0
        double precision,allocatable,dimension(:) :: dxm,dym,dzm
	INTEGER,allocatable,dimension(:,:) :: I_nr_V,J_nr_V,K_nr_V
	INTEGER,allocatable,dimension(:,:) :: I_nr_U,J_nr_U,K_nr_U
	INTEGER,allocatable,dimension(:,:) :: I_nr_W,J_nr_W,K_nr_W
        INTEGER,allocatable,dimension(:,:) :: I_nr_Sp,J_nr_Sp,K_nr_Sp
        INTEGER,allocatable,dimension(:,:) :: I_nr_T,J_nr_T,K_nr_T
	INTEGER,allocatable,dimension(:) :: imbnumber,imb_shape,IBMnum
	INTEGER,allocatable,dimension(:) :: kmaxU,kmaxV,kmaxW,nodes
        INTEGER,allocatable,dimension(:) :: kmaxSp,kmaxT
	INTEGER,allocatable,dimension(:) :: imb_proc,imb_block,turax
	INTEGER,allocatable,dimension(:) :: lag_bod_loc,cmax,linfin
        INTEGER,allocatable,dimension(:) :: domtemp,imb_block_loc,axis
        INTEGER,allocatable,dimension(:) :: imbinblock_loc,rott_loc
        double precision,allocatable,dimension(:) :: rdiv_imb
        integer,allocatable,dimension(:) :: IBip,IBjp,IBkp !Aleks 04/23 - arrays for storing ijk values of IB points
        CHARACTER*32, allocatable, dimension (:) :: filepoints

      end module imb