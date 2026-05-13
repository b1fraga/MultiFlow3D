!##########################################################################
module vars
  use, intrinsic :: iso_fortran_env, only: dp => real64
!##########################################################################
          SAVE
          real(dp) :: g_dx,g_dy,g_dz,dens,Re,eps,fac,Pr,beta,Th,Tc
          real(dp) :: qzero,qstpn,forcn,gx,gy,gz,Sc_t,rrey,Tinit
          real(dp) :: flomas,rmax,alfapr,resor,fric,TI_SEM
          real(dp) :: ctime,dt,dtavg,dtsum,safety_factor,noise,Mdef
          real(dp) :: t_start_averaging1,t_start_averaging2,ubulk
          real(dp) :: xst,xen,yst,yen,zst,zen
          integer :: alfabc,niter,nswp(4),nsweep,iter,nsweeps,ntime
          integer :: ngg,nge,ngc,n_unstpt,OMP_threads,iaddinlet,ireadinlet
          integer :: ngrid_input,ngrd_gl,maxcy,iproln,irestr
          integer :: itmax,sweeps,numfile,numfile1,numfile2
          integer :: conv_sch,diff_sch,sgs_model,solver,mg_itrsch
          integer :: bc_w,bc_e,bc_s,bc_n,bc_b,bc_t,UPROF_SEM
          integer :: Tbc_w,Tbc_e,Tbc_s,Tbc_n,Tbc_b,Tbc_t
          integer :: itime,itime_start,itime_end,n_out,ITMAX_PI
          integer :: ipref,jpref,kpref,prefdom,ITMAX_SEM,NE_SEM
          integer :: pl,pl_ex,differencing,LMR,normal_inter,order
          integer :: tsteps_pt,count,jtime  !Aleks 02/2023
          character(len=80) :: keyword,L_n
          logical :: LRESTART,LIMB,SGS,PERIODIC,LENERGY,LROUGH
          logical :: pressureforce,time_averaging,reinitmean
          logical :: LPT,save_inflow,read_inflow,L_dt,LSCALAR,LAS,L_LSM
          logical :: LTRANSIENT,LNonNewt,L_LSMbase
          integer :: np
          real(dp) :: ntav1_count,ntav2_count,ntav_restart  !Aleks 04/24
!=========================================================================
      end module vars
!##########################################################################
