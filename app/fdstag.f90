program fdstag
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use multidata, only: id_unst, dom, i_unst, dom_id, dom_indid, imbinblk, j_unst, dom_ad, &
                       k_unst, nbp, nbpmax, num_domains
  use multiflow3d_mpi, only: ierr, mpi_comm_world, myrank, &
                             init_parallelisation, end_parallelisation
#if USE_JSON == 1
  use vars, only: bc_b, bc_t, bc_e, bc_s, bc_n, bc_w, beta, conv_sch, &
                  dens, diff_sch, differencing, dt, eps, fric, g_dx, g_dy, &
                  g_dz, gx, gy, gz, iproln, irestr, itime_end, itmax_pi, &
                  itmax_sem, keyword, l_dt, l_LSM, l_LSMbase, las, lenergy, &
                  limb, lmr, lnonnewt, lpt, lrough, lscalar, lrestart, &
                  ltransient, maxcy, mg_itrsch, n_out, niter, n_unstpt, &
                  ngrid_input, noise, np, numfile, pl_ex, pr, pressureforce, &
                  rrey, safety_factor, save_inflow, sgs, sgs_model, solver, &
                  sweeps, t_start_averaging1, t_start_averaging2, tbc_b, &
                  tbc_e, tbc_n, tbc_s, tbc_t, tc, th, ti_sem, time_averaging, &
                  tinit, tsteps_pt, ubulk, uprof_sem, l_n, re, read_inflow, &
                  reinitmean, sc_t, tbc_w, nswp
  use io, only: read_mdmap, read_control_file
  use multidata, only: dom, id_unst, i_unst, j_unst, k_unst, nbp
#else
  use vars, only: l_lsm, limb, lpt, lrestart, lrough, noise, &
       numfile, solver, tc, th, ti_sem, time_averaging, tinit, &
       tsteps_pt, ubulk, uprof_sem, nswp, niter, n_out, sweeps, tbc_w, &
       eps,  t_start_averaging1, t_start_averaging2, tbc_b, tbc_e, tbc_n, &
       tbc_s, tbc_t, reinitmean, rrey, safety_factor, save_inflow, sc_t, sgs, &
       sgs_model, solver, n_unstpt, ngrid_input, noise, np, pl_ex, pr, &
       pressureforce, re, read_inflow, lmr, lnonnewt, lpt, lrestart, lrough, &
       lscalar, ltransient, maxcy, mg_itrsch, itmax_pi, itmax_sem, keyword, &
       l_dt, l_lsm, l_lsmbase, l_n, las, lenergy, limb, differencing, dt, &
       fric, g_dx, g_dy, g_dz, gx, gy, gz, iproln, irestr, itime_end, &
       numfile, bc_b, bc_e, bc_n, bc_s, bc_t, bc_w, beta, conv_sch, dens, diff_sch
  use io, only: read_mdmap, read_control
#endif
  implicit none
#if USE_JSON == 1
  integer :: ib
#endif
  call init_parallelisation

  call read_mdmap(numfile,dom_id, dom_indid,dom_ad,imbinblk,dom,nbp,num_domains)
#if USE_JSON == 1
  call read_control_file("inputs/control.cin",Keyword,L_n,&
                         g_dx,g_dy,g_dz,ubulk,rrey,Pr,Sc_t,beta,&
                         gx,gy,gz,dens,conv_sch,diff_sch,differencing,&
                         solver,ngrid_input,mg_itrsch,&
                         maxcy,irestr,iproln,&
                         dt,safety_factor,eps,fric,&
                         L_dt,lrestart,reinitmean,LTRANSIENT,&
                         sweeps,itime_end,n_out,tsteps_pt,niter,&
                         nswp(1),nswp(2),nswp(3),nswp(4),bc_w,&
                         bc_e,bc_s,bc_n,&
                         bc_b,bc_t,&
                         save_inflow,time_averaging,SGS,&
                         ITMAX_PI,UPROF_SEM,ITMAX_SEM,&
                         TI_SEM,t_start_averaging1,t_start_averaging2,&
                         noise,Th,Tc,Tinit,SGS_model,LMR,pl_ex,&
                         LIMB,LENERGY,LROUGH,LPT,l_LSM,L_LSMbase,LSCALAR,LAS,LNonNewt,&
                         Tbc_w,Tbc_e,Tbc_s,Tbc_n,&
                         Tbc_b,Tbc_t,n_unstpt,&
                         id_unst,i_unst,j_unst,k_unst)
  Re = 1.0_dp/rrey
  if (bc_w==5) pressureforce=.true.
  if (.not.LPT) np=0

  if (trim(L_n)=="n") fric=fric**0.33_dp


  do ib=1,nbp
    dom(ib)%bc_west=bc_w
    dom(ib)%bc_east=bc_e
    dom(ib)%bc_south=bc_s
    dom(ib)%bc_north=bc_n
    dom(ib)%bc_bottom=bc_b
    dom(ib)%bc_top=bc_t
    dom(ib)%Tbc_west=Tbc_w
    dom(ib)%Tbc_east=Tbc_e
    dom(ib)%Tbc_south=Tbc_s
    dom(ib)%Tbc_north=Tbc_n
    dom(ib)%Tbc_bottom=Tbc_b
    dom(ib)%Tbc_top=Tbc_t
    dom(ib)%ngrid=ngrid_input
    if (dom(ib)%bc_west==7) read_inflow=.true.
  end do

  if(differencing==3 .and. pl_ex/=2) then
     pl_ex=2
     print*,"error: you select WENO but do not assign", &
          "  correct number of ghost planes,now it is corrected to 2"
  end if

  if(SGS .and. sgs_model==3 .and. pl_ex/=2) then
     pl_ex=2
     print*,"error: you select 1-EQN model but do not assign", &
          "  correct number of ghost planes,now it is corrected to 2"
  end if

  if (L_LSM .and. solver==1) then
     if (myrank==0) then
        print*,"Error: SIP solver not presently compatible with LSM"
     end if
     stop
  end if

  if (L_LSM .and. differencing/=3) then
     if (myrank==0) then
        print*,"Error: WENO differencing must be used with LSM"
     end if
     stop
  end if
#else
  call read_control(eps,n_out,niter,sweeps,tbc_w,nswp,th,tc,ti_sem, &
       time_averaging,tinit,tsteps_pt,ubulk,UPROF_SEM,t_start_averaging1, &
       t_start_averaging2,Tbc_e,Tbc_s,Tbc_n,Tbc_b,Tbc_t,reinitmean,rrey, &
       safety_factor,save_inflow,Sc_t,SGS,sgs_model,solver,n_unstpt, &
       ngrid_input,noise,np,pl_ex,pr,pressureforce,re,read_inflow,LMR,LNonNewt, &
       LPT,LRESTART,LROUGH,LSCALAR,LTRANSIENT,maxcy,mg_itrsch,ITMAX_PI, &
       ITMAX_SEM,keyword,L_dt,L_LSM,L_LSMbase,L_n,LAS,LENERGY,LIMB, &
       differencing,dt,fric,g_dx,g_dy,g_dz,gx,gy,gz,iproln,irestr,itime_end, &
       bc_w,bc_e,bc_s,bc_n,bc_b,bc_t,beta,conv_sch,dens,diff_sch,id_unst, &
       i_unst,j_unst,k_unst,nbp,dom)
#endif
  call read_infodom

  call alloc_dom

  call localparameters

  call initial

  call initflowfield

  if (LROUGH)  then                             !Richard2015
     if (.not.LRESTART) then
        call init_rough
     else
        call rough_restart
     end if
  end if

  if (LIMB) call imb_initial
  if (LIMB) call PartLocMPI                         !Pablo2015

  call iniflux

  if(.not.LRESTART) then
     if (time_averaging) then
        call update_mean
        if (noise>0.0_dp) call add_noise(noise)
     end if
  end if

  if (LPT) then                     !Brunho2013
     if (myrank==0)  open(unit=202, file="particle_log")
     call init_particle
  end if

  if ((solver==2).and.(.not.L_LSM)) call coeff

  call MPI_BARRIER (MPI_COMM_WORLD,ierr)
  if(myrank==0) then
     write (numfile,*) "============START ITERATIONS========="
     write (6,*) "============START ITERATIONS========="
  end if

  call flosol

  call end_parallelisation

end program fdstag
