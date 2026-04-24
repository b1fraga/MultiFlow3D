!##########################################################################
      subroutine log_law(bound,ib)
!           Bruño Fraga Bugallo
!           Cardiff 2016
!           Schumman type boundary conditions
!##########################################################################
          use vars
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          integer :: i,j,k,ib,bound,icont,tkmax
          real(dp) :: delta,n_x,n_y,n_z,vnor,vtan
          real(dp) :: uc,vc,wc,small
          real(dp) :: aaa,bbb,const1,const2,const3,const4
          real(dp) :: ustar,yplus,ustarold,conv,Ecte,kappa

          small = 1.e-30_dp
          kappa = 0.41_dp
          Ecte  = 9.0_dp


          if (LAS.or.L_LSM) then                                            !variable density
          do i=1,dom(ib)%ttc_i
          do j=1,dom(ib)%ttc_j
          do k=1,dom(ib)%ttc_k
                      rrey=dom(ib)%mu(i,j,k)/dom(ib)%dens(i,j,k)
                  end do
                  end do
                  end do
          end if

          SELECT CASE (bound)

            CASE (1)
              delta=dom(ib)%dx
              n_x=1.0_dp
              n_y=0.0_dp
              n_z=0.0_dp
              i=dom(ib)%isp
              do k=dom(ib)%ksp-1,dom(ib)%kep+1
                  do j=dom(ib)%jsp-1,dom(ib)%jep+1
                      uc=0.5_dp*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                      vc=0.5_dp*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                      wc=0.5_dp*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                      vnor=n_x*uc+n_y*vc+n_z*wc
                      vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))
                      if (dom(ib)%bc_west==62) then               !rough wall
                      ustar = kappa*vtan / log(30.d0*delta/fric)
                      yplus = delta*ustar/rrey
                      else if (dom(ib)%bc_west==61) then               !smooth wall
                      conv  = 1.0_dp
                      ustar = 1.0_dp
                      icont = 0
                      tkmax = 10
                      do while ((conv>(1.d-3)).and. &
                (icont<tkmax))
                          icont=icont+1
                          ustarold = ustar
                          yplus = max(1.0000001_dp,delta*ustar/rrey)
                          ustar = vtan*kappa/log(Ecte*yplus)
                          conv  = abs((ustar-ustarold)/ustar)
                      end do
                      if (icont==tkmax) then
                      print*,"USTAR DOESNT CONVERGE"
                      print*,"Tang. vel. = ",vtan
                      stop
                      end if
                      end if
                      if (yplus<11.067d0) then
                      dom(ib)%tauww(j,k)=rrey * vtan / delta
                      else
                      dom(ib)%tauww(j,k)=ustar**2.0_dp             !tau_1
                      dom(ib)%tauww2(j,k)=dom(ib)%tauww(j,k)/vtan   !needs to be multiplied by a velocity component to provide tau_1j
                      end if
                  end do
              end do

            CASE (2)
              delta=dom(ib)%dx
              n_x=-1.0_dp
              n_y=0.0_dp
              n_z=0.0_dp
              i=dom(ib)%iep
              do k=dom(ib)%ksp-1,dom(ib)%kep+1
                  do j=dom(ib)%jsp-1,dom(ib)%jep+1
                      uc=0.5_dp*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                      vc=0.5_dp*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                      wc=0.5_dp*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                      vnor=n_x*uc+n_y*vc+n_z*wc
                      vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))
                      if (dom(ib)%bc_east==62) then               !rough wall
                      ustar = kappa*vtan / log(30.d0*delta/fric)
                      yplus = delta*ustar/rrey
                      else if (dom(ib)%bc_east==61) then               !smooth wall
                      conv  = 1.0_dp
                      ustar = 1.0_dp
                      icont = 0
                      tkmax = 10
                      do while ((conv>(1.d-3)).and. &
                (icont<tkmax))
                          icont=icont+1
                          ustarold = ustar
                          yplus = max(1.0000001_dp,delta*ustar/rrey)
                          ustar = vtan*kappa/log(Ecte*yplus)
                          conv  = abs((ustar-ustarold)/ustar)
                      end do
                      if (icont==tkmax) then
                      print*,"USTAR DOESNT CONVERGE"
                      print*,"Tang. vel. = ",vtan
                      stop
                      end if
                      end if
                      if (yplus<11.067d0) then
                      dom(ib)%tauwe(j,k)=rrey * vtan / delta
                      else
                      dom(ib)%tauwe(j,k)=ustar**2.0_dp             !tau_1
                      dom(ib)%tauwe2(j,k)=dom(ib)%tauwe(j,k)/vtan       !units m/s
                      end if
                  end do
              end do

            CASE(3)
              delta=dom(ib)%dy
              n_x=0.0_dp
              n_y=1.0_dp
              n_z=0.0_dp
              j=dom(ib)%jsp

              do k=dom(ib)%ksp-1,dom(ib)%kep+1
                  do i=dom(ib)%isp-1,dom(ib)%iep+1
                      uc=0.5_dp*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                      vc=0.5_dp*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                      wc=0.5_dp*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                      vnor=n_x*uc+n_y*vc+n_z*wc
                      vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))
                      if (dom(ib)%bc_south==62) then              !rough wall
                      ustar = kappa*vtan / log(30.d0*delta/fric)
                      yplus = delta*ustar/rrey
                      else if (dom(ib)%bc_south==61) then              !smooth wall
                      conv  = 1.0_dp
                      ustar = 1.0_dp
                      icont = 0
                      tkmax = 10
                      do while ((conv>(1.d-3)).and. &
                (icont<tkmax))
                          icont=icont+1
                          ustarold = ustar
                          yplus = max(1.0000001_dp,delta*ustar/rrey)
                          ustar = vtan*kappa/log(Ecte*yplus)
                          conv  = abs((ustar-ustarold)/ustar)
                      end do
                      if (icont==tkmax) then
                      print*,"USTAR DOESNT CONVERGE"
                      print*,"Tang. vel. = ",vtan
                      stop
                      end if
                      end if
                      if (yplus<11.067d0) then
                      dom(ib)%tauws(i,k)=rrey * vtan / delta
                      else
                      dom(ib)%tauws(i,k)=ustar**2.0_dp             !tau_1
                      dom(ib)%tauws2(i,k)=dom(ib)%tauws(i,k)/vtan
                      end if
                  end do
              end do


            CASE (4)
              delta=dom(ib)%dy
              n_x=0.0_dp
              n_y=-1.0_dp
              n_z=0.0_dp
              j=dom(ib)%jep
              do k=dom(ib)%ksp-1,dom(ib)%kep+1
                  do i=dom(ib)%isp-1,dom(ib)%iep+1
                      uc=0.5_dp*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                      vc=0.5_dp*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                      wc=0.5_dp*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                      vnor=n_x*uc+n_y*vc+n_z*wc
                      vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))
                      if (dom(ib)%bc_north==62) then              !rough wall
                      ustar = kappa*vtan / log(30.d0*delta/fric)
                      yplus = delta*ustar/rrey
                      else if (dom(ib)%bc_north==61) then              !smooth wall
                      conv  = 1.0_dp
                      ustar = 1.0_dp
                      icont = 0
                      tkmax = 10
                      do while ((conv>(1.d-3)).and. &
                (icont<tkmax))
                          icont=icont+1
                          ustarold = ustar
                          yplus = max(1.0000001_dp,delta*ustar/rrey)
                          ustar = vtan*kappa/log(Ecte*yplus)
                          conv  = abs((ustar-ustarold)/ustar)
                      end do
                      if (icont==tkmax) then
                      print*,"USTAR DOESNT CONVERGE"
                      print*,"Tang. vel. = ",vtan
                      stop
                      end if
                      end if
                      if (yplus<11.067d0) then
                      dom(ib)%tauwn(i,k)=rrey * vtan / delta
                      else
                      dom(ib)%tauwn(i,k)=ustar**2.0_dp             !tau_1
                      dom(ib)%tauwn2(i,k)=dom(ib)%tauwn(i,k)/vtan
                      end if
                  end do
              end do

            CASE (5)

              delta=dom(ib)%dz
              n_x=0.0_dp
              n_y=0.0_dp
              n_z=1.0_dp
              k=dom(ib)%ksp

              do j=dom(ib)%jsp-1,dom(ib)%jep+1
                  do i=dom(ib)%isp-1,dom(ib)%iep+1
                      uc=0.5_dp*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                      vc=0.5_dp*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                      wc=0.5_dp*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                      vnor=n_x*uc+n_y*vc+n_z*wc
                      vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))
                      if (dom(ib)%bc_bottom==62) then             !rough wall
                      ustar = kappa*vtan / log(30.d0*delta/fric)
                      yplus = delta*ustar/rrey
                      else if (dom(ib)%bc_bottom==61) then             !smooth wall
                      conv  = 1.0_dp
                      ustar = 1.0_dp
                      icont = 0
                      tkmax = 10

                      do while ((conv>(1.d-3)).and. &
                (icont<tkmax))
                          icont=icont+1
                          ustarold = ustar
                          yplus = max(1.0000001_dp,delta*ustar/rrey)
                          ustar = vtan*kappa/log(Ecte*yplus)
                          conv  = abs((ustar-ustarold)/ustar)
                      end do
!               if (dom_id(ib).eq.0.and.j.eq.20)
!     &         write(6,*)ustar,yplus,delta,rrey,kappa
                      if (icont==tkmax) then
                      print*,"USTAR DOESNT CONVERGE"
                      print*,"Tang. vel. = ",vtan
                      stop
                      end if
                      end if
                      if (yplus<11.067d0) then
                      dom(ib)%tauwb(i,j)=rrey * vtan / delta
                      else
                      dom(ib)%tauwb(i,j)=ustar**2.0_dp             !tau_1
                      dom(ib)%tauwb2(i,j)=dom(ib)%tauwb(i,j)/vtan
                      !dom(ib)%vis(i,j,k)=dom(ib)%tauwb2(i,j)*delta
                      end if
                  end do
              end do


            CASE (6)
              delta=dom(ib)%dz
              n_x=0.0_dp
              n_y=0.0_dp
              n_z=-1.0_dp
              k=dom(ib)%kep
              do j=dom(ib)%jsp-1,dom(ib)%jep+1
                  do i=dom(ib)%isp-1,dom(ib)%iep+1
                      uc=0.5_dp*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                      vc=0.5_dp*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                      wc=0.5_dp*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                      vnor=n_x*uc+n_y*vc+n_z*wc
                      vtan=sqrt(abs(uc*uc+vc*vc+wc*wc-vnor*vnor+small))
                      if (dom(ib)%bc_top==62) then                !rough wall
                      ustar = kappa*vtan / log(30.d0*delta/fric)
                      yplus = delta*ustar/rrey
                      else if (dom(ib)%bc_top==61) then                !smooth wall
                      conv  = 1.0_dp
                      ustar = 1.0_dp
                      icont = 0
                      tkmax = 10
                      do while ((conv>(1.d-3)).and. &
                (icont<tkmax))
                          icont=icont+1
                          ustarold = ustar
                          yplus = max(1.0000001_dp,delta*ustar/rrey)
                          ustar = vtan*kappa/log(Ecte*yplus)
                          conv  = abs((ustar-ustarold)/ustar)
                      end do
                      if (icont==tkmax) then
                      print*,"USTAR DOESNT CONVERGE"
                      print*,"Tang. vel. = ",vtan
                      stop
                      end if
                      end if
                      if (yplus<11.63d0) then
                      dom(ib)%tauwt(i,j)=rrey * vtan / delta
                      else
                      dom(ib)%tauwt(i,j)=ustar**2.0_dp             !tau_1
                      dom(ib)%tauwt2(i,j)=dom(ib)%tauwt(i,j)/vtan
                      end if
                  end do
              end do

          end select

          return
      end subroutine log_law
!##########################################################################
