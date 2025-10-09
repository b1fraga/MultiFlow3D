!##########################################################################
      subroutine wall_function(bound,ib)
!           Bruño Fraga Bugallo
!           Cardiff 2014
!           werner/wengle type boundary conditions
!##########################################################################
          use vars
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64

          implicit none
          integer i,j,k,ib,bound,cond
          double precision delta,n_x,n_y,n_z,vnor,vtan,dvtan,sub
          double precision uc,vc,wc,small,dycell,rycell,vtankr
          double precision tausub,taupow
          double precision aaa,bbb,const1,const2,const3,const4


          if (LAS.or.L_LSM) then                                            !variable density
          do i=1,dom(ib)%ttc_i ;do j=1,dom(ib)%ttc_j;do k=1,dom(ib)%ttc_k
                      rrey=dom(ib)%mu(i,j,k)/dom(ib)%dens(i,j,k)
                  enddo;enddo;enddo
          endif
          small = 1.e-30_dp

          SELECT CASE (bound)
            CASE (1)
              cond=dom(ib)%bc_west
            CASE (2)
              cond=dom(ib)%bc_east
            CASE (3)
              cond=dom(ib)%bc_south
            CASE (4)
              cond=dom(ib)%bc_north
            CASE (5)
              cond=dom(ib)%bc_bottom
            CASE (6)
              cond=dom(ib)%bc_top
          end select

SELECT CASE (cond)

  CASE (63)
!.....specify constants for 1/6 power law ..............................
    aaa = 8.3_dp
    bbb = 0.1666666666_dp

  CASE (64)
!.....specify constants for 1/7 power law ..............................
    aaa = 8.3_dp
    bbb = 0.1428571429_dp

  CASE (65)
!.....specify constants for 1/8 power law ..............................
    aaa = 8.3_dp
    bbb = 0.125_dp

END SELECT


!.....constant factors .................................................

          const1 = 0.5_dp * (1.0_dp - bbb) * aaa ** ((1.0_dp + bbb) / (1.0_dp - bbb))
          const2 = (1.0_dp + bbb) / aaa
          const3 = aaa ** (2.0_dp / (1.0_dp - bbb))
          const4 = 2.0_dp / (1.0_dp + bbb)

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

                      dycell = 2.0_dp * delta
                      rycell = 1.0_dp / dycell
                      vtankr = 0.5_dp * rrey * rycell * const3
                      dvtan  = vtankr - vtan
                      sub    = MAX (SIGN(1.0_dp,dvtan),0.0_dp)

                      tausub   = rrey * vtan / delta
                      taupow   = ( const1 * (rrey * rycell)**(1.0_dp+bbb) + &
                           ( const2 * (rrey * rycell)**bbb) * vtan) &
                              ** const4
                      dom(ib)%tauww(j,k)=(sub*tausub+(1.0_dp-sub)*taupow)   !tau_1
                      dom(ib)%tauww2(j,k)=dom(ib)%tauww(j,k)/vtan       !needs to be multiplied by a velocity component to provide tau_1j
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

                      dycell = 2.0_dp * delta
                      rycell = 1.0_dp / dycell
                      vtankr = 0.5_dp * rrey * rycell * const3
                      dvtan  = vtankr - vtan
                      sub    = MAX (SIGN(1.0_dp,dvtan),0.0_dp)

                      tausub   = rrey * vtan / delta
                      taupow   = ( const1 * (rrey * rycell)**(1.0_dp+bbb) + &
                           ( const2 * (rrey * rycell)**bbb) * vtan) &
                              ** const4
                      dom(ib)%tauwe(j,k)=(sub*tausub+(1.0_dp-sub)*taupow)   !units m2/s2
                      dom(ib)%tauwe2(j,k)=dom(ib)%tauwe(j,k)/vtan       !units m/s
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

                      dycell = 2.0_dp * delta
                      rycell = 1.0_dp / dycell
                      vtankr = 0.5_dp * rrey * rycell * const3
                      dvtan  = vtankr - vtan
                      sub    = MAX (SIGN(1.0_dp,dvtan),0.0_dp)

                      tausub   = rrey * vtan / delta
                      taupow   = ( const1 * (rrey * rycell)**(1.0_dp+bbb) + &
                           ( const2 * (rrey * rycell)**bbb) * vtan) &
                              ** const4
                      dom(ib)%tauws(i,k)=(sub*tausub+(1.0_dp-sub)*taupow)
                      dom(ib)%tauws2(i,k)=dom(ib)%tauws(i,k)/vtan
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

                      dycell = 2.0_dp * delta
                      rycell = 1.0_dp / dycell
                      vtankr = 0.5_dp * rrey * rycell * const3
                      dvtan  = vtankr - vtan
                      sub    = MAX (SIGN(1.0_dp,dvtan),0.0_dp)

                      tausub   = rrey * vtan / delta
                      taupow   = ( const1 * (rrey * rycell)**(1.0_dp+bbb) + &
                           ( const2 * (rrey * rycell)**bbb) * vtan) &
                              ** const4
                      dom(ib)%tauwn(i,k)=(sub*tausub+(1.0_dp-sub)*taupow)
                      dom(ib)%tauwn2(i,k)=dom(ib)%tauwn(i,k)/vtan
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

                      dycell = 2.0_dp * delta
                      rycell = 1.0_dp / dycell
                      vtankr = 0.5_dp * rrey * rycell * const3
                      dvtan  = vtankr - vtan
                      sub    = MAX (SIGN(1.0_dp,dvtan),0.0_dp)

                      tausub   = rrey * vtan / delta
                      taupow   = ( const1 * (rrey * rycell)**(1.0_dp+bbb) + &
                           ( const2 * (rrey * rycell)**bbb) * vtan) &
                              ** const4
                      dom(ib)%tauwb(i,j)=(sub*tausub+(1.0_dp-sub)*taupow)
                      dom(ib)%tauwb2(i,j)=dom(ib)%tauwb(i,j)/vtan
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

                      dycell = 2.0_dp * delta
                      rycell = 1.0_dp / dycell
                      vtankr = 0.5_dp * rrey * rycell * const3
                      dvtan  = vtankr - vtan
                      sub    = MAX (SIGN(1.0_dp,dvtan),0.0_dp)

                      tausub   = rrey * vtan / delta
                      taupow   = ( const1 * (rrey * rycell)**(1.0_dp+bbb) + &
                           ( const2 * (rrey * rycell)**bbb) * vtan) &
                              ** const4
                      dom(ib)%tauwt(i,j)=(sub*tausub+(1.0_dp-sub)*taupow)
                      dom(ib)%tauwt2(i,j)=dom(ib)%tauwt(i,j)/vtan
                  end do
              end do

          end select

          return
      end subroutine wall_function
!##########################################################################
