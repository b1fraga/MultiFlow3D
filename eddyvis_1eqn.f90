!##########################################################################
      subroutine eddyv_1eqn
!##########################################################################
          use vars
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          integer :: i,j,k,rk
          integer :: ib,is,ie,js,je,ks,ke
          double precision :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
          double precision :: dukdx,dvkdy,dwkdz,kp,km,ku,kc,kd,b_r
          double precision :: awT,aeT,asT,anT,abT,atT,apT
          double precision :: dxx,dyy,dzz,vsgs
          double precision :: visc_w,visc_e,visc_s,visc_n,visc_b,visc_t
          double precision :: vr_a,vr_b,vr_c,vr_d
          double precision :: s12,s13,s23
          double precision :: ss,cons_k,cons_eps
          double precision :: conv,diff,prod,other
          double precision :: delta_grid
          double precision :: alfark(3)
          double precision ::dx,dy,dz
          integer :: sn
          character*8 :: chb1
          character*25 :: gf

          alfark(1)=1.0_dp/3.
          alfark(2)=0.5_dp
          alfark(3)=1.0_dp

          cons_k=0.05; cons_eps=1.00

          do rk=1,3

              do ib=1,nbp

                  do k=1,dom(ib)%ttc_k
                      do i=1,dom(ib)%ttc_i
                          do j=1,dom(ib)%ttc_j
                              dom(ib)%ksgso(i,j,k)=dom(ib)%ksgs(i,j,k)
                          end do
                      end do
                  end do

                  dxx=dom(ib)%dx*dom(ib)%dx
                  dyy=dom(ib)%dy*dom(ib)%dy
                  dzz=dom(ib)%dz*dom(ib)%dz

                  delta_grid=(dom(ib)%dx*dom(ib)%dy*dom(ib)%dz)**(1.0_dp/3.0_dp)

                  is=dom(ib)%isp; ie=dom(ib)%iep
                  js=dom(ib)%jsp; je=dom(ib)%jep
                  ks=dom(ib)%ksp; ke=dom(ib)%kep

                  do k=ks-1,ke+1
                      do j=js-1,je+1
                          do i=is-1,ie+1

                              if (L_LSM)  rrey=dom(ib)%mu(i,j,k)/dom(ib)%dens(i,j,k)

!====================================================
                              vr_a = 0.25_dp*( dom(ib)%u(i,j,k)   + dom(ib)%u(i-1,j,k) + &
                        dom(ib)%u(i,j+1,k) + dom(ib)%u(i-1,j+1,k) )
                              vr_b = 0.25_dp*( dom(ib)%u(i,j,k)   + dom(ib)%u(i-1,j,k) + &
                        dom(ib)%u(i,j-1,k) + dom(ib)%u(i-1,j-1,k) )

                              vr_c = 0.25_dp*( dom(ib)%u(i,j,k)   + dom(ib)%u(i-1,j,k) + &
                        dom(ib)%u(i,j,k+1) + dom(ib)%u(i-1,j,k+1) )
                              vr_d = 0.25_dp*( dom(ib)%u(i,j,k)   + dom(ib)%u(i-1,j,k) + &
                        dom(ib)%u(i,j,k-1) + dom(ib)%u(i-1,j,k-1) )

                              dudx = ( dom(ib)%u(i,j,k) - dom(ib)%u(i-1,j,k) )/dom(ib)%dx
                              dudy = ( vr_a - vr_b )/dom(ib)%dy
                              dudz = ( vr_c - vr_d )/dom(ib)%dz
!====================================================
                              vr_a = 0.25_dp*( dom(ib)%v(i,j,k)   + dom(ib)%v(i,j-1,k) + &
                        dom(ib)%v(i+1,j,k) + dom(ib)%v(i+1,j-1,k) )
                              vr_b = 0.25_dp*( dom(ib)%v(i,j,k)   + dom(ib)%v(i,j-1,k) + &
                        dom(ib)%v(i-1,j,k) + dom(ib)%v(i-1,j-1,k) )

                              vr_c = 0.25_dp*( dom(ib)%v(i,j,k)   + dom(ib)%v(i,j-1,k) + &
                        dom(ib)%v(i,j,k+1) + dom(ib)%v(i,j-1,k+1) )
                              vr_d = 0.25_dp*( dom(ib)%v(i,j,k)   + dom(ib)%v(i,j-1,k) + &
                        dom(ib)%v(i,j,k-1) + dom(ib)%v(i,j-1,k-1) )

                              dvdy = ( dom(ib)%v(i,j,k) - dom(ib)%v(i,j-1,k) )/dom(ib)%dy
                              dvdx = ( vr_a - vr_b )/dom(ib)%dx
                              dvdz = ( vr_c - vr_d )/dom(ib)%dz
!====================================================
                              vr_a = 0.25_dp*( dom(ib)%w(i,j,k)   + dom(ib)%w(i,j,k-1) + &
                        dom(ib)%w(i+1,j,k) + dom(ib)%w(i+1,j,k-1) )
                              vr_b = 0.25_dp*( dom(ib)%w(i,j,k)   + dom(ib)%w(i,j,k-1) + &
                        dom(ib)%w(i-1,j,k) + dom(ib)%w(i-1,j,k-1) )

                              vr_c = 0.25_dp*( dom(ib)%w(i,j,k)   + dom(ib)%w(i,j,k-1) + &
                        dom(ib)%w(i,j+1,k) + dom(ib)%w(i,j+1,k-1) )
                              vr_d = 0.25_dp*( dom(ib)%w(i,j,k)   + dom(ib)%w(i,j,k-1) + &
                        dom(ib)%w(i,j-1,k) + dom(ib)%w(i,j-1,k-1) )

                              dwdz = ( dom(ib)%w(i,j,k) - dom(ib)%w(i,j,k-1) )/dom(ib)%dz
                              dwdx = ( vr_a - vr_b )/dom(ib)%dx
                              dwdy = ( vr_c - vr_d )/dom(ib)%dy


                              s12 = 0.5_dp * (dudy + dvdx)
                              s13 = 0.5_dp * (dudz + dwdx)
                              s23 = 0.5_dp * (dvdz + dwdy)
                              ss  = ( dudx*dudx + dvdy*dvdy   + dwdz*dwdz  + &
                            2.0_dp*s12*s12   + 2.0_dp*s13*s13 + 2.0_dp*s23*s23 )

                              vsgs=cons_k*delta_grid*sqrt(dom(ib)%ksgso(i,j,k))

                              prod=2.0_dp*vsgs*ss

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                              if(dom(ib)%u(i-1,j,k)>0.0_dp) then
                              ku=dom(ib)%ksgso(i-2,j,k)
                              kc=dom(ib)%ksgso(i-1,j,k)
                              kd=dom(ib)%ksgso(i,j,k)
                              b_r=max(0.0_dp, &
                        min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                              km=kc+0.5_dp*b_r*(kc-ku)
                              else if(dom(ib)%u(i-1,j,k)<0.0_dp) then
                              ku=dom(ib)%ksgso(i+1,j,k)
                              kc=dom(ib)%ksgso(i,j,k)
                              kd=dom(ib)%ksgso(i-1,j,k)
                              b_r=max(0.0_dp, &
                        min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                              km=kc+0.5_dp*b_r*(kc-ku)
                              else
                              km=0.5_dp*(dom(ib)%ksgso(i,j,k)+dom(ib)%ksgso(i-1,j,k))
                              end if
                              if(dom(ib)%u(i,j,k)>0.0_dp) then
                              ku=dom(ib)%ksgso(i-1,j,k)
                              kc=dom(ib)%ksgso(i,j,k)
                              kd=dom(ib)%ksgso(i+1,j,k)
                              b_r=max(0.0_dp, &
                        min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                              kp=kc+0.5_dp*b_r*(kc-ku)
                              else if(dom(ib)%u(i,j,k)<0.0_dp) then
                              ku=dom(ib)%ksgso(i+2,j,k)
                              kc=dom(ib)%ksgso(i+1,j,k)
                              kd=dom(ib)%ksgso(i,j,k)
                              b_r=max(0.0_dp, &
                        min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                              kp=kc+0.5_dp*b_r*(kc-ku)
                              else
                              kp=0.5_dp*(dom(ib)%ksgso(i,j,k)+dom(ib)%ksgso(i+1,j,k))
                              end if
                              dukdx=(dom(ib)%u(i,j,k)*kp-dom(ib)%u(i-1,j,k)*km)/dom(ib)%dx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              if(dom(ib)%v(i,j-1,k)>0.0_dp) then
                              ku=dom(ib)%ksgso(i,j-2,k)
                              kc=dom(ib)%ksgso(i,j-1,k)
                              kd=dom(ib)%ksgso(i,j,k)
                              b_r=max(0.0_dp, &
                        min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                              km=kc+0.5_dp*b_r*(kc-ku)
                              else if(dom(ib)%v(i,j-1,k)<0.0_dp) then
                              ku=dom(ib)%ksgso(i,j+1,k)
                              kc=dom(ib)%ksgso(i,j,k)
                              kd=dom(ib)%ksgso(i,j-1,k)
                              b_r=max(0.0_dp, &
                        min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                              km=kc+0.5_dp*b_r*(kc-ku)
                              else
                              km=0.5_dp*(dom(ib)%ksgso(i,j,k)+dom(ib)%ksgso(i,j-1,k))
                              end if
                              if(dom(ib)%v(i,j,k)>0.0_dp) then
                              ku=dom(ib)%ksgso(i,j-1,k)
                              kc=dom(ib)%ksgso(i,j,k)
                              kd=dom(ib)%ksgso(i,j+1,k)
                              b_r=max(0.0_dp, &
                        min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                              kp=kc+0.5_dp*b_r*(kc-ku)
                              else if(dom(ib)%v(i,j,k)<0.0_dp) then
                              ku=dom(ib)%ksgso(i,j+2,k)
                              kc=dom(ib)%ksgso(i,j+1,k)
                              kd=dom(ib)%ksgso(i,j,k)
                              b_r=max(0.0_dp, &
                        min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                              kp=kc+0.5_dp*b_r*(kc-ku)
                              else
                              kp=0.5_dp*(dom(ib)%ksgso(i,j,k)+dom(ib)%ksgso(i,j+1,k))
                              end if
                              dvkdy=(dom(ib)%v(i,j,k)*kp-dom(ib)%v(i,j-1,k)*km)/dom(ib)%dy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              if(dom(ib)%w(i,j,k-1)>0.0_dp) then
                              ku=dom(ib)%ksgso(i,j,k-2)
                              kc=dom(ib)%ksgso(i,j,k-1)
                              kd=dom(ib)%ksgso(i,j,k)
                              b_r=max(0.0_dp, &
                        min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                              km=kc+0.5_dp*b_r*(kc-ku)
                              else if(dom(ib)%w(i,j,k-1)<0.0_dp) then
                              ku=dom(ib)%ksgso(i,j,k+1)
                              kc=dom(ib)%ksgso(i,j,k)
                              kd=dom(ib)%ksgso(i,j,k-1)
                              b_r=max(0.0_dp, &
                        min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                              km=kc+0.5_dp*b_r*(kc-ku)
                              else
                              km=0.5_dp*(dom(ib)%ksgso(i,j,k)+dom(ib)%ksgso(i,j,k-1))
                              end if
                              if(dom(ib)%w(i,j,k)>0.0_dp) then
                              ku=dom(ib)%ksgso(i,j,k-1)
                              kc=dom(ib)%ksgso(i,j,k)
                              kd=dom(ib)%ksgso(i,j,k+1)
                              b_r=max(0.0_dp, &
                        min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                              kp=kc+0.5_dp*b_r*(kc-ku)
                              else if(dom(ib)%w(i,j,k)<0.0_dp) then
                              ku=dom(ib)%ksgso(i,j,k+2)
                              kc=dom(ib)%ksgso(i,j,k+1)
                              kd=dom(ib)%ksgso(i,j,k)
                              b_r=max(0.0_dp, &
                        min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                              kp=kc+0.5_dp*b_r*(kc-ku)
                              else
                              kp=0.5_dp*(dom(ib)%ksgso(i,j,k)+dom(ib)%ksgso(i,j,k+1))
                              end if
                              dwkdz=(dom(ib)%w(i,j,k)*kp-dom(ib)%w(i,j,k-1)*km)/dom(ib)%dz
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

                              conv=(dukdx+dvkdy+dwkdz)


                              visc_w=rrey+0.5_dp*cons_k*delta_grid* &
                        (sqrt(dom(ib)%ksgso(i,j,k))+sqrt(dom(ib)%ksgso(i-1,j,k)))
                              visc_e=rrey+0.5_dp*cons_k*delta_grid* &
                        (sqrt(dom(ib)%ksgso(i,j,k))+sqrt(dom(ib)%ksgso(i+1,j,k)))
                              visc_s=rrey+0.5_dp*cons_k*delta_grid* &
                        (sqrt(dom(ib)%ksgso(i,j,k))+sqrt(dom(ib)%ksgso(i,j-1,k)))
                              visc_n=rrey+0.5_dp*cons_k*delta_grid* &
                        (sqrt(dom(ib)%ksgso(i,j,k))+sqrt(dom(ib)%ksgso(i,j+1,k)))
                              visc_b=rrey+0.5_dp*cons_k*delta_grid* &
                        (sqrt(dom(ib)%ksgso(i,j,k))+sqrt(dom(ib)%ksgso(i,j,k-1)))
                              visc_t=rrey+0.5_dp*cons_k*delta_grid* &
                        (sqrt(dom(ib)%ksgso(i,j,k))+sqrt(dom(ib)%ksgso(i,j,k+1)))

                              awT=visc_w/dxx; aeT=visc_e/dxx
                              anT=visc_n/dyy; asT=visc_s/dyy
                              atT=visc_t/dzz; abT=visc_b/dzz
                              apT = -1.0_dp*(awT+aeT+asT+anT+abT+atT)
                              diff=(apT*dom(ib)%ksgso(i,j,k)+ &
                        anT*dom(ib)%ksgso(i,j+1,k) + asT*dom(ib)%ksgso(i,j-1,k)+ &
                        aeT*dom(ib)%ksgso(i+1,j,k) + awT*dom(ib)%ksgso(i-1,j,k)+ &
                        atT*dom(ib)%ksgso(i,j,k+1) + abT*dom(ib)%ksgso(i,j,k-1))

                              other=cons_eps*sqrt(dom(ib)%ksgso(i,j,k)**3)/delta_grid

                              dom(ib)%ksgs(i,j,k)=(dom(ib)%ksgso(i,j,k)+ &
                        alfark(rk)*dt*(diff-conv+prod-other))

                              if(dom(ib)%ksgs(i,j,k)<0.0_dp) then
                              print*,'ERRORRRR in the 1-EQN model & STOP'
                              write (6,*) 'ERRORRRR in the 1-EQN model & STOP'
                              stop
                              end if

                              dom(ib)%vis(i,j,k) = rrey + &
                        cons_k*delta_grid*sqrt(dom(ib)%ksgs(i,j,k))

                          end do
                      end do
                  end do

              end do

              call boundksgs_1eqn
              call exchange(6)

          end do


          do ib=1,nbp
              is=dom(ib)%isp; ie=dom(ib)%iep
              js=dom(ib)%jsp; je=dom(ib)%jep
              ks=dom(ib)%ksp; ke=dom(ib)%kep

              dx=dom(ib)%dx
              dy=dom(ib)%dy
              dz=dom(ib)%dz
!..........................................................................
!=== West ===> ..  4=wall  ..   1=Inlet
!..........................................................................
              if (dom(ib)%iprev<0) then
              if (dom(ib)%bc_west==4) then
              do k=ks-1,ke+1
                  do j=js-1,je+1
                      dom(ib)%vis(is,j,k)= rrey
                      dom(ib)%vis(is-1,j,k)= dom(ib)%vis(is,j,k)
                  end do
              end do
              else if (dom(ib)%bc_west==1) then
              do k=ks-1,ke+1
                  do j=js-1,je+1
                      dom(ib)%vis(is-1,j,k)= dom(ib)%vis(is,j,k)
                  end do
              end do
              else
              do k=ks-1,ke+1
                  do j=js-1,je+1
                      dom(ib)%vis(is-1,j,k)= dom(ib)%vis(is,j,k)
                  end do
              end do
              end if
              else
              do k=ks-1,ke+1
                  do j=js-1,je+1
                      dom(ib)%vis(is-1,j,k)= dom(ib)%vis(is,j,k)
                  end do
              end do
              end if
!..........................................................................
!=== East ===> ..  4=wall  ..   2=Outflow
!..........................................................................
              if (dom(ib)%inext<0) then
              if (dom(ib)%bc_east==4) then
              do k=ks-1,ke+1
                  do j=js-1,je+1
                      dom(ib)%vis(ie,j,k)= rrey
                      dom(ib)%vis(ie+1,j,k)= dom(ib)%vis(ie,j,k)
                  end do
              end do
              else if (dom(ib)%bc_east==2) then
              do k=ks-1,ke+1
                  do j=js-1,je+1
                      dom(ib)%vis(ie+1,j,k)= dom(ib)%vis(ie,j,k)
                  end do
              end do
              else
              do k=ks-1,ke+1
                  do j=js-1,je+1
                      dom(ib)%vis(ie+1,j,k)= dom(ib)%vis(ie,j,k)
                  end do
              end do
              end if
              else
              do k=ks-1,ke+1
                  do j=js-1,je+1
                      dom(ib)%vis(ie+1,j,k)= dom(ib)%vis(ie,j,k)
                  end do
              end do
              end if
!..........................................................................
!=== South ===> ..   4=wall ..   3=Symmetry
!..........................................................................
              if (dom(ib)%jprev<0) then
              if (dom(ib)%bc_south==4) then
              do k=ks-1,ke+1
                  do i=is-1,ie+1
                      dom(ib)%vis(i,js,k)= rrey
                      dom(ib)%vis(i,js-1,k)= dom(ib)%vis(i,js,k)
                  end do
              end do
              else
              do k=ks-1,ke+1
                  do i=is-1,ie+1
                      dom(ib)%vis(i,js-1,k)= dom(ib)%vis(i,js,k)
                  end do
              end do
              end if
              else
              do k=ks-1,ke+1
                  do i=is-1,ie+1
                      dom(ib)%vis(i,js-1,k)= dom(ib)%vis(i,js,k)
                  end do
              end do
              end if
!..........................................................................
!=== North ===>  ..   4=wall ..   44=moving wall ..  3=Symmetry
!..........................................................................
              if (dom(ib)%jnext<0) then
              if (dom(ib)%bc_north==4) then
              do k=ks-1,ke+1
                  do i=is-1,ie+1
                      dom(ib)%vis(i,je,k) = rrey
                      dom(ib)%vis(i,je+1,k) =  dom(ib)%vis(i,je,k)
                  end do
              end do
              else
              do k=ks-1,ke+1
                  do i=is-1,ie+1
                      dom(ib)%vis(i,je+1,k) =  dom(ib)%vis(i,je,k)
                  end do
              end do
              end if
              else
              do k=ks-1,ke+1
                  do i=is-1,ie+1
                      dom(ib)%vis(i,je+1,k) =  dom(ib)%vis(i,je,k)
                  end do
              end do
              end if
!..........................................................................
!=== Bottom ===> ..   4=wall ..   3=Symmetry
!..........................................................................
              if (dom(ib)%kprev<0) then
              if (dom(ib)%bc_bottom==4) then
              do j=js-1,je+1
                  do i=is-1,ie+1
                      dom(ib)%vis(i,j,ks)= rrey
                      dom(ib)%vis(i,j,ks-1)= dom(ib)%vis(i,j,ks)
                  end do
              end do
              else
              do j=js-1,je+1
                  do i=is-1,ie+1
                      dom(ib)%vis(i,j,ks-1)= dom(ib)%vis(i,j,ks)
                  end do
              end do
              end if
              else
              do j=js-1,je+1
                  do i=is-1,ie+1
                      dom(ib)%vis(i,j,ks-1)= dom(ib)%vis(i,j,ks)
                  end do
              end do
              end if
!..........................................................................
!=== Top ===>  ..   4=wall ..     3=Symmetry
!..........................................................................
              if (dom(ib)%knext<0) then
              if (dom(ib)%bc_top==4) then
              do j=js-1,je+1
                  do i=is-1,ie+1
                      dom(ib)%vis(i,j,ke)   = rrey
                      dom(ib)%vis(i,j,ke+1)   = dom(ib)%vis(i,j,ke)
                  end do
              end do
              else
              do j=js-1,je+1
                  do i=is-1,ie+1
                      dom(ib)%vis(i,j,ke+1)   = dom(ib)%vis(i,j,ke)
                  end do
              end do
              end if
              else
              do j=js-1,je+1
                  do i=is-1,ie+1
                      dom(ib)%vis(i,j,ke+1)   = dom(ib)%vis(i,j,ke)
                  end do
              end do
              end if

          end do

          return
      end subroutine eddyv_1eqn
!##########################################################################
