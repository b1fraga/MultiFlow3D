!=======================================================================
!               Scalar transport
!           Yan Liu
!           Bruño Fraga Bugallo
!           Cardiff 2015-2024
!=======================================================================
!##########################################################################
      subroutine sediment_init
!##########################################################################

          use vars
          use multiflow3d_mpi
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          integer :: ib,i,j,k,tti,ttj,ttk
          integer :: is,ie,js,je,ks,ke
          !real(dp) :: Vcell,temp1,temp2

          do ib=1,nbp

              tti=dom(ib)%ttc_i
              ttj=dom(ib)%ttc_j
              ttk=dom(ib)%ttc_k

              do k=1,ttk
                  do j=1,ttj
                      do i=1,tti

                          if (dom(ib)%z(k)<0.05_dp) then  !sludge
                          dom(ib)%S(i,j,k) = 1
                          else
                          dom(ib)%S(i,j,k) = 0  !fresh water
                          endif
                      enddo; enddo; enddo
          enddo

      end subroutine sediment_init

!##########################################################################
      subroutine sediment_4thtest
!##########################################################################
          use vars
          use multiflow3d_mpi
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          integer :: i,j,k,ib
          real(dp) :: dxx,dyy,dzz
          real(dp) :: conv,diff
          real(dp) :: duSdx,dvSdy,dwSdz,dwsSdz,dSdt
          real(dp) :: awS,aeS,asS,anS,ab_S,atS,apS
          real(dp) :: masout
          real(dp) :: ues,uws,vns,vss,wts,wbs,wsbs
          real(dp) :: masconv,masdiff,masws
          real(dp) :: kp,km,ku,kc,kd,b_r,ws,Vcell

          do ib=1,nbp

              Vcell = dom(ib)%dx*dom(ib)%dy*dom(ib)%dz

!     if (dom_id(ib).eq.9) then                         !release point for bubble plume
!       dom(ib)%So(dom(ib)%iep,dom(ib)%jep,30)=
!      &    dom(ib)%So(dom(ib)%iep,dom(ib)%jep,30) + (3.41d-6/Vcell)*dt
!       dom(ib)%So(dom(ib)%iep,dom(ib)%jep,31)=
!      &    dom(ib)%So(dom(ib)%iep,dom(ib)%jep,31) + (3.41d-6/Vcell)*dt
! !     write(6,*)'1',dom(ib)%S(dom(ib)%iep,dom(ib)%jep,29)
! !     write(6,*)dom(ib)%S(dom(ib)%iep,dom(ib)%jep,28)
! !     write(6,*)Vcell,dt
!     elseif (dom_id(ib).eq.10) then
!       dom(ib)%So(dom(ib)%isp,dom(ib)%jep,30)=
!      &    dom(ib)%So(dom(ib)%isp,dom(ib)%jep,30) + (3.41d-6/Vcell)*dt
!       dom(ib)%So(dom(ib)%isp,dom(ib)%jep,31)=
!      &    dom(ib)%So(dom(ib)%isp,dom(ib)%jep,31) + (3.41d-6/Vcell)*dt
! !     write(6,*)'2',dom(ib)%S(dom(ib)%isp,dom(ib)%jep,29)
! !     write(6,*)dom(ib)%S(dom(ib)%isp,dom(ib)%jep,28)
! !     write(6,*)Vcell,dt
!     elseif (dom_id(ib).eq.17) then
!       dom(ib)%So(dom(ib)%iep,dom(ib)%jsp,30)=
!      &    dom(ib)%So(dom(ib)%iep,dom(ib)%jsp,30) + (3.41d-6/Vcell)*dt
!       dom(ib)%So(dom(ib)%iep,dom(ib)%jsp,31)=
!      &    dom(ib)%So(dom(ib)%iep,dom(ib)%jsp,31) + (3.41d-6/Vcell)*dt
! !     write(6,*)'3',dom(ib)%S(dom(ib)%iep,dom(ib)%jsp,29)
! !     write(6,*)dom(ib)%S(dom(ib)%iep,dom(ib)%jsp,28)
! !     write(6,*)Vcell,dt
!     elseif (dom_id(ib).eq.18) then
!       dom(ib)%So(dom(ib)%isp,dom(ib)%jsp,30)=
!      &    dom(ib)%So(dom(ib)%isp,dom(ib)%jsp,30) + (3.41d-6/Vcell)*dt
!       dom(ib)%So(dom(ib)%isp,dom(ib)%jsp,31)=
!      &    dom(ib)%So(dom(ib)%isp,dom(ib)%jsp,31) + (3.41d-6/Vcell)*dt
! !     write(6,*)'4',dom(ib)%S(dom(ib)%isp,dom(ib)%jsp,29)
! !     write(6,*)dom(ib)%S(dom(ib)%isp,dom(ib)%jsp,28)
! !     write(6,*)Vcell,dt
!         endif


              if (itime == itime_start) then
              dom(ib)%sfactor = 1.0_dp
              end if

              dxx=dom(ib)%dx*dom(ib)%dx
              dyy=dom(ib)%dy*dom(ib)%dy
              dzz=dom(ib)%dz*dom(ib)%dz

              do k=dom(ib)%ksp,dom(ib)%kep
                  do i=dom(ib)%isp,dom(ib)%iep
                      do j=dom(ib)%jsp,dom(ib)%jep
!-------Convection
                          if(dom(ib)%u(i-1,j,k)>0.0_dp) then
                          ku=dom(ib)%So(i-2,j,k)
                          kc=dom(ib)%So(i-1,j,k)
                          kd=dom(ib)%So(i,j,k)
                          b_r=max(0.0_dp, &
                     min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                          km=(kc+0.5_dp*b_r*(kc-ku))
                          else if(dom(ib)%u(i-1,j,k)<0.0_dp) then
                          ku=dom(ib)%So(i+1,j,k)
                          kc=dom(ib)%So(i,j,k)
                          kd=dom(ib)%So(i-1,j,k)
                          b_r=max(0.0_dp, &
                     min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                          km=(kc+0.5_dp*b_r*(kc-ku))
                          else
                          km=0.5_dp*(dom(ib)%So(i,j,k)+dom(ib)%So(i-1,j,k))
                          end if

                          if(dom(ib)%u(i,j,k)>0.0_dp) then
                          ku=dom(ib)%So(i-1,j,k)
                          kc=dom(ib)%So(i,j,k)
                          kd=dom(ib)%So(i+1,j,k)
                          b_r=max(0.0_dp, &
                     min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                          kp=(kc+0.5_dp*b_r*(kc-ku))
                          else if(dom(ib)%u(i,j,k)<0.0_dp) then
                          ku=dom(ib)%So(i+2,j,k)
                          kc=dom(ib)%So(i+1,j,k)
                          kd=dom(ib)%So(i,j,k)
                          b_r=max(0.0_dp, &
                    min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                          kp=(kc+0.5_dp*b_r*(kc-ku))
                          else
                          kp=0.5_dp*(dom(ib)%So(i,j,k)+dom(ib)%So(i+1,j,k))
                          end if
                          duSdx=(dom(ib)%u(i,j,k)*kp-dom(ib)%u(i-1,j,k)*km)/dom(ib)%dx
!------
                          if(dom(ib)%v(i,j-1,k)>0.0_dp) then
                          ku=dom(ib)%So(i,j-2,k)
                          kc=dom(ib)%So(i,j-1,k)
                          kd=dom(ib)%So(i,j,k)
                          b_r=max(0.0_dp, &
                    min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                          km=(kc+0.5_dp*b_r*(kc-ku))
                          else if(dom(ib)%v(i,j-1,k)<0.0_dp) then
                          ku=dom(ib)%So(i,j+1,k)
                          kc=dom(ib)%So(i,j,k)
                          kd=dom(ib)%So(i,j-1,k)
                          b_r=max(0.0_dp, &
                    min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                          km=(kc+0.5_dp*b_r*(kc-ku))
                          else
                          km=0.5_dp*(dom(ib)%So(i,j,k)+dom(ib)%So(i,j-1,k))
                          end if

                          if(dom(ib)%v(i,j,k)>0.0_dp) then
                          ku=dom(ib)%So(i,j-1,k)
                          kc=dom(ib)%So(i,j,k)
                          kd=dom(ib)%So(i,j+1,k)
                          b_r=max(0.0_dp, &
                    min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                          kp=(kc+0.5_dp*b_r*(kc-ku))

                          else if(dom(ib)%v(i,j,k)<0.0_dp) then
                          ku=dom(ib)%So(i,j+2,k)
                          kc=dom(ib)%So(i,j+1,k)
                          kd=dom(ib)%So(i,j,k)
                          b_r=max(0.0_dp, &
                    min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                          kp=(kc+0.5_dp*b_r*(kc-ku))
                          else
                          kp=0.5_dp*(dom(ib)%So(i,j,k)+dom(ib)%So(i,j+1,k))
                          end if
                          dvSdy=(dom(ib)%v(i,j,k)*kp-dom(ib)%v(i,j-1,k)*km)/dom(ib)%dy
!-------
                          if(dom(ib)%w(i,j,k-1)>0.0_dp) then
                          ku=dom(ib)%So(i,j,k-2)
                          kc=dom(ib)%So(i,j,k-1)
                          kd=dom(ib)%So(i,j,k)
                          b_r=max(0.0_dp, &
                    min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                          km=(kc+0.5_dp*b_r*(kc-ku))

                          else if(dom(ib)%w(i,j,k-1)<0.0_dp) then
                          ku=dom(ib)%So(i,j,k+1)
                          kc=dom(ib)%So(i,j,k)
                          kd=dom(ib)%So(i,j,k-1)
                          b_r=max(0.0_dp, &
                    min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                          km=(kc+0.5_dp*b_r*(kc-ku))
                          else
                          km=0.5_dp*(dom(ib)%So(i,j,k)+dom(ib)%So(i,j,k-1))
                          end if
                          if(dom(ib)%w(i,j,k)>0.0_dp) then
                          ku=dom(ib)%So(i,j,k-1)
                          kc=dom(ib)%So(i,j,k)
                          kd=dom(ib)%So(i,j,k+1)
                          b_r=max(0.0_dp, &
                    min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                          kp=(kc+0.5_dp*b_r*(kc-ku))

                          else if(dom(ib)%w(i,j,k)<0.0_dp) then
                          ku=dom(ib)%So(i,j,k+2)
                          kc=dom(ib)%So(i,j,k+1)
                          kd=dom(ib)%So(i,j,k)
                          b_r=max(0.0_dp, &
                    min(2.0_dp*((kd-kc)/(kc-ku)),0.75_dp*((kd-kc)/(kc-ku))+0.25_dp,4.0_dp))
                          kp=(kc+0.5_dp*b_r*(kc-ku))
                          else
                          kp=0.5_dp*(dom(ib)%So(i,j,k)+dom(ib)%So(i,j,k+1))
                          end if
                          dwSdz=(dom(ib)%w(i,j,k)*kp-dom(ib)%w(i,j,k-1)*km)/dom(ib)%dz

                          conv=(duSdx+dvSdy+dwSdz)

!-------Diffusion
                          awS=-dom(ib)%vis(i,j,k)/(dxx*Pr)
                          aeS=-dom(ib)%vis(i,j,k)/(dxx*Pr)
                          anS=-dom(ib)%vis(i,j,k)/(dyy*Pr)
                          asS=-dom(ib)%vis(i,j,k)/(dyy*Pr)
                          atS=-dom(ib)%vis(i,j,k)/(dzz*Pr)
                          ab_S=-dom(ib)%vis(i,j,k)/(dzz*Pr)

                          apS = -1.0_dp*(awS+aeS+asS+anS+ab_S+atS)

                          diff=(apS*dom(ib)%So(i,j,k)+ &
                    anS*dom(ib)%So(i,j+1,k) + asS*dom(ib)%So(i,j-1,k)+ &
                    aeS*dom(ib)%So(i+1,j,k) + awS*dom(ib)%So(i-1,j,k)+ &
                    atS*dom(ib)%So(i,j,k+1) + ab_S*dom(ib)%So(i,j,k-1))

                          dom(ib)%S(i,j,k)=dom(ib)%So(i,j,k)-dt*(conv+diff)
                          dom(ib)%dens(i,j,k)=(0.007587_dp*dom(ib)%S(i,j,k)+0.9947_dp)*1000.0_dp

!   if (dom(ib)%S(i,j,k) .lt. 0.0_dp) then
!   write (81,*) dom(ib)%S(i,j,k), dom(ib)%So(i,j,k)
!   write (81,*) i,j,k
!   write (81,*) conv,diff
!   write (81,*) duSdx,dvSdy,dwSdz
!   end if

                      end do
                  end do
              end do

          end do

          do ib =1, nbp
              do k=dom(ib)%ksp,dom(ib)%kep
                  do i=dom(ib)%isp,dom(ib)%iep
                      do j=dom(ib)%jsp,dom(ib)%jep
                          if (dom(ib)%S(i,j,k) > 100) then
!   call tecplot_S(itime)
                          write(6,*)'ERROR: scalar too big'
                          stop
                          end if
                      end do
                  end do
              end do
          end do

          call exchange(8)  !S
          call exchange(10)  !dens

          call boundS

          return
      end subroutine sediment_4thtest
!##########################################################################
      subroutine boundS
!##########################################################################
          use multiflow3d_mpi
          use vars
          use multidata
          implicit none
          integer :: i,j,k,ib,ni,nj,nk,ly
          integer :: is,ie,js,je,ks,ke
          real(dp) :: absz,absy

          if (PERIODIC) call exchange_bc(8,pl_ex)

          do ly=0,pl_ex

              do ib=1,nbp
                  ni=dom(ib)%ttc_i; nj=dom(ib)%ttc_j; nk=dom(ib)%ttc_k
                  is=dom(ib)%isp; ie=dom(ib)%iep
                  js=dom(ib)%jsp; je=dom(ib)%jep
                  ks=dom(ib)%ksp; ke=dom(ib)%kep

! Boundary Conditions for S - ASSUMING ALL ADIABATIC BCs
!..............................................................................
!=== West ===>   ..  4=wall  ..    1=Inflow
!..............................................................................
                  if (dom(ib)%iprev<0) then
                  !      if (dom(ib)%Tbc_west.eq.4) then
                  !   do k=ks-1,ke+1; do j=js-1,je+1
                  do k=1,nk; do j=1,nj
                          dom(ib)%S(is-1-ly,j,k)= dom(ib)%S(is+ly,j,k)
                          dom(ib)%dens(is-1-ly,j,k)=dom(ib)%dens(is+ly,j,k)
                      end do; end do

                  !      else if (dom(ib)%Tbc_west.eq.1) then                   !CHANGE
                  !         do k=ks-1,ke+1; do j=js-1,je+1
                  !            dom(ib)%S(is-1-ly,j,k)= 1.0_dp
                  !         end do; end do
                  !    else if (dom(ib)%Tbc_west.eq. 11) then
                  !       do k=ks-1,ke+1; do j=js-1,je+1
                  !       absz=abs(dom(ib)%zc(k)-5.22)
                  !       absy=abs(dom(ib)%yc(j)-4.32)
                  !       if (absz .le. 0.5_dp*dom(ib)%dz) then
                  !         if (absy .le. 0.5_dp*dom(ib)%dy) then
                  !        dom(ib)%S(is-1-ly,j,k)= 1.0_dp
                  !         else
                  !            dom(ib)%S(is-1-ly,j,k)= 0.0_dp
                  !         end if
                  !       else
                  !          dom(ib)%S(is-1-ly,j,k)= 0.0_dp
                  !       end if
                  !       end do; end do
                  !      end if
                  end if
!...............................................................................
!=== East ===>   ..  4=wall  ..    2=Outflow
!...............................................................................
                  if (dom(ib)%inext<0) then
                  !      if (dom(ib)%Tbc_east.eq.4) then
                  !   do k=ks-1,ke+1; do j=js-1,je+1
                  do k=1,nk; do j=1,nj
                          dom(ib)%S(ie+1+ly,j,k)= dom(ib)%S(ie-ly,j,k)
                          dom(ib)%dens(ie+1+ly,j,k)= dom(ib)%dens(ie-ly,j,k)
                      end do; end do

                  !      else if (dom(ib)%Tbc_east.eq.2_dp) then
                  !         do k=ks-1,ke+1; do j=js-1,je+1
                  !            dom(ib)%S(ie+1+ly,j,k)= dom(ib)%S(ie-ly,j,k)
                  !         end do; end do
                  !      end if
                  end if
!...............................................................................
!=== South ===>  ..  4=wall  ..
!...............................................................................
                  if (dom(ib)%jprev<0) then
                  !      else if (dom(ib)%Tbc_east.eq.2) then
                  !   do k=ks-1,ke+1; do i=is-1,ie+1
                  do k=1,nk; do i=1,ni
                          dom(ib)%S(i,js-1-ly,k)= dom(ib)%S(i,js+ly,k)
                          dom(ib)%dens(i,js-1-ly,k)= dom(ib)%dens(i,js+ly,k)
                      end do; end do

                  !      end if
                  end if
!.............................................................................
!=== North ===>  ..  4=wall  ..
!.............................................................................
                  if (dom(ib)%jnext<0) then
                  !      if (dom(ib)%Tbc_north.eq.4) then
                  !   do k=ks-1,ke+1; do i=is-1,ie+1
                  do k=1,nk; do i=1,ni
                          dom(ib)%S(i,je+1+ly,k) = dom(ib)%S(i,je-ly,k)
                          dom(ib)%dens(i,je+1+ly,k) = dom(ib)%dens(i,je-ly,k)
                      end do; end do

                  !      end if
                  end if
!...............................................................................
!=== Bottom ===>  ..  6=Net deposition  ..   7=Erosion
!...............................................................................
                  if (dom(ib)%kprev<0) then
                  !      if (dom(ib)%Tbc_bottom.eq.6) then
                  !         do j=js-1,je+1; do i=is-1,ie+1
                  !            dom(ib)%S(i,j,ks-1-ly)= 0.0_dp
                  !         end do; end do

                  !      else if (dom(ib)%Tbc_bottom.eq.7) then
                  !   do j=js-1,je+1; do i=is-1,ie+1
                  do j=1,nj; do i=1,ni
                          dom(ib)%S(i,j,ks-1-ly)= dom(ib)%S(i,j,ks+ly)
                          dom(ib)%dens(i,j,ks-1-ly)= dom(ib)%dens(i,j,ks+ly)
                      end do; end do
                  !      end if
                  end if
!.............................................................................
!=== Top ===>  ..  8=free surface
!.............................................................................
                  if (dom(ib)%knext<0) then
                  !      if (dom(ib)%Tbc_top.eq.8) then
                  !   do j=js-1,je+1; do i=is-1,ie+1
                  do j=1,nj; do i=1,ni
                          dom(ib)%S(i,j,ke+1+ly) = dom(ib)%S(i,j,ke-ly)
                          dom(ib)%dens(i,j,ke+1+ly) = dom(ib)%dens(i,j,ke-ly)
                      end do; end do

                  !      end if
                  end if

!==============================================================================
              end do

          end do  ! ly

      end subroutine boundS

!##########################################################################
      subroutine NonNewtonian
!##########################################################################

          use vars
          use multiflow3d_mpi
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none

          integer :: ib,i,j,k
          real(dp) :: strain,n
          !constitutive relationship is for sludge
          do ib=1,nbp
              do i=dom(ib)%isp,dom(ib)%iep
                  do j=dom(ib)%jsp,dom(ib)%jep
                      do k=dom(ib)%ksp,dom(ib)%kep
                          n=0.6894_dp+0.0046831_dp*(dom(ib)%T(i,j,j)-273) &
                    -0.042813_dp*dom(ib)%S(i,j,k)
                          if (strain(i,j,k)>1d-8) then
                          dom(ib)%mu(i,j,k)= rrey*dens*strain(i,j,k)**(n-1.d0)
                          else
                          dom(ib)%mu(i,j,k)= rrey*dens
                          endif
                          dom(ib)%vis(i,j,k)=dom(ib)%mu(i,j,k)/dom(ib)%dens(i,j,k)

                      enddo;enddo;enddo
          enddo

      end subroutine NonNewtonian

!##########################################################################
      subroutine Active_scalar
!##########################################################################
          use vars
          use multiflow3d_mpi
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none

          integer :: ib,i,j,k
          !constitutive relationship is for density
          do ib=1,nbp
              do i=1,dom(ib)%ttc_i;do j=1,dom(ib)%ttc_j;do k=1,dom(ib)%ttc_k
                          dom(ib)%dens(i,j,k)=0.0367_dp*dom(ib)%S(i,j,k)**3.d0 &       !based on sludge on this case
                                -2.38_dp*dom(ib)%S(i,j,k)**2.d0 &
                                +14.6_dp*dom(ib)%S(i,j,k)+1000

                      enddo;enddo;enddo
          enddo

      end subroutine Active_scalar
