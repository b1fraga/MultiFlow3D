!#############################################################################
      subroutine boundksgs(cmu)
!           Bruño Fraga Bugallo
!           Cardiff 2016
!##############################################################################
          use vars
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          integer :: i,j,k,ib,ly
          integer :: is,ie,js,je,ks,ke
          real(dp) :: cmu

          if (PERIODIC) call exchange_bc(6,pl_ex)

          do ly=0,pl_ex
              do ib=1,nbp

                  is=dom(ib)%isp
                  ie=dom(ib)%iep
                  js=dom(ib)%jsp
                  je=dom(ib)%jep
                  ks=dom(ib)%ksp
                  ke=dom(ib)%kep

! Boundary Conditions for ksgs
!..............................................................................
!=== West ===>
!..............................................................................
                  if (dom(ib)%iprev<0) then
                  if (dom(ib)%bc_west==4) then
                  do k=ks-1,ke+1
                  do j=js-1,je+1
                          dom(ib)%ksgs(is-1-ly,j,k)= 0.d0                   !k=0 at the boundary?
                      end do
                      end do
                  else if (dom(ib)%bc_west>=63) then                  !Wall functions
                  if (ly==0) call wall_function(1,ib)
                  do k=ks-1,ke+1
                  do j=js-1,je+1
                          dom(ib)%ksgs(is,j,k)=dom(ib)%tauww(j,k)/sqrt(cmu)     !k=ustar**2/cmu**0.5_dp
                          dom(ib)%ksgs(is-1-ly,j,k)=0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_west==61 &
            .or.dom(ib)%bc_west==62) then
                  if (ly==0) call log_law(1,ib)
                  do k=ks-1,ke+1
                  do j=js-1,je+1
                          dom(ib)%ksgs(is,j,k)=dom(ib)%tauww(j,k)/sqrt(cmu)
                          dom(ib)%ksgs(is-1-ly,j,k)=0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_west==1) then
                  do k=ks-1,ke+1
                  do j=js-1,je+1
                          dom(ib)%ksgs(is-1-ly,j,k)= &
                    (3.d0/2.d0)*(ubulk*0.1_dp)**2.0_dp
                      end do
                      end do
                  else if (dom(ib)%bc_west/=5) then                   !if 5->exchange, if 2,3 -> dk/dn=0
                  do k=ks-1,ke+1
                  do j=js-1,je+1
                          dom(ib)%ksgs(is-1-ly,j,k)= dom(ib)%ksgs(is,j,k)
                      end do
                      end do
                  end if
                  end if
!...............................................................................
!=== East ===>
!...............................................................................
                  if (dom(ib)%inext<0) then
                  if (dom(ib)%bc_east==4) then
                  do k=ks-1,ke+1
                  do j=js-1,je+1
                          dom(ib)%ksgs(ie+1+ly,j,k)= 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_east>=63) then                  !Wall functions
                  if (ly==0) call wall_function(2,ib)
                  do k=ks-1,ke+1
                  do j=js-1,je+1
                          dom(ib)%ksgs(ie,j,k)=dom(ib)%tauwe(j,k)/sqrt(cmu)
                          dom(ib)%ksgs(ie+1+ly,j,k)=0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_east==61 &
            .or.dom(ib)%bc_east==62) then
                  if (ly==0) call log_law(2,ib)
                  do k=ks-1,ke+1
                  do j=js-1,je+1
                          dom(ib)%ksgs(ie,j,k)=dom(ib)%tauwe(j,k)/sqrt(cmu)
                          dom(ib)%ksgs(ie+1+ly,j,k)=0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_east/=5) then
                  do k=ks-1,ke+1
                  do j=js-1,je+1
                          dom(ib)%ksgs(ie+1+ly,j,k)= dom(ib)%ksgs(ie,j,k)
                      end do
                      end do
                  end if
                  end if
!...............................................................................
!=== South ===>
!...............................................................................
                  if (dom(ib)%jprev<0) then
                  if (dom(ib)%bc_south==4) then
                  do k=ks-1,ke+1
                  do i=is-1,ie+1
                          dom(ib)%ksgs(i,js-1-ly,k)= 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_south>=63) then                 !Wall functions
                  if (ly==0) call wall_function(3,ib)
                  do k=ks-1,ke+1
                  do i=is-1,ie+1
                          dom(ib)%ksgs(i,js,k)=dom(ib)%tauws(i,k)/sqrt(cmu)
                          dom(ib)%ksgs(i,js-1-ly,k)= 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_south==61 &
            .or.dom(ib)%bc_south==62) then
                  if (ly==0) call log_law(3,ib)
                  do k=ks-1,ke+1
                  do i=is-1,ie+1
                          dom(ib)%ksgs(i,js,k)=dom(ib)%tauws(i,k)/sqrt(cmu)
                          dom(ib)%ksgs(i,js-1-ly,k)= 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_south/=5) then
                  do k=ks-1,ke+1
                  do i=is-1,ie+1
                          dom(ib)%ksgs(i,js-1-ly,k)= dom(ib)%ksgs(i,js,k)
                      end do
                      end do
                  end if
                  end if
!.............................................................................
!=== North ===>
!.............................................................................
                  if (dom(ib)%jnext<0) then
                  if (dom(ib)%bc_north==4) then
                  do k=ks-1,ke+1
                  do i=is-1,ie+1
                          dom(ib)%ksgs(i,je+1+ly,k) = 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_north>=63) then                 !Wall functions
                  if (ly==0) call wall_function(4,ib)
                  do k=ks-1,ke+1
                  do i=is-1,ie+1
                          dom(ib)%ksgs(i,je,k)=dom(ib)%tauwn(i,k)/sqrt(cmu)
                          dom(ib)%ksgs(i,je+1+ly,k) = 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_north==61 &
            .or.dom(ib)%bc_north==62) then
                  if (ly==0) call log_law(4,ib)
                  do k=ks-1,ke+1
                  do i=is-1,ie+1
                          dom(ib)%ksgs(i,je,k)=dom(ib)%tauwn(i,k)/sqrt(cmu)
                          dom(ib)%ksgs(i,je+1+ly,k) = 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_north/=5) then
                  do k=ks-1,ke+1
                  do i=is-1,ie+1
                          dom(ib)%ksgs(i,je+1+ly,k) = dom(ib)%ksgs(i,je,k)
                      end do
                      end do
                  end if
                  end if
!...............................................................................
!=== Bottom ===>
!...............................................................................
                  if (dom(ib)%kprev<0) then
                  if (dom(ib)%bc_bottom==4) then
                  do j=js-1,je+1
                  do i=is-1,ie+1
                          dom(ib)%ksgs(i,j,ks-1-ly)= 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_bottom>=63) then                    !Wall functions
                  if (ly==0) call wall_function(5,ib)
                  do j=js-1,je+1
                  do i=is-1,ie+1
                          dom(ib)%ksgs(i,j,ks)=dom(ib)%tauwb(i,j)/sqrt(cmu)
                          dom(ib)%ksgs(i,j,ks-1-ly)= 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_bottom==61 &
            .or.dom(ib)%bc_bottom==62) then
                  if (ly==0) then
                  call log_law(5,ib)
                  end if
                  do j=js-1,je+1
                  do i=is-1,ie+1
                          dom(ib)%ksgs(i,j,ks)=dom(ib)%tauwb(i,j)/sqrt(cmu)
                          dom(ib)%ksgs(i,j,ks-1-ly)= 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_bottom/=5) then
                  do j=js-1,je+1
                  do i=is-1,ie+1
                          dom(ib)%ksgs(i,j,ks-1-ly)= dom(ib)%ksgs(i,j,ks)
                      end do
                      end do
                  end if
                  end if
!.............................................................................
!=== Top ===>
!.............................................................................
                  if (dom(ib)%knext<0) then
                  if (dom(ib)%bc_top==4) then
                  do j=js-1,je+1
                  do i=is-1,ie+1
                          dom(ib)%ksgs(i,j,ke+1+ly) = 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_top>=63) then                   !Wall functions
                  if (ly==0) call wall_function(6,ib)
                  do j=js-1,je+1
                  do i=is-1,ie+1
                          dom(ib)%ksgs(i,j,ke)=dom(ib)%tauwt(i,j)/sqrt(cmu)
                          dom(ib)%ksgs(i,j,ke+1+ly) =0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_top==61 &
            .or.dom(ib)%bc_top==62) then
                  if (ly==0) call log_law(6,ib)
                  do j=js-1,je+1
                  do i=is-1,ie+1
                          dom(ib)%ksgs(i,j,ke)=dom(ib)%tauwt(i,j)/sqrt(cmu)
                          dom(ib)%ksgs(i,j,ke+1+ly) =0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_top/=5) then
                  do j=js-1,je+1
                  do i=is-1,ie+1
                          dom(ib)%ksgs(i,j,ke+1+ly) = dom(ib)%ksgs(i,j,ke)
                      end do
                      end do
                  end if
                  end if

!==============================================================================
              end do
          end do
      end subroutine boundksgs
!#############################################################################
      subroutine boundeps
!##############################################################################
          use vars
          use multidata
          use module_LSM
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          integer :: i,j,k,ib,ly
          integer :: is,ie,js,je,ks,ke
          real(dp) :: rk1,drkdy,lz
          real(dp) :: dx,dy,dz,kappa,delta

          kappa=0.41_dp
          lz=zen-zst
          if (L_LSM) lz=length
          if (L_LSMbase) lz=length

          if (PERIODIC) call exchange_bc(9,pl_ex)

          do ly=0,pl_ex
              do ib=1,nbp

                  dx=dom(ib)%dx
                  dy=dom(ib)%dy
                  dz=dom(ib)%dz

                  is=dom(ib)%isp
                  ie=dom(ib)%iep
                  js=dom(ib)%jsp
                  je=dom(ib)%jep
                  ks=dom(ib)%ksp
                  ke=dom(ib)%kep

! Boundary Conditions for epsylon
!..............................................................................
!=== West ===>
!..............................................................................
                  if (dom(ib)%iprev<0) then
                  if (dom(ib)%bc_west==4) then
                  do k=ks-1,ke+1
                  do j=js-1,je+1
                          rk1 = sqrt(dom(ib)%ksgs(is,j,k))
                          drkdy = rk1 /(0.5_dp*dom(ib)%dx)
                          dom(ib)%eps(is-1-ly,j,k)= 2.d0*rrey*drkdy**2.0_dp
                      end do
                      end do
                  else if (dom(ib)%bc_west>=63) then                  !Wall functions
                  if (ly==0) call wall_function(1,ib)
                  do k=ks-1,ke+1
                  do j=js-1,je+1
                          dom(ib)%eps(is,j,k)= &
                    dom(ib)%tauww(j,k)**1.5_dp/(kappa*dx)
                          dom(ib)%eps(is-1-ly,j,k)= 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_west==61 &
            .or.dom(ib)%bc_west==62) then
                  if (ly==0) call log_law(1,ib)
                  do k=ks-1,ke+1
                  do j=js-1,je+1
                          dom(ib)%eps(is,j,k)= &
                    dom(ib)%tauww(j,k)**1.5_dp/(kappa*dx)
                          dom(ib)%eps(is-1-ly,j,k)= 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_west==1) then
                  do k=ks-1,ke+1
                  do j=js-1,je+1
                          dom(ib)%eps(is-1-ly,j,k) = &
                    0.09_dp**0.75_dp*dom(ib)%ksgs(is-1-ly,j,k)**1.5_dp/(0.07_dp*lz)
                      end do
                      end do
                  else if (dom(ib)%bc_west/=5) then
                  do k=ks-1,ke+1
                  do j=js-1,je+1
                          dom(ib)%eps(is-1-ly,j,k)= dom(ib)%eps(is,j,k)
                      end do
                      end do
                  end if
                  end if
!...............................................................................
!=== East ===>
!...............................................................................
                  if (dom(ib)%inext<0) then
                  if (dom(ib)%bc_east==4) then
                  do k=ks-1,ke+1
                  do j=js-1,je+1
                          rk1 = sqrt(dom(ib)%ksgs(ie,j,k))
                          drkdy = rk1 /(0.5_dp*dom(ib)%dx)
                          dom(ib)%eps(ie+1+ly,j,k)= 2.d0*rrey*drkdy**2.0_dp
                      end do
                      end do
                  else if (dom(ib)%bc_east>=63) then                  !Wall functions
                  if (ly==0) call wall_function(2,ib)
                  do k=ks-1,ke+1
                  do j=js-1,je+1
                          dom(ib)%eps(ie,j,k)= &
                    dom(ib)%tauwe(j,k)**1.5_dp/(kappa*dx)
                          dom(ib)%eps(ie+1+ly,j,k)= 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_east==61 &
            .or.dom(ib)%bc_east==62) then
                  if (ly==0) call log_law(2,ib)
                  do k=ks-1,ke+1
                  do j=js-1,je+1
                          dom(ib)%eps(ie,j,k)= &
                    dom(ib)%tauwe(j,k)**1.5_dp/(kappa*dx)
                          dom(ib)%eps(ie+1+ly,j,k)= 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_east/=5) then
                  do k=ks-1,ke+1
                  do j=js-1,je+1
                          dom(ib)%eps(ie+1+ly,j,k)= dom(ib)%eps(ie,j,k)
                      end do
                      end do
                  end if
                  end if
!...............................................................................
!=== South ===>
!...............................................................................
                  if (dom(ib)%jprev<0) then
                  if (dom(ib)%bc_south==4) then
                  do k=ks-1,ke+1
                  do i=is-1,ie+1
                          rk1 = sqrt(dom(ib)%ksgs(i,js,k))
                          drkdy = rk1 /(0.5_dp*dom(ib)%dy)
                          dom(ib)%eps(i,js-1-ly,k)= 2.d0*rrey*drkdy**2.0_dp
                      end do
                      end do
                  else if (dom(ib)%bc_south>=63) then                 !Wall functions
                  if (ly==0) call wall_function(3,ib)
                  do k=ks-1,ke+1
                  do i=is-1,ie+1
                          dom(ib)%eps(i,js,k)= &
                    dom(ib)%tauws(i,k)**1.5_dp/(kappa*dy)
                          dom(ib)%eps(i,js-1-ly,k)= 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_south==61 &
            .or.dom(ib)%bc_south==62) then
                  if (ly==0) call log_law(3,ib)
                  do k=ks-1,ke+1
                  do i=is-1,ie+1
                          dom(ib)%eps(i,js,k)= &
                    dom(ib)%tauws(i,k)**1.5_dp/(kappa*dy)
                          dom(ib)%eps(i,js-1-ly,k)= 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_south/=5) then
                  do k=ks-1,ke+1
                  do i=is-1,ie+1
                          dom(ib)%eps(i,js-1-ly,k)= dom(ib)%eps(i,js,k)
                      end do
                      end do
                  end if
                  end if
!.............................................................................
!=== North ===>
!.............................................................................
                  if (dom(ib)%jnext<0) then
                  if (dom(ib)%bc_north==4) then
                  do k=ks-1,ke+1
                  do i=is-1,ie+1
                          rk1 = sqrt(dom(ib)%ksgs(i,je,k))
                          drkdy = rk1 /(0.5_dp*dom(ib)%dy)
                          dom(ib)%eps(i,je+1+ly,k)= 2.d0*rrey*drkdy**2.0_dp
                      end do
                      end do
                  else if (dom(ib)%bc_north>=63) then                 !Wall functions
                  if (ly==0) call wall_function(4,ib)
                  do k=ks-1,ke+1
                  do i=is-1,ie+1
                          dom(ib)%eps(i,je,k)= &
                    dom(ib)%tauwn(i,k)**1.5_dp/(kappa*dy)
                          dom(ib)%eps(i,je+1+ly,k) = 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_north==61 &
            .or.dom(ib)%bc_north==62) then
                  if (ly==0) call log_law(4,ib)
                  do k=ks-1,ke+1
                  do i=is-1,ie+1
                          dom(ib)%eps(i,je,k)= &
                    dom(ib)%tauwn(i,k)**1.5_dp/(kappa*dy)
                          dom(ib)%eps(i,je+1+ly,k) = 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_north/=5) then
                  do k=ks-1,ke+1
                  do i=is-1,ie+1
                          dom(ib)%eps(i,je+1+ly,k) = dom(ib)%eps(i,je,k)
                      end do
                      end do
                  end if
                  end if
!...............................................................................
!=== Bottom ===>
!...............................................................................
                  if (dom(ib)%kprev<0) then
                  if (dom(ib)%bc_bottom==4) then
                  do j=js-1,je+1
                  do i=is-1,ie+1
                          rk1 = sqrt(dom(ib)%ksgs(i,j,ks))
                          drkdy = rk1 /(0.5_dp*dom(ib)%dz)
                          dom(ib)%eps(i,j,ks-1-ly)= 2.d0*rrey*drkdy**2.0_dp
                      end do
                      end do
                  else if (dom(ib)%bc_bottom>=63) then                    !Wall functions
                  if (ly==0) call wall_function(5,ib)
                  do j=js-1,je+1
                  do i=is-1,ie+1
                          dom(ib)%eps(i,j,ks)= &
                    dom(ib)%tauwb(i,j)**1.5_dp/(kappa*dz)
                          dom(ib)%eps(i,j,ks-1-ly)= 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_bottom==61 &
            .or.dom(ib)%bc_bottom==62) then
                  if (ly==0) then
                  call log_law(5,ib)
                  end if
                  do j=js-1,je+1
                  do i=is-1,ie+1
                          dom(ib)%eps(i,j,ks)= &
                    dom(ib)%tauwb(i,j)**1.5_dp/(kappa*dz)
                          dom(ib)%eps(i,j,ks-1-ly)= 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_bottom/=5) then
                  do j=js-1,je+1
                  do i=is-1,ie+1
                          dom(ib)%eps(i,j,ks-1-ly)= dom(ib)%eps(i,j,ks)
                      end do
                      end do
                  end if
                  end if
!.............................................................................
!=== Top ===>
!.............................................................................
                  if (dom(ib)%knext<0) then
                  if (dom(ib)%bc_top==4) then
                  do j=js-1,je+1
                  do i=is-1,ie+1
                          rk1 = sqrt(dom(ib)%ksgs(i,j,ke))
                          drkdy = rk1 /(0.5_dp*dom(ib)%dz)
                          dom(ib)%eps(i,j,ke+1+ly)= 2.d0*rrey*drkdy**2.0_dp
                      end do
                      end do
                  else if (dom(ib)%bc_top>=63) then                   !Wall functions
                  if (ly==0) call wall_function(6,ib)
                  do j=js-1,je+1
                  do i=is-1,ie+1
                          dom(ib)%eps(i,j,ke)= &
                    dom(ib)%tauwt(i,j)**1.5_dp/(kappa*dz)
                          dom(ib)%eps(i,j,ke+1+ly) = 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_top==61 &
            .or.dom(ib)%bc_top==62) then
                  if (ly==0) call log_law(6,ib)
                  do j=js-1,je+1
                  do i=is-1,ie+1
                          dom(ib)%eps(i,j,ke)= &
                    dom(ib)%tauwt(i,j)**1.5_dp/(kappa*dz)
                          dom(ib)%eps(i,j,ke+1+ly) = 0.d0
                      end do
                      end do
                  else if (dom(ib)%bc_top/=5) then
                  do j=js-1,je+1
                  do i=is-1,ie+1
                          dom(ib)%eps(i,j,ke+1+ly) = dom(ib)%eps(i,j,ke)
                      end do
                      end do
                  end if
                  end if

!==============================================================================
              end do
          end do
      end subroutine boundeps
!#############################################################################
