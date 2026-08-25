!======================================================================!
!                 LAGRANGIAN PARTICLE TRACKING                         !
!----------------------------------------------------------------------!
!                           Bruño Fraga                                !
!                      Cardiff Uni 2013-2017                           !
!                        Stanford Uni 2018                             !
!                    Uni of Birmingham 2019-2024                       !
!======================================================================!
!######################################################################!
module multiflow3d_LPT
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use mpi_f08, only : mpi_comm_world, mpi_double_precision
  use multidata, only: multidom
  use multiflow3d_collison, only: collision_particle, collision_walls
  use multiflow3d_mpi, only : ierr, nprocs
  use omp_lib, only : omp_get_num_threads, &
                      omp_get_thread_num, &
                      omp_set_num_threads
  implicit none
  private

  public :: particle_tracking, final_LPT

contains

  subroutine particle_tracking(wop_pt, vop_pt, id, Lcol,Lcolwall, PSIcell, &
       rhop_loc, np_loc, xp_loc,yp_loc,zp_loc, uop_loc,vop_loc, &
       wop_loc,dp_loc, uop_pt, up_pt,vp_pt,wp_pt, Fpu,Fpv,Fpw, order, &
       pl, Re,las,lenergy,alfapr,dens,dt,gx,gy,gz,nbp,dom_id,dom,&
       xpg_loc,ypg_loc,zpg_loc,vopg_loc, wopg_loc,rhopg_loc,npg_loc,&
       uopg_loc,dpg_loc,k_n,yst,yen,zst,zen,bc_s,bc_n,bc_b,bc_t)
  !
  !     Calculates particles' velocities and the resulting source terms
  !
  !######################################################################!

    real(dp), dimension(:), intent(in) :: uop_pt, wop_pt, vop_pt
    integer, allocatable, dimension(:), intent(inout) ::  id
    real(dp), allocatable, dimension(:), intent(in) :: rhop_loc
    logical , intent(in) :: PSIcell, Lcol,Lcolwall
    integer , intent(in) :: np_loc
    real(dp), allocatable, dimension(:), intent(inout):: xp_loc,yp_loc,zp_loc
    real(dp),  allocatable,  dimension(:), intent(inout) :: uop_loc,vop_loc, wop_loc
    real(dp),  allocatable, dimension(:), intent(in) :: dp_loc
    real(dp), allocatable, dimension(:), intent(inout) :: up_pt,vp_pt,wp_pt
    real(dp), allocatable, dimension(:), intent(inout) :: Fpu,Fpv,Fpw
    integer , intent(in) :: order, pl
    real(dp) , intent(in) :: Re
    logical , intent(in) :: las, lenergy
    real(dp) , intent(in) :: alfapr, dens
    real(dp) , intent(in) :: dt,gx,gy,gz
    integer , intent(in) :: nbp
    integer, allocatable, dimension(:), intent(in) :: dom_id
    type (multidom), pointer, dimension(:), intent(in) :: dom

    real(dp), allocatable, dimension(:),intent(inout):: xpg_loc,ypg_loc,zpg_loc, &
                                          vopg_loc, wopg_loc,rhopg_loc
    integer,intent(in) :: npg_loc
    real(dp), allocatable, dimension(:),intent(inout):: uopg_loc,dpg_loc
    real(dp) ,intent(in):: k_n
    real(dp) ,intent(in):: yst,yen,zst,zen
    integer ,intent(in):: bc_s,bc_n,bc_b,bc_t

    integer :: i,j,k,l
    integer :: ib,is,ie,js,je,ks,ke
    integer :: nt,m
    integer :: iballs_u,iballe_u,jballs_u,jballe_u,kballs_u,kballe_u
    integer :: iballs_v,iballe_v,jballs_v,jballe_v,kballs_v,kballe_v
    integer :: iballs_w,iballe_w,jballs_w,jballe_w,kballs_w,kballe_w
    real(dp) :: REp,rx,ry,rz
    real(dp) :: a,b,c,wx,wy,wz,Cd,ao,bo,co
    real(dp) :: dwdy,dvdz,dudz,dwdx,dvdx,dudy
    real(dp) :: dh,delta,gamma_p,ddh,ddelta
    real(dp) :: Vcell,Vp  !,Vball
    real(dp), allocatable, dimension(:):: ui_pt,vi_pt,wi_pt
    real(dp), allocatable, dimension(:):: uoi_pt,voi_pt,woi_pt
    integer,allocatable,dimension(:)::  ip,jp,kp,ipu,jpv,kpw

    allocate (ui_pt(np_loc),vi_pt(np_loc),wi_pt(np_loc))
    allocate (uoi_pt(np_loc),voi_pt(np_loc),woi_pt(np_loc))
    allocate (ip(np_loc),jp(np_loc),kp(np_loc))
    allocate (ipu(np_loc),jpv(np_loc),kpw(np_loc))
    allocate (up_pt(np_loc),vp_pt(np_loc),wp_pt(np_loc))
    allocate (Fpu(np_loc),Fpv(np_loc),Fpw(np_loc))

    if (np_loc<=100) nt = 1

    call OMP_SET_NUM_THREADS(nt)

    select case (order)
      case (1)
         m = 1  !1.5d0
      case (2)
         m = 2  !2.5d0
      case (3)
         m = 2  !2.d0
      case (4)
         m = 2  !2.5d0
      case (5)
         m = 1  !1.5d0
      case (6)
         m = 2  !2.d0
    end select

    !loop in domains
    do ib=1,nbp

       Vcell = dom(ib)%dx*dom(ib)%dy*dom(ib)%dz

       !computational domain limits (one-axis index)
       is = dom(ib)%isp
       ie = dom(ib)%iep
       js = dom(ib)%jsp
       je = dom(ib)%jep
       ks = dom(ib)%ksp
       ke = dom(ib)%kep

       !loop in particles
       !$OMP       PARALLEL DEFAULT (SHARED), PRIVATE(i,j,k,l,&
       !$OMP      iballs_u,iballe_u,jballs_u,jballe_u,kballs_u,kballe_u,&
       !$OMP      iballs_v,iballe_v,jballs_v,jballe_v,kballs_v,kballe_v,&
       !$OMP      iballs_w,iballe_w,jballs_w,jballe_w,kballs_w,kballe_w,&
       !$OMP      REp,rx,ry,rz,Vp,delta,gamma_p,&
       !$OMP      a,b,c,ao,bo,co,Cd,wx,wy,wz,&
       !$OMP      dwdy,dvdz,dudz,dvdx,dudy,dwdx)

       !$OMP DO SCHEDULE (DYNAMIC,1)
       do l=1,np_loc

          if (id(l)==dom_id(ib)) then       !particle belongs to THIS block

             Vp = 3.1416_dp*dp_loc(l)**3.0d0/6.0d0

             ip(l)=INT((xp_loc(l)-dom(ib)%x(is-1)-1.0d-12)/dom(ib)%dx)+1+pl
             jp(l)=INT((yp_loc(l)-dom(ib)%y(js-1)-1.0d-12)/dom(ib)%dy)+1+pl
             kp(l)=INT((zp_loc(l)-dom(ib)%z(ks-1)-1.0d-12)/dom(ib)%dz)+1+pl


             !locate the u,v and w nodes

             if (xp_loc(l)>dom(ib)%xc(ip(l))) then

                ipu(l) = ip(l)

             else if (xp_loc(l)<=dom(ib)%xc(ip(l))) then

                ipu(l) = ip(l) - 1

             end if

             if (yp_loc(l)>dom(ib)%yc(jp(l))) then

                jpv(l) = jp(l)

             else if (yp_loc(l)<=dom(ib)%yc(jp(l))) then

                jpv(l) = jp(l) - 1

             end if

             if (zp_loc(l)>dom(ib)%zc(kp(l))) then

                kpw(l) = kp(l)

             else if (zp_loc(l)<=dom(ib)%zc(kp(l))) then

                kpw(l) = kp(l) - 1

             end if


             rx = max(dp_loc(l),dom(ib)%dx)
             ry = max(dp_loc(l),dom(ib)%dy)
             rz = max(dp_loc(l),dom(ib)%dz)

             !Ball

             if (order==3.or.order==6) then
                if (ipu(l)==ip(l)) then
                   iballs_u = ipu(l) - 1 * NINT(rx/dom(ib)%dx)
                   iballe_u = ipu(l) + m * NINT(rx/dom(ib)%dx)
                   iballs_v = ip(l) - 1 * NINT(rx/dom(ib)%dx)
                   iballe_v = ip(l) + m * NINT(rx/dom(ib)%dx)
                   iballs_w = ip(l) - 1 * NINT(rx/dom(ib)%dx)
                   iballe_w = ip(l) + m * NINT(rx/dom(ib)%dx)
                else
                   iballs_u = ipu(l) - m * NINT(rx/dom(ib)%dx)
                   iballe_u = ipu(l) + 1 * NINT(rx/dom(ib)%dx)
                   iballs_v = ip(l) - m * NINT(rx/dom(ib)%dx)
                   iballe_v = ip(l) + 1 * NINT(rx/dom(ib)%dx)
                   iballs_w = ip(l) - m * NINT(rx/dom(ib)%dx)
                   iballe_w = ip(l) + 1 * NINT(rx/dom(ib)%dx)
                end if
                if (jpv(l)==jp(l)) then
                   jballs_u = jp(l) - 1 * NINT(ry/dom(ib)%dy)
                   jballe_u = jp(l) + m * NINT(ry/dom(ib)%dy)
                   jballs_v = jpv(l) - 1 * NINT(ry/dom(ib)%dy)
                   jballe_v = jpv(l) + m * NINT(ry/dom(ib)%dy)
                   jballs_w = jp(l) - 1 * NINT(ry/dom(ib)%dy)
                   jballe_w = jp(l) + m * NINT(ry/dom(ib)%dy)
                else
                   jballs_u = jp(l) - m * NINT(ry/dom(ib)%dy)
                   jballe_u = jp(l) + 1 * NINT(ry/dom(ib)%dy)
                   jballs_v = jpv(l) - m * NINT(ry/dom(ib)%dy)
                   jballe_v = jpv(l) + 1 * NINT(ry/dom(ib)%dy)
                   jballs_w = jp(l) - m * NINT(ry/dom(ib)%dy)
                   jballe_w = jp(l) + 1 * NINT(ry/dom(ib)%dy)
                end if
                if (kpw(l)==kp(l)) then
                   kballs_u = kp(l) - 1 * NINT(rz/dom(ib)%dz)
                   kballe_u = kp(l) + m * NINT(rz/dom(ib)%dz)
                   kballs_v = kp(l) - 1 * NINT(rz/dom(ib)%dz)
                   kballe_v = kp(l) + m * NINT(rz/dom(ib)%dz)
                   kballs_w = kpw(l) - 1 * NINT(rz/dom(ib)%dz)
                   kballe_w = kpw(l) + m * NINT(rz/dom(ib)%dz)
                else
                   kballs_u = kp(l) - m * NINT(rz/dom(ib)%dz)
                   kballe_u = kp(l) + 1 * NINT(rz/dom(ib)%dz)
                   kballs_v = kp(l) - m * NINT(rz/dom(ib)%dz)
                   kballe_v = kp(l) + 1 * NINT(rz/dom(ib)%dz)
                   kballs_w = kpw(l) - m * NINT(rz/dom(ib)%dz)
                   kballe_w = kpw(l) + 1 * NINT(rz/dom(ib)%dz)
                end if
             else
                iballs_u = ipu(l) - m * NINT(rx/dom(ib)%dx)
                iballe_u = ipu(l) + m * NINT(rx/dom(ib)%dx)
                jballs_u = jp(l) - m * NINT(ry/dom(ib)%dy)
                jballe_u = jp(l) + m * NINT(ry/dom(ib)%dy)
                kballs_u = kp(l) - m * NINT(rz/dom(ib)%dz)
                kballe_u = kp(l) + m * NINT(rz/dom(ib)%dz)

                iballs_v = ip(l) - m * NINT(rx/dom(ib)%dx)
                iballe_v = ip(l) + m * NINT(rx/dom(ib)%dx)
                jballs_v = jpv(l) - m * NINT(ry/dom(ib)%dy)
                jballe_v = jpv(l) + m * NINT(ry/dom(ib)%dy)
                kballs_v = kp(l) - m * NINT(rz/dom(ib)%dz)
                kballe_v = kp(l) + m * NINT(rz/dom(ib)%dz)

                iballs_w = ip(l) - m * NINT(rx/dom(ib)%dx)
                iballe_w = ip(l) + m * NINT(rx/dom(ib)%dx)
                jballs_w = jp(l) - m * NINT(ry/dom(ib)%dy)
                jballe_w = jp(l) + m * NINT(ry/dom(ib)%dy)
                kballs_w = kpw(l) - m * NINT(rz/dom(ib)%dz)
                kballe_w = kpw(l) + m * NINT(rz/dom(ib)%dz)
             end if

             iballs_u = max(iballs_u,1)
             iballe_u = min(iballe_u,dom(ib)%ttc_i)
             jballs_u = max(jballs_u,1)
             jballe_u = min(jballe_u,dom(ib)%ttc_j)
             kballs_u = max(kballs_u,1)
             kballe_u = min(kballe_u,dom(ib)%ttc_k)
             iballs_v = max(iballs_v,1)
             iballe_v = min(iballe_v,dom(ib)%ttc_i)
             jballs_v = max(jballs_v,1)
             jballe_v = min(jballe_v,dom(ib)%ttc_j)
             kballs_v = max(kballs_v,1)
             kballe_v = min(kballe_v,dom(ib)%ttc_k)
             iballs_w = max(iballs_w,1)
             iballe_w = min(iballe_w,dom(ib)%ttc_i)
             jballs_w = max(jballs_w,1)
             jballe_w = min(jballe_w,dom(ib)%ttc_j)
             kballs_w = max(kballs_w,1)
             kballe_w = min(kballe_w,dom(ib)%ttc_k)

             uoi_pt(l) = 0.0d0
             voi_pt(l) = 0.0d0
             woi_pt(l) = 0.0d0

             do i=iballs_u,iballe_u
                do j=jballs_u,jballe_u
                   do k=kballs_u,kballe_u

                      uoi_pt(l) = uoi_pt(l) + dom(ib)%uoo(i,j,k)* &
                           dh(rx,ry,rz,dom(ib)%x(i),dom(ib)%yc(j) &
                           ,dom(ib)%zc(k),xp_loc(l),yp_loc(l),zp_loc(l),order)

                   end do
                end do
             end do

             do i=iballs_v,iballe_v
                do j=jballs_v,jballe_v
                   do k=kballs_v,kballe_v

                      voi_pt(l) = voi_pt(l) + dom(ib)%voo(i,j,k)* &
                           dh(rx,ry,rz,dom(ib)%xc(i),dom(ib)%y(j) &
                           ,dom(ib)%zc(k),xp_loc(l),yp_loc(l),zp_loc(l),order)

                   end do
                end do
             end do

             do i=iballs_w,iballe_w
                do j=jballs_w,jballe_w
                   do k=kballs_w,kballe_w

                      woi_pt(l) = woi_pt(l) + dom(ib)%woo(i,j,k)* &
                           dh(rx,ry,rz,dom(ib)%xc(i),dom(ib)%yc(j) &
                           ,dom(ib)%z(k),xp_loc(l),yp_loc(l),zp_loc(l),order)

                   end do
                end do
             end do

             up_pt(l) = uop_loc(l)
             vp_pt(l) = vop_loc(l)
             wp_pt(l) = wop_loc(l)

             ui_pt(l) = 0.0d0
             vi_pt(l) = 0.0d0
             wi_pt(l) = 0.0d0

             do i=iballs_u,iballe_u
                do j=jballs_u,jballe_u
                   do k=kballs_u,kballe_u

                      delta = dh(rx,ry,rz,dom(ib)%x(i),dom(ib)%yc(j) &
                           ,dom(ib)%zc(k),xp_loc(l),yp_loc(l),zp_loc(l),order)

                      ui_pt(l) = ui_pt(l) + dom(ib)%ustar(i,j,k) * delta

                   end do
                end do
             end do

             do i=iballs_v,iballe_v
                do j=jballs_v,jballe_v
                   do k=kballs_v,kballe_v

                      delta = dh(rx,ry,rz,dom(ib)%xc(i),dom(ib)%y(j) &
                           ,dom(ib)%zc(k),xp_loc(l),yp_loc(l),zp_loc(l),order)

                      vi_pt(l) = vi_pt(l) + dom(ib)%vstar(i,j,k) * delta

                   end do
                end do
             end do

             do i=iballs_w,iballe_w
                do j=jballs_w,jballe_w
                   do k=kballs_w,kballe_w

                      delta = dh(rx,ry,rz,dom(ib)%xc(i),dom(ib)%yc(j) &
                           ,dom(ib)%z(k),xp_loc(l),yp_loc(l),zp_loc(l),order)

                      wi_pt(l) = wi_pt(l) + dom(ib)%wstar(i,j,k) * delta

                   end do
                end do
             end do

             !Slip vel components
             a = up_pt(l)-ui_pt(l)
             b = vp_pt(l)-vi_pt(l)
             c = wp_pt(l)-wi_pt(l)

             REp = dp_loc(l)* (sqrt((uop_loc(l)-ui_pt(l))**2.0d0+(vop_loc(l) &
                  -vi_pt(l))**2.0d0+(wop_loc(l)-wi_pt(l))**2.0d0))/(1.0d0/Re)

             if (REp<=800) Cd = 24.0d0*(1.0d0+0.15d0*(REp**0.687d0))/REp
             if (REp>800) Cd = 0.44d0

             !Vorticity calculation
             dudy = 0.0_dp
             do i=iballs_u,iballe_u
                do j=jballs_u,jballe_u
                   do k=kballs_u,kballe_u

                      ddelta = ddh(rx,ry,rz,dom(ib)%x(i),dom(ib)%yc(j) &
                           ,dom(ib)%zc(k),xp_loc(l),yp_loc(l),zp_loc(l),order,2)

                      dudy = dudy + dom(ib)%uoo(i,j,k)*ddelta

                   end do
                end do
             end do


             dudz = 0.0_dp
             do i=iballs_u,iballe_u
                do j=jballs_u,jballe_u
                   do k=kballs_u,kballe_u

                      ddelta = ddh(rx,ry,rz,dom(ib)%x(i),dom(ib)%yc(j) &
                           ,dom(ib)%zc(k),xp_loc(l),yp_loc(l),zp_loc(l),order,3)

                      dudz = dudz + dom(ib)%uoo(i,j,k)*ddelta

                   end do
                end do
             end do


             dvdx = 0.0_dp
             do i=iballs_v,iballe_v
                do j=jballs_v,jballe_v
                   do k=kballs_v,kballe_v

                      ddelta = ddh(rx,ry,rz,dom(ib)%xc(i),dom(ib)%y(j) &
                           ,dom(ib)%zc(k),xp_loc(l),yp_loc(l),zp_loc(l),order,1)

                      dvdx = dvdx + dom(ib)%voo(i,j,k)*ddelta

                   end do
                end do
             end do


             dvdz = 0.0_dp
             do i=iballs_v,iballe_v
                do j=jballs_v,jballe_v
                   do k=kballs_v,kballe_v

                      ddelta = ddh(rx,ry,rz,dom(ib)%xc(i),dom(ib)%y(j) &
                           ,dom(ib)%zc(k),xp_loc(l),yp_loc(l),zp_loc(l),order,3)

                      dvdz = dvdz + dom(ib)%voo(i,j,k)*ddelta

                   end do
                end do
             end do


             dwdx = 0.0_dp
             do i=iballs_w,iballe_w
                do j=jballs_w,jballe_w
                   do k=kballs_w,kballe_w

                      ddelta = ddh(rx,ry,rz,dom(ib)%xc(i),dom(ib)%yc(j) &
                           ,dom(ib)%z(k),xp_loc(l),yp_loc(l),zp_loc(l),order,1)

                      dwdx = dwdx + dom(ib)%woo(i,j,k)*ddelta

                   end do
                end do
             end do


             dwdy = 0.0_dp
             do i=iballs_w,iballe_w
                do j=jballs_w,jballe_w
                   do k=kballs_w,kballe_w

                      ddelta = ddh(rx,ry,rz,dom(ib)%xc(i),dom(ib)%yc(j) &
                           ,dom(ib)%z(k),xp_loc(l),yp_loc(l),zp_loc(l),order,2)

                      dwdy = dwdy + dom(ib)%woo(i,j,k)*ddelta

                   end do
                end do
             end do

             wx = dwdy-dvdz
             wy = dudz-dwdx
             wz = dvdx-dudy

             if (LENERGY.or.LAS) then                                        !variable density form
                gamma_p=rhop_loc(l)/dom(ib)%dens(ip(l),jp(l),kp(l))               !variable density
             else
                gamma_p=rhop_loc(l)/dens                                          ! constant density
             end if

             if ((dp_loc(l))<0.00001_dp) then
                !Particles with dp<10um treated as passive Aleks 05/2022
                up_pt(l) = ui_pt(l)
                vp_pt(l) = vi_pt(l)
                wp_pt(l) = wi_pt(l)
             else
                up_pt(l) = uop_loc(l) + dt * &
                     (gx*(gamma_p-1.0d0)/(gamma_p+0.5_dp)+&                       !Buoyancy
                     (((1.0_dp+0.5_dp)/(gamma_p+0.5_dp))*((ui_pt(l)-uoi_pt(l))/dt) &
                     -(3.0d0/(4.0d0*dp_loc(l)*(gamma_p+0.5_dp))) &
                     *Cd*sqrt(a**2.0d0+b**2.0d0+c**2.0d0)*a &
                     -(1.0_dp/(gamma_p+0.5_dp))*0.53d0*(b*wz-c*wy)))


                vp_pt(l) = vop_loc(l) + dt* &
                     (gy*(gamma_p-1.0d0)/(gamma_p+0.5_dp)+&                       !Buoyancy
                     (((1.0_dp+0.5_dp)/(gamma_p+0.5_dp))*((vi_pt(l)-voi_pt(l))/dt) &
                     -(3.0d0/(4.0d0*dp_loc(l)*(gamma_p+0.5_dp))) &
                     *Cd*sqrt(a**2.0d0+b**2.0d0+c**2.0d0)*b &
                     -(1.0_dp/(gamma_p+0.5_dp))*0.53d0*(c*wx-a*wz)))


                wp_pt(l) = wop_loc(l) + dt* &
                     (gz*(gamma_p-1.0d0)/(gamma_p+0.5_dp)+&                      !Buoyancy
                     ((1.0_dp+0.5_dp)/(gamma_p+0.5_dp))*((wi_pt(l)-woi_pt(l))/dt)&!Fluid stress
                     -(3.0d0/(4.0d0*dp_loc(l)*(gamma_p+0.5_dp)))&           !Added Mass and drag
                     *Cd*sqrt(a**2.0d0+b**2.0d0+c**2.0d0)*c&                !Added Mass and drag
                     -(1.0_dp/(gamma_p+0.5_dp))*0.53d0*(a*wy-b*wx))                     !Lift
             end if

             if (Lcolwall) then
                 !updating particle velocities based on collisions with walls
                call collision_walls(l,rhop_loc,yp_loc,zp_loc,dp_loc,&
                             up_pt,vp_pt,wp_pt,k_n,yst,yen,zst,zen,bc_s,bc_n,bc_b,bc_t,dt)
             end if
             if (Lcol) then
                 !updating particle velocities based on p2p collisions
                 call collision_particle(l,xpg_loc,ypg_loc,zpg_loc,vopg_loc, wopg_loc,rhopg_loc,&
                                xp_loc,yp_loc,zp_loc,uop_loc,vop_loc, wop_loc,id,rhop_loc, &
                                dp_loc,up_pt,vp_pt,wp_pt,np_loc,npg_loc,uopg_loc,dpg_loc,k_n,dt,nbp,&
                                dom_id)
             end if


             if ((dp_loc(l))>=0.00001_dp) then  !only do calcs if dp>=10um
                !Update slip velocity
                a = up_pt(l)-ui_pt(l)
                b = vp_pt(l)-vi_pt(l)
                c = wp_pt(l)-wi_pt(l)

                  Fpu(l) = -(((1.0_dp+0.5_dp)/(gamma_p+0.5_dp))*((ui_pt(l)-uoi_pt(l))/dt) &
                  -(3.0d0/(4.0d0*dp_loc(l)*(gamma_p+0.5_dp))) &
                  *Cd*sqrt(a**2.0d0+b**2.0d0+c**2.0d0)*a &
                  -(1.0_dp/(gamma_p+0.5_dp))*0.53d0*(b*wz-c*wy))

                  Fpv(l) = -(((1.0_dp+0.5_dp)/(gamma_p+0.5_dp))*((vi_pt(l)-voi_pt(l))/dt) &
                  -(3.0d0/(4.0d0*dp_loc(l)*(gamma_p+0.5_dp))) &
                  *Cd*sqrt(a**2.0d0+b**2.0d0+c**2.0d0)*b &
                  -(1.0_dp/(gamma_p+0.5_dp))*0.53d0*(c*wx-a*wz))

                  !Fpw(l) =-(((1.0_dp-gamma_p)/(gamma_p+0.5_dp))*9.81d0+                   !Buoyancy
               Fpw(l) =-(((1.0_dp+0.5_dp)/(gamma_p+0.5_dp))*((wi_pt(l)-woi_pt(l))/dt)& !Fluid stress
                  -(3.0d0/(4.0d0*dp_loc(l)*(gamma_p+0.5_dp)))&  !Added Mass and drag
                  *Cd*sqrt(a**2.0d0+b**2.0d0+c**2.0d0)*c&       !Added Mass and drag
                  -(1.0_dp/(gamma_p+0.5_dp))*0.53d0*(a*wy-b*wx))                     !Lift

                  !$OMP CRITICAL
                  if (PSIcell) then

                     dom(ib)%ustar(ipu(l),jp(l),kp(l)) = &
                          dom(ib)%ustar(ipu(l),jp(l),kp(l)) + dt * alfapr * Fpu(l) * &
                          Vp/Vcell

                     dom(ib)%vstar(ip(l),jpv(l),kp(l)) = &
                          dom(ib)%vstar(ip(l),jpv(l),kp(l)) + dt * alfapr * Fpv(l) * &
                          Vp/Vcell

                     dom(ib)%wstar(ip(l),jp(l),kpw(l)) = &
                          dom(ib)%wstar(ip(l),jp(l),kpw(l)) + dt * alfapr * Fpw(l) * &
                          Vp/Vcell
                  else
                     do i=iballs_u,iballe_u
                        do j=jballs_u,jballe_u
                           do k=kballs_u,kballe_u

                              delta = dh(rx,ry,rz,dom(ib)%x(i),dom(ib)%yc(j) &
                                   ,dom(ib)%zc(k),xp_loc(l),yp_loc(l),zp_loc(l),order)

                              dom(ib)%ustar(i,j,k) = &
                                   dom(ib)%ustar(i,j,k) + dt * alfapr * Fpu(l) * delta * &
                                   Vp/Vcell

                           end do
                        end do
                     end do


                     do i=iballs_v,iballe_v
                        do j=jballs_v,jballe_v
                           do k=kballs_v,kballe_v

                              delta = dh(rx,ry,rz,dom(ib)%xc(i),dom(ib)%y(j) &
                                   ,dom(ib)%zc(k),xp_loc(l),yp_loc(l),zp_loc(l),order)

                              dom(ib)%vstar(i,j,k) = &
                                   dom(ib)%vstar(i,j,k) + dt * alfapr * Fpv(l) * delta * &
                                   Vp/Vcell

                           end do
                        end do
                     end do


                     do i=iballs_w,iballe_w
                        do j=jballs_w,jballe_w
                           do k=kballs_w,kballe_w

                              delta = dh(rx,ry,rz,dom(ib)%xc(i),dom(ib)%yc(j) &
                                   ,dom(ib)%z(k),xp_loc(l),yp_loc(l),zp_loc(l),order)

                              dom(ib)%wstar(i,j,k) = &
                                   dom(ib)%wstar(i,j,k) + dt * alfapr * Fpw(l) * delta * &
                                   Vp/Vcell

                           end do
                        end do
                     end do

                  end if

                  !$OMP END CRITICAL
               end if

               !     Actualizar velocidad paso previo
               if    (abs(up_pt(l))>10.0d0*abs(uop_pt(l))) then
               else if (abs(vp_pt(l))>10.0d0*abs(vop_pt(l))) then
               else if (abs(wp_pt(l))>10.0d0*abs(wop_pt(l))) then
               end if
               uop_loc(l) = up_pt(l)
               vop_loc(l) = vp_pt(l)
               wop_loc(l) = wp_pt(l)

               !     Actualizar posicion de particula
               xp_loc(l)=xp_loc(l)+up_pt(l)*dt
               yp_loc(l)=yp_loc(l)+vp_pt(l)*dt
               zp_loc(l)=zp_loc(l)+wp_pt(l)*dt

            end if   !if the particle belongs to the block

         end do  !end of loop in particles
         !$OMP ENDDO
         !$OMP END PARALLEL


      end do      !end loop in domains


      deallocate (ui_pt,vi_pt,wi_pt,uoi_pt,voi_pt,woi_pt)
      deallocate (ip,jp,kp,ipu,jpv,kpw)
      deallocate (up_pt,vp_pt,wp_pt)
      deallocate (id)

      return
    end subroutine particle_tracking

    !##########################################################################
    subroutine final_LPT(Wop_pt, vop_pt, ptsinproc, rhop_loc, rho_pt, xp_pt,yp_pt,zp_pt, &
         dp_pt, np_loc, Fu,Fv,Fw, xp_loc,yp_loc,zp_loc, uop_loc,vop_loc,wop_loc, dp_loc, &
         uop_pt, Fpu,Fpv,Fpw)
    !     sends backp(l) to master processor
    !#########################################################################

      real(dp), dimension(:), intent(in) :: wop_pt, vop_pt
      integer, dimension(:), intent(in) ::    ptsinproc
      real(dp), allocatable, dimension(:), intent(inout) :: rhop_loc
      real(dp), dimension(:), intent(in) :: rho_pt
      real(dp), dimension(:), intent(in) :: xp_pt,yp_pt,zp_pt, uop_pt
      real(dp), dimension(:), intent(in) :: dp_pt
      integer , intent(in) :: np_loc
      real(dp), dimension(:), intent(in) :: Fu,Fv,Fw
      real(dp), allocatable, dimension(:), intent(inout) :: xp_loc,yp_loc,zp_loc
      real(dp), allocatable, dimension(:), intent(inout) :: uop_loc,vop_loc,wop_loc, dp_loc
      real(dp), allocatable, dimension(:), intent(inout) :: Fpu,Fpv,Fpw
      integer,dimension(nprocs) :: strider
      integer :: s

      strider(1) = 0
      do s=2,nprocs
         strider(s) = ptsinproc(s-1) + strider(s-1)
      end do

      call MPI_BARRIER (MPI_COMM_WORLD,ierr)

      call MPI_GATHERV(xp_loc,np_loc,MPI_DOUBLE_PRECISION,xp_pt &
           ,ptsinproc,strider,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_GATHERV(yp_loc,np_loc,MPI_DOUBLE_PRECISION,yp_pt &
           ,ptsinproc,strider,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_GATHERV(zp_loc,np_loc,MPI_DOUBLE_PRECISION,zp_pt &
           ,ptsinproc,strider,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      call MPI_GATHERV(uop_loc,np_loc,MPI_DOUBLE_PRECISION,uop_pt &
           ,ptsinproc,strider,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_GATHERV(vop_loc,np_loc,MPI_DOUBLE_PRECISION,vop_pt &
           ,ptsinproc,strider,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_GATHERV(wop_loc,np_loc,MPI_DOUBLE_PRECISION,wop_pt &
           ,ptsinproc,strider,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      call MPI_GATHERV(Fpu,np_loc,MPI_DOUBLE_PRECISION,Fu &
           ,ptsinproc,strider,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_GATHERV(Fpv,np_loc,MPI_DOUBLE_PRECISION,Fv &
           ,ptsinproc,strider,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_GATHERV(Fpw,np_loc,MPI_DOUBLE_PRECISION,Fw &
           ,ptsinproc,strider,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_GATHERV(dp_loc,np_loc,MPI_DOUBLE_PRECISION,dp_pt &
           ,ptsinproc,strider,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_GATHERV(rhop_loc,np_loc,MPI_DOUBLE_PRECISION,rho_pt &
           ,ptsinproc,strider,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      if (np_loc>0) then
         deallocate (xp_loc,yp_loc,zp_loc)
         deallocate (uop_loc,vop_loc,wop_loc)
         deallocate (Fpu,Fpv,Fpw)
         deallocate (dp_loc,rhop_loc)
      end if

      return
    end subroutine final_LPT

end module multiflow3d_LPT

