!======================================================================!
!                 LAGRANGIAN PARTICLE TRACKING                         !
!----------------------------------------------------------------------!
!                          Boyang Chen                                 !
!                          Bruño Fraga                                 !
!                   Uni of Birmingham 2018-2025                        !
!======================================================================!
!######################################################################!
module multiflow3d_collison
  use, intrinsic :: iso_fortran_env, only: dp => real64

  implicit none
  private

  public :: collision_particle, collision_walls

contains

  subroutine collision_particle(l,xpg_loc,ypg_loc,zpg_loc,vopg_loc, wopg_loc,rhopg_loc,&
                                xp_loc,yp_loc,zp_loc,uop_loc,vop_loc, wop_loc,id,rhop_loc, &
                                dp_loc,up_pt,vp_pt,wp_pt,np_loc,npg_loc,uopg_loc,dpg_loc,k_n,dt,nbp,&
                                dom_id)
  !     Soft-sphere collision model                                    !
  !!###################################################################!

    !NOTE: make Lcol and Lcolwalls arrays for every fraction
    real(dp), allocatable, dimension(:),intent(inout):: xpg_loc,ypg_loc,zpg_loc, &
                                          vopg_loc, wopg_loc,rhopg_loc,dpg_loc,uopg_loc
    real(dp), allocatable, dimension(:),intent(in):: xp_loc,yp_loc,zp_loc,&
                                          uop_loc,vop_loc, wop_loc
    integer,allocatable,dimension(:),intent(in)::  id
    real(dp), allocatable, dimension(:),intent(in):: rhop_loc
    real(dp) ,intent(in):: dt
    integer ,intent(in):: nbp
    integer,allocatable,dimension(:) ,intent(in):: dom_id

    integer :: tot_np,ib
    integer,intent(in) :: l
    integer :: l2,ls
    real(dp) :: dis_x,dis_y,dis_z,dis_dd
    real(dp) :: dif_uvw,dis_xyz
    real(dp) :: dif1_uvw
    real(dp) :: dif2_uvw
    real(dp) :: overlap
    real(dp) :: lambda_p
    real(dp),allocatable,dimension(:):: up_sv,vp_sv,wp_sv
    real(dp),allocatable,dimension(:):: xp_sv,yp_sv,zp_sv
    real(dp),allocatable,dimension(:):: upg_sv,vpg_sv,wpg_sv
    real(dp),allocatable,dimension(:):: xpg_sv,ypg_sv,zpg_sv
    real(dp),allocatable,dimension(:):: dp_sv,dpg_sv
    real(dp) :: collision_t
    real(dp) :: collision_x,collision_y,collision_z
    real(dp) :: vector_x,vector_y,vector_z
    real(dp) :: collision_tx,collision_ty,collision_tz
    real(dp), allocatable, dimension(:) ,intent(in):: dp_loc
    real(dp), allocatable, dimension(:),intent(inout):: up_pt,vp_pt,wp_pt
    integer,intent(in) :: np_loc,npg_loc
    real(dp),intent(in) :: k_n

    real(dp) :: theta_col,e_col,mp

    !2. Damping
    e_col=1.0d0
    mp=rhop_loc(l)*(4/3)*3.1416_dp*(0.5_dp*dp_loc(l))**3
    theta_col=-2*log(e_col)*(mp*k_n)**0.5_dp/ &
              (3.1416_dp**2.0_dp+(log(e_col))**2.0_dp)

    tot_np = np_loc+npg_loc

    allocate(xp_sv(np_loc),yp_sv(np_loc),zp_sv(np_loc))
    allocate(up_sv(np_loc),vp_sv(np_loc),wp_sv(np_loc))
    allocate(xpg_sv(tot_np),ypg_sv(tot_np))
    allocate(upg_sv(tot_np),vpg_sv(tot_np))
    allocate(zpg_sv(tot_np),wpg_sv(tot_np))
    allocate(dp_sv(tot_np),dpg_sv(tot_np))

    xp_sv = 0.0_dp
    zp_sv = 0.0_dp
    zp_sv = 0.0_dp
    up_sv = 0.0_dp
    vp_sv = 0.0_dp
    wp_sv = 0.0_dp
    xpg_sv = 0.0_dp
    zpg_sv = 0.0_dp
    zpg_sv = 0.0_dp
    upg_sv = 0.0_dp
    vpg_sv = 0.0_dp
    wpg_sv = 0.0_dp
    dp_sv = 0.0_dp
    dpg_sv = 0.0_dp

    do ls=1,np_loc                ! save real particles
       xp_sv(ls) = xp_loc(ls)
       yp_sv(ls) = yp_loc(ls)
       zp_sv(ls) = zp_loc(ls)
       up_sv(ls) = uop_loc(ls)
       vp_sv(ls) = vop_loc(ls)
       wp_sv(ls) = wop_loc(ls)
       dp_sv(ls) = dp_loc(ls)
       xpg_sv(ls) = xp_loc(ls)
       ypg_sv(ls) = yp_loc(ls)
       zpg_sv(ls) = zp_loc(ls)
       upg_sv(ls) = uop_loc(ls)
       vpg_sv(ls) = vop_loc(ls)
       wpg_sv(ls) = wop_loc(ls)
       dpg_sv(ls) = dp_loc(ls)
    end do
    if (npg_loc>0) then
       do ls=1,npg_loc               ! save ghost particles
          xpg_sv(ls+np_loc) = xpg_loc(ls)
          ypg_sv(ls+np_loc) = ypg_loc(ls)
          zpg_sv(ls+np_loc) = zpg_loc(ls)
          upg_sv(ls+np_loc) = uopg_loc(ls)
          vpg_sv(ls+np_loc) = vopg_loc(ls)
          wpg_sv(ls+np_loc) = wopg_loc(ls)
          dpg_sv(ls+np_loc) = dpg_loc(ls)
       end do
    end if

    do ib=1,nbp
       ! ====================> p2p collision
       if (id(l)==dom_id(ib)) then
          do l2 = 1,tot_np
             dis_x = xpg_sv(l2)-xp_sv(l)                           ! difference on coordinate in x
             dis_y = ypg_sv(l2)-yp_sv(l)                           ! difference on coordinate in y
             dis_z = zpg_sv(l2)-zp_sv(l)                           ! difference on coordinate in z
             dis_dd = (dp_loc(l)+dpg_sv(l2))*0.5_dp                   ! sum up Radius
             dis_xyz = sqrt(dis_x**2+dis_y**2+dis_z**2)
             lambda_p = 0.375_dp*0.2_dp*(dp_sv(l)*0.5_dp+dpg_sv(l2)*0.5_dp)
             !  CFL particle-particle
             if ((dis_xyz/=0.0d0).and.(dis_xyz<(dis_dd+lambda_p))) then          !

                dif1_uvw = up_sv(l)*dis_x/dis_xyz+vp_sv(l)*dis_y/dis_xyz &
                           +wp_sv(l)*dis_z/dis_xyz
                dif2_uvw = upg_sv(l2)*dis_x/dis_xyz+vpg_sv(l2)*dis_y/dis_xyz &
                           +wpg_sv(l2)*dis_z/dis_xyz

                dif_uvw = dif1_uvw - dif2_uvw                     ! difference on velocity(vector)
                overlap = MAX((dis_dd-abs(dis_xyz)),0.0d0)
                collision_x = -k_n * overlap * dis_x/dis_xyz &
                              - theta_col * dif_uvw * dis_x/dis_xyz
                collision_y = -k_n * overlap * dis_y/dis_xyz &
                              - theta_col * dif_uvw * dis_y/dis_xyz
                collision_z = -k_n * overlap * dis_z/dis_xyz &
                              - theta_col * dif_uvw * dis_z/dis_xyz

                !---------------------- Particle rotation ----------------------------------------

                collision_t = -0.1_dp*sqrt(collision_x**2+ &
                             collision_y**2+collision_z**2)   !  uf*abs(Fcoln)

                vector_x = (up_sv(l)-upg_sv(l2))*(1-dis_x**2/dis_xyz**2)
                vector_y = (vp_sv(l)-vpg_sv(l2))*(1-dis_y**2/dis_xyz**2)
                vector_z = (wp_sv(l)-wpg_sv(l2))*(1-dis_z**2/dis_xyz**2)

                collision_tx = collision_t*vector_x/(sqrt(vector_x**2+ &
                               vector_y**2+vector_z**2)+1d-12)
                collision_ty = collision_t*vector_y/(sqrt(vector_x**2+ &
                               vector_y**2+vector_z**2)+1d-12)
                collision_tz = collision_t*vector_z/(sqrt(vector_x**2+ &
                               vector_y**2+vector_z**2)+1d-12)

                ! --------------------------------------------------------------------------------
                ! breakdown collision force

                up_pt(l) = up_pt(l)+dt*(collision_x+collision_tx)/mp
                vp_pt(l) = vp_pt(l)+dt*(collision_y+collision_ty)/mp
                wp_pt(l) = wp_pt(l)+dt*(collision_z+collision_tz)/mp
             end if                         ! MPI block
          end do                         ! end search loop
       end if                         ! distance
    end do                         ! end p2p loop

    if (npg_loc>0) then
       deallocate (xpg_loc,ypg_loc,zpg_loc)
       deallocate (uopg_loc,vopg_loc,wopg_loc)
       deallocate (dpg_loc,rhopg_loc)
    end if

    return
  end subroutine collision_particle

  !######################################################################!
  subroutine collision_walls(l,rhop_loc,yp_loc,zp_loc,dp_loc,&
                             up_pt,vp_pt,wp_pt,k_n,yst,yen,zst,zen,bc_s,bc_n,bc_b,bc_t,dt)   !
  !     Calculates collisions with walls and boundaries                  !
  !######################################################################!
    integer,intent(in) :: l
    real(dp) :: fcol_n,fcol_t,mu_f
    real(dp) :: lambda_w,lambda_u,lambda_v
    real(dp) :: theta_col,e_col,mp  !,k_t
    real(dp) :: deltap

    real(dp), allocatable, dimension(:),intent(in):: rhop_loc
    real(dp), allocatable, dimension(:),intent(in):: yp_loc,zp_loc
    real(dp), allocatable, dimension(:),intent(in):: dp_loc
    real(dp), allocatable, dimension(:),intent(inout):: up_pt,vp_pt,wp_pt
    real(dp),intent(in) :: k_n
    real(dp),intent(in) :: yst,yen,zst,zen
    integer,intent(in) :: bc_s,bc_n,bc_b,bc_t
    real(dp),intent(in) :: dt

    mu_f=9.2d-2

    !1.Define force range
    lambda_u=0.75_dp*up_pt(l)*dt
    lambda_v=0.75_dp*vp_pt(l)*dt
    lambda_w=0.75_dp*wp_pt(l)*dt

    !2.0_dp Spring stiffness

    !3. Damping
    e_col=1.0d0
    mp=rhop_loc(l)*(4/3)*3.1416_dp*(0.5_dp*dp_loc(l))**3
    theta_col=-2*log(e_col)*(mp*k_n)**0.5_dp/ &
               (3.1416_dp**2.0_dp+(log(e_col))**2.0_dp)

    ! ----------------------- collisions with bottom wall ----------------------------------
    if (zp_loc(l)<lambda_w+0.5_dp*dp_loc(l)) then

       !a. overlap
       deltap=max((zp_loc(l)-dp_loc(l)/2)-zst,0.0d0)
       !b. normal force
       fcol_n=-k_n*deltap-theta_col*wp_pt(l)
       wp_pt(l) = wp_pt(l) + dt*fcol_n/mp
       !c. tangential force
       fcol_t=mu_f*fcol_n
       if (bc_b/=3) then                 !slip condition
          up_pt(l) = up_pt(l) + dt*fcol_t/mp
          vp_pt(l) = vp_pt(l) + dt*fcol_t/mp
       end if

    end if

    ! ----------------------- collisions with top wall ----------------------------------
    if (zp_loc(l)>zen-(lambda_w+0.5_dp*dp_loc(l))) then

       !a. overlap
       deltap=max(zp_loc(l)-(zen+0.5_dp*dp_loc(l)),0.0d0)
       !b. normal force
       fcol_n=-k_n*deltap-theta_col*wp_pt(l)
       wp_pt(l) = wp_pt(l) + dt*fcol_n/mp
       !c. tangential force
       fcol_t=mu_f*fcol_n
       if (bc_t/=3) then                 !slip condition
          up_pt(l) = up_pt(l) + dt*fcol_t/mp
          vp_pt(l) = vp_pt(l) + dt*fcol_t/mp
       end if

    end if

    ! ----------------------- collisions with south wall ----------------------------------
    if (yp_loc(l)<lambda_v+0.5_dp*dp_loc(l)) then

       !a. overlap
       deltap=max((yp_loc(l)-dp_loc(l)/2)-yst,0.0d0)
       !b. normal force
       fcol_n=-k_n*deltap-theta_col*vp_pt(l)
       vp_pt(l) = vp_pt(l) + dt*fcol_n/mp
       !c. tangential force
       fcol_t=mu_f*fcol_n
       if (bc_s/=3) then                 !slip condition
          up_pt(l) = up_pt(l) + dt*fcol_t/mp
          wp_pt(l) = wp_pt(l) + dt*fcol_t/mp
       end if
    end if

    ! ----------------------- collisions with north wall ----------------------------------
    if (yp_loc(l)>yen-(lambda_v+0.5_dp*dp_loc(l))) then

       !a. overlap
       deltap=max(yp_loc(l)-(yen+0.5_dp*dp_loc(l)),0.0d0)
       !b. normal force
       fcol_n=-k_n*deltap-theta_col*vp_pt(l)
       vp_pt(l) = vp_pt(l) + dt*fcol_n/mp

       !c. tangential force
       fcol_t=mu_f*fcol_n

       if (bc_n/=3) then                 !slip condition
          up_pt(l) = up_pt(l) + dt*fcol_t/mp
          wp_pt(l) = wp_pt(l) + dt*fcol_t/mp
       end if
    end if

  end subroutine collision_walls

end module multiflow3d_collison
