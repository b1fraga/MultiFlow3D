!======================================================================!
!                 LAGRANGIAN PARTICLE TRACKING                         !
!----------------------------------------------------------------------!
!                          Boyang Chen                                 ! 
!                          Bruño Fraga                                 !
!                   Uni of Birmingham 2018-2025                        !
!======================================================================!
C !######################################################################!
       Subroutine collision_particle(l)                   
      !     Soft-sphere collision model                                      !
C !######################################################################!

       !NOTE: make Lcol and Lcolwalls arrays for every fraction 

       use multidata
       use mpi
       use vars   
       use vars_pt

       implicit none

       integer tot_np,ib
       integer l,l1,l2,ls,ls1,ls2
       double precision dis_x,dis_y,dis_z,dis_dd,max_dis
       double precision cita_xy,cita_xyz,dis_xy,dif_uvw,dis_xyz
       double precision dif1_uv,dis1_xy,dif1_uvw,dis1_xyz
       double precision dif2_uv,dis2_xy,dif2_uvw,dis2_xyz
       double precision ip1,jp1,kp1,ip2,jp2,kp2
       double precision damp,stiffness,collision,overlap
       double precision lambda_p,lambda_wb,lambda_ww,lambda_wt
       double precision lambda_we,lambda_ws,lambda_wn
       double precision,allocatable,dimension(:):: up_sv,vp_sv,wp_sv
       double precision,allocatable,dimension(:):: xp_sv,yp_sv,zp_sv 
       double precision,allocatable,dimension(:):: upg_sv,vpg_sv,wpg_sv
       double precision,allocatable,dimension(:):: xpg_sv,ypg_sv,zpg_sv 
       double precision,allocatable,dimension(:):: dp_sv,dpg_sv
       double precision dif1_uvw_t,dif2_uvw_t,collision_t
       double precision collision_x,collision_y,collision_z
       double precision vector_x,vector_y,vector_z
       double precision collision_tx,collision_ty,collision_tz
       logical,allocatable,dimension(:):: collide_pt

       double precision k_n,theta_col,e_col,mp

      !1. Spring stiffness
       k_n=1.72d7
      
      !2. Damping
       e_col=1.d0
       mp=rhop_loc(l)*(4/3)*3.1416*(0.5*dp_loc(l))**3
       theta_col=-2*alog(e_col)*(mp*k_n)**0.5/
     &      (3.1416**2+(alog(e_col))**2)

       tot_np = np_loc+npg_loc

       allocate(xp_sv(np_loc),yp_sv(np_loc),zp_sv(np_loc))
       allocate(up_sv(np_loc),vp_sv(np_loc),wp_sv(np_loc))
       allocate(xpg_sv(tot_np),ypg_sv(tot_np))
       allocate(upg_sv(tot_np),vpg_sv(tot_np))
       allocate(zpg_sv(tot_np),wpg_sv(tot_np))
       allocate(dp_sv(tot_np),dpg_sv(tot_np))

       xp_sv = 0.0 ; zp_sv = 0.0 ; zp_sv = 0.0
       up_sv = 0.0 ; vp_sv = 0.0 ; wp_sv = 0.0
       xpg_sv = 0.0 ; zpg_sv = 0.0 ; zpg_sv = 0.0
       upg_sv = 0.0 ; vpg_sv = 0.0 ; wpg_sv = 0.0
       dp_sv = 0.0 ; dpg_sv = 0.0 

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
       enddo
       do ls=1,npg_loc               ! save ghost particles
              xpg_sv(ls+np_loc) = xpg_loc(ls)
              ypg_sv(ls+np_loc) = ypg_loc(ls)
              zpg_sv(ls+np_loc) = zpg_loc(ls)
              upg_sv(ls+np_loc) = uopg_loc(ls)
              vpg_sv(ls+np_loc) = vopg_loc(ls)
              wpg_sv(ls+np_loc) = wopg_loc(ls)
              dpg_sv(ls+np_loc) = dp_loc(ls)
       enddo

       do ib=1,nbp
 ! ====================> p2p collision
       IF (id(l).eq.dom_id(ib)) then 
       do l2 = 1,tot_np
             dis_x = xpg_sv(l2)-xp_sv(l)                           ! difference on coordinate in x
             dis_y = ypg_sv(l2)-yp_sv(l)                           ! difference on coordinate in y
             dis_z = zpg_sv(l2)-zp_sv(l)                           ! difference on coordinate in z
             dis_dd = (dp_loc(l)+dpg_sv(l2))*0.5                   ! sum up Radius 
             dis_xyz = sqrt(dis_x**2+dis_y**2+dis_z**2)
             lambda_p = 0.375*0.2*(dp_sv(l)*0.5+dpg_sv(l2)*0.5)          ! CFL for particle-particle
       if ((dis_xyz.ne.0.d0).and.(dis_xyz.lt.(dis_dd+lambda_p))) then          ! 
            
 !           write(myrank+700,*) l1,zp_sv(l1),wp_sv(l1)
       dif1_uvw = up_sv(l)*dis_x/dis_xyz+vp_sv(l)*dis_y/dis_xyz
     &      +wp_sv(l)*dis_z/dis_xyz
       dif2_uvw = upg_sv(l2)*dis_x/dis_xyz+vpg_sv(l2)*dis_y/dis_xyz
     &      +wpg_sv(l2)*dis_z/dis_xyz           
                                                                              
             dif_uvw = dif1_uvw - dif2_uvw                               ! difference on velocity(vector)    
             overlap = MAX((dis_dd-abs(dis_xyz)),0.d0)
             collision_x = -k_n * overlap * dis_x/dis_xyz
     & - theta_col * dif_uvw * dis_x/dis_xyz
             collision_y = -k_n * overlap * dis_y/dis_xyz
     & - theta_col * dif_uvw * dis_y/dis_xyz
             collision_z = -k_n * overlap * dis_z/dis_xyz
     & - theta_col * dif_uvw * dis_z/dis_xyz
 !           write(myrank+800,*) collision_x,collision_y,collision_z
 !---------------------- Particle rotation ----------------------------------------
            
c           dif1_uvw_t = (up_sv(l1)*dis_x/dis_xyz+vp_sv(l1)*dis_y/dis_xyz)
c     &  *dis_z/dis_xyz+wp_sv(l1)*sqrt(dis_x**2+dis_y**2)/dis_xyz
c           dif2_uvw_t = (up_sv(l2)*dis_x/dis_xyz+vp_sv(l2)*dis_y/dis_xyz)
c     &  *dis_z/dis_xyz+wp_sv(l2)*sqrt(dis_x**2+dis_y**2)/dis_xyz
             collision_t = -0.1*sqrt(collision_x**2+
     &      collision_y**2+collision_z**2)   !  uf*abs(Fcoln)
            
             vector_x = (up_sv(l)-upg_sv(l2))*(1-dis_x**2/dis_xyz**2)
             vector_y = (vp_sv(l)-vpg_sv(l2))*(1-dis_y**2/dis_xyz**2)
             vector_z = (wp_sv(l)-wpg_sv(l2))*(1-dis_z**2/dis_xyz**2)
 !           write(myrank+800,*) vector_x,vector_y,vector_z
             collision_tx = collision_t*vector_x/(sqrt(vector_x**2+
     &  vector_y**2+vector_z**2)+1d-12) 
             collision_ty = collision_t*vector_y/(sqrt(vector_x**2+
     &  vector_y**2+vector_z**2)+1d-12)
             collision_tz = collision_t*vector_z/(sqrt(vector_x**2+
     &  vector_y**2+vector_z**2)+1d-12) 
 !           write(myrank+800,*) collision_tz,collision_ty,collision_tx,collision_t              
 ! --------------------------------------------------------------------------------
 ! breakdown collision force         
      
       up_pt(l) = up_pt(l)+dt*(collision_x+collision_tx) 
       vp_pt(l) = vp_pt(l)+dt*(collision_y+collision_ty)
       wp_pt(l) = wp_pt(l)+dt*(collision_z+collision_tz)
       endif                         ! MPI block 
       enddo                         ! end search loop
       ENDIF                         ! distance
       enddo                         ! end p2p loop 

       if (npg_loc.gt.0) then
              deallocate (xpg_loc,ypg_loc,zpg_loc)
              deallocate (uopg_loc,vopg_loc,wopg_loc)
              deallocate (dpg_loc,rhopg_loc)
       endif

       return
       end subroutine

!######################################################################!      
      Subroutine collision_walls(l)                !
!     Calculates collisions with walls and boundaries                  !
!######################################################################!
      use multidata
      use mpi
      use vars   
      use vars_pt

      implicit none

      integer l
      double precision fcol_n,fcol_t,mu_f
      double precision lambda_w,lambda_u,lambda_v
      double precision k_n,theta_col,e_col,mp!,k_t
      double precision deltap

       mu_f=9.2d-2

      !1.Define force range
       lambda_u=0.75*up_pt(l)*dt  
       lambda_v=0.75*vp_pt(l)*dt
       lambda_w=0.75*wp_pt(l)*dt
      
      !2. Spring stiffness
       k_n=1.72d7
!       k_t=1.48d7  
      
      !3. Damping
       e_col=1.d0
       mp=rhop_loc(l)*(4/3)*3.1416*(0.5*dp_loc(l))**3
       theta_col=-2*alog(e_col)*(mp*k_n)**0.5/
     &      (3.1416**2+(alog(e_col))**2)
      
! ----------------------- collisions with bottom wall ----------------------------------                          
       if (zp_loc(l).lt.lambda_w+0.5*dp_loc(l)) then
      
      !a. overlap
       deltap=max((zp_loc(l)-dp_loc(l)/2)-zst,0.d0)
      !b. normal force
       fcol_n=-k_n*deltap-theta_col*wp_pt(l)
       wp_pt(l) = wp_pt(l) + dt*fcol_n/mp 
      !c. tangential force
       fcol_t=mu_f*fcol_n
       if (bc_b.ne.3) then                 !slip condition
            up_pt(l) = up_pt(l) + dt*fcol_t/mp 
            vp_pt(l) = vp_pt(l) + dt*fcol_t/mp 
       endif 
      !write (6,*)l,wp_pt(l),zp_loc(l),fcol_n,fcol_t
       endif 
! ----------------------- collisions with top wall ----------------------------------                          
       if (zp_loc(l).gt.zen-lambda_w+0.5*dp_loc(l)) then
      
            !a. overlap
            deltap=max(zp_loc(l)-zen+0.5*dp_loc(l),0.d0) 
            !b. normal force
            fcol_n=-k_n*deltap-theta_col*wp_pt(l)
            wp_pt(l) = wp_pt(l) + dt*fcol_n/mp 
            !c. tangential force
            fcol_t=mu_f*fcol_n
            if (bc_t.ne.3) then                 !slip condition
                  up_pt(l) = up_pt(l) + dt*fcol_t/mp 
                  vp_pt(l) = vp_pt(l) + dt*fcol_t/mp 
            endif 
            !write (6,*)l,wp_pt(l),zp_loc(l),fcol_n,fcol_t
            
       endif 
! ----------------------- collisions with south wall ----------------------------------                          
       if (yp_loc(l).lt.lambda_v+0.5*dp_loc(l)) then
      
      !a. overlap
       deltap=max((yp_loc(l)-dp_loc(l)/2)-yst,0.d0)
      !b. normal force
       fcol_n=-k_n*deltap-theta_col*vp_pt(l)
       vp_pt(l) = vp_pt(l) + dt*fcol_n/mp 
      !c. tangential force
       fcol_t=mu_f*fcol_n
       if (bc_s.ne.3) then                 !slip condition
            up_pt(l) = up_pt(l) + dt*fcol_t/mp 
            wp_pt(l) = wp_pt(l) + dt*fcol_t/mp 
       endif             
       endif 
! ----------------------- collisions with north wall ----------------------------------                          
       if (yp_loc(l).gt.yen-lambda_v+0.5*dp_loc(l)) then
      
!a. overlap
       deltap=max(yp_loc(l)-yen+0.5*dp_loc(l),0.d0) 
!b. normal force
       fcol_n=-k_n*deltap-theta_col*vp_pt(l)
       vp_pt(l) = vp_pt(l) + dt*fcol_n/mp 
            
!c. tangential force
       fcol_t=mu_f*fcol_n
            
       if (bc_n.ne.3) then                 !slip condition
            up_pt(l) = up_pt(l) + dt*fcol_t/mp 
            wp_pt(l) = wp_pt(l) + dt*fcol_t/mp 
       endif 
       endif 
       return
       end