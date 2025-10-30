!======================================================================!
!                 LAGRANGIAN PARTICLE TRACKING                         !
!----------------------------------------------------------------------!
!                          Boyang Chen                                 ! 
!                          Bruño Fraga                                 !
!                   Uni of Birmingham 2018-2024                        !
!======================================================================!
C !######################################################################!
C       Double precision function collision_particle(l)                   
C !     Soft-sphere collision model                                      !
C !######################################################################!

C       !NOTE: make Lcol and Lcolwalls arrays for every fraction 

C       use multidata
C       use mpi
C       use vars   
C       use vars_pt

C       implicit none

C       integer tot_np,ib
C       integer l,l1,l2,ls,ls1,ls2
C       double precision dis_x,dis_y,dis_z,dis_dd,max_dis
C       double precision cita_xy,cita_xyz,dis_xy,dif_uvw,dis_xyz
C       double precision dif1_uv,dis1_xy,dif1_uvw,dis1_xyz
C       double precision dif2_uv,dis2_xy,dif2_uvw,dis2_xyz
C       double precision ip1,jp1,kp1,ip2,jp2,kp2
C       double precision damp,stiffness,collision,overlap
C       double precision lambda_p,lambda_wb,lambda_ww,lambda_wt
C       double precision lambda_we,lambda_ws,lambda_wn
C       double precision,allocatable,dimension(:):: up_sv,vp_sv,wp_sv
C       double precision,allocatable,dimension(:):: xp_sv,yp_sv,zp_sv 
C       double precision,allocatable,dimension(:):: upg_sv,vpg_sv,wpg_sv
C       double precision,allocatable,dimension(:):: xpg_sv,ypg_sv,zpg_sv 
C       double precision,allocatable,dimension(:):: dp_sv,dpg_sv
C       double precision dif1_uvw_t,dif2_uvw_t,collision_t
C       double precision collision_x,collision_y,collision_z
C       double precision vector_x,vector_y,vector_z
C       double precision collision_tx,collision_ty,collision_tz
C       logical,allocatable,dimension(:):: collide_pt

C       tot_np = np_loc+npg_loc

C       allocate(xp_sv(np_loc),yp_sv(np_loc),zp_sv(np_loc))
C       allocate(up_sv(np_loc),vp_sv(np_loc),wp_sv(np_loc))
C       allocate(xpg_sv(tot_np),ypg_sv(tot_np))
C       allocate(upg_sv(tot_np),vpg_sv(tot_np))
C       allocate(zpg_sv(tot_np),wpg_sv(tot_np))
C       allocate(dp_sv(tot_np),dpg_sv(tot_np))

C       xp_sv = 0.0 ; zp_sv = 0.0 ; zp_sv = 0.0
C       up_sv = 0.0 ; vp_sv = 0.0 ; wp_sv = 0.0
C       xpg_sv = 0.0 ; zpg_sv = 0.0 ; zpg_sv = 0.0
C       upg_sv = 0.0 ; vpg_sv = 0.0 ; wpg_sv = 0.0
C       dp_sv = 0.0 ; dpg_sv = 0.0 

C       do ls=1,np_loc                ! save real + ghost particles
C              xp_sv(ls) = xp_loc(ls)
C              yp_sv(ls) = yp_loc(ls)
C              zp_sv(ls) = zp_loc(ls)
C              up_sv(ls) = uop_loc(ls)
C              vp_sv(ls) = vop_loc(ls)
C              wp_sv(ls) = wop_loc(ls)
C              dp_sv(ls) = dp_loc(ls)
C              xpg_sv(ls) = xp_loc(ls)
C              ypg_sv(ls) = yp_loc(ls)
C              zpg_sv(ls) = zp_loc(ls)
C              upg_sv(ls) = uop_loc(ls)
C              vpg_sv(ls) = vop_loc(ls)
C              wpg_sv(ls) = wop_loc(ls)
C              dpg_sv(ls) = dp_loc(ls)
C       enddo
C       do ls=1,npg_loc               ! save ghost particles+bubbles
C              xpg_sv(ls+np_loc) = xpg_loc(ls)
C              ypg_sv(ls+np_loc) = ypg_loc(ls)
C              zpg_sv(ls+np_loc) = zpg_loc(ls)
C              upg_sv(ls+np_loc) = uopg_loc(ls)
C              vpg_sv(ls+np_loc) = vopg_loc(ls)
C              wpg_sv(ls+np_loc) = wopg_loc(ls)
C              dpg_sv(ls+np_loc) = dp_loc(ls)
C       enddo

C       do ib=1,nbp
C ! ====================> collision looping for particle-particle
C       IF (id(l).eq.dom_id(ib)) then 
C       do l2 = 1,tot_np
C             dis_x = xpg_sv(l2)-xp_sv(l)                           ! difference on coordinate in x
C             dis_y = ypg_sv(l2)-yp_sv(l)                           ! difference on coordinate in y
C             dis_z = zpg_sv(l2)-zp_sv(l)                           ! difference on coordinate in z
C             dis_dd = (dp_loc(l)+dpg_sv(l2))*0.5                   ! sum up Radius 
C             dis_xyz = sqrt(dis_x**2+dis_y**2+dis_z**2)
C             lambda_p = 0.375*0.2*(dp_sv(l)*0.5+dpg_sv(l2)*0.5)          ! CFL for particle-particle
C       if ((dis_xyz.ne.0.d0).and.(dis_xyz.lt.(dis_dd+lambda_p))) then          ! 
            
C !           write(myrank+700,*) l1,zp_sv(l1),wp_sv(l1)
C             dif1_uvw = up_sv(l)*dis_x/dis_xyz+vp_sv(l)*dis_y/dis_xyz
C      &      +wp_sv(l)*dis_z/dis_xyz
C             dif2_uvw = upg_sv(l2)*dis_x/dis_xyz+vpg_sv(l2)*dis_y/dis_xyz
C      &      +wpg_sv(l2)*dis_z/dis_xyz           
                                                                              
C             dif_uvw = dif1_uvw - dif2_uvw                               ! difference on velocity(vector)    
C             overlap = MAX((dis_dd-abs(dis_xyz)),0.d0)
C             collision_x = -ea_k * overlap * dis_x/dis_xyz
C      & - ea_yita * dif_uvw * dis_x/dis_xyz
C             collision_y = -ea_k * overlap * dis_y/dis_xyz
C      & - ea_yita * dif_uvw * dis_y/dis_xyz
C             collision_z = -ea_k * overlap * dis_z/dis_xyz
C      & - ea_yita * dif_uvw * dis_z/dis_xyz
C !           write(myrank+800,*) collision_x,collision_y,collision_z
C !---------------------- Particle rotation ----------------------------------------
            
C c           dif1_uvw_t = (up_sv(l1)*dis_x/dis_xyz+vp_sv(l1)*dis_y/dis_xyz)
C c     &  *dis_z/dis_xyz+wp_sv(l1)*sqrt(dis_x**2+dis_y**2)/dis_xyz
C c           dif2_uvw_t = (up_sv(l2)*dis_x/dis_xyz+vp_sv(l2)*dis_y/dis_xyz)
C c     &  *dis_z/dis_xyz+wp_sv(l2)*sqrt(dis_x**2+dis_y**2)/dis_xyz
C             collision_t = -0.1*sqrt(collision_x**2+
C      &                  collision_y**2+collision_z**2)   !  uf*abs(Fcoln)
            
C             vector_x = (up_sv(l)-upg_sv(l2))*(1-dis_x**2/dis_xyz**2)
C             vector_y = (vp_sv(l)-vpg_sv(l2))*(1-dis_y**2/dis_xyz**2)
C             vector_z = (wp_sv(l)-wpg_sv(l2))*(1-dis_z**2/dis_xyz**2)
C !           write(myrank+800,*) vector_x,vector_y,vector_z
C             collision_tx = collision_t*vector_x/(sqrt(vector_x**2+
C      &  vector_y**2+vector_z**2)+1d-12) 
C             collision_ty = collision_t*vector_y/(sqrt(vector_x**2+
C      &  vector_y**2+vector_z**2)+1d-12)
C             collision_tz = collision_t*vector_z/(sqrt(vector_x**2+
C      &  vector_y**2+vector_z**2)+1d-12) 
C !           write(myrank+800,*) collision_tz,collision_ty,collision_tx,collision_t              
C ! --------------------------------------------------------------------------------
C             wp_pt(l) = wp_pt(l)+dt*(collision_z+collision_tz)                 ! breakdown collision force         
      
C             up_pt(l) = up_pt(l)+dt*(collision_x+collision_tx) 
      
C             vp_pt(l) = vp_pt(l)+dt*(collision_y+collision_ty)
C             endif                         ! end yourself looping 
C c     endif 
C       enddo                         ! end loop for search
C c     endif 
C       ENDIF
C c     enddo     ! end loop for particle-particle collision 

C       enddo

C       if (npg_loc.gt.0) then
C       deallocate (xpg_loc,ypg_loc,zpg_loc)
C       deallocate (uopg_loc,vopg_loc,wopg_loc)
C       deallocate (dpg_loc)
C       endif

C       return
C       end function

!######################################################################!      
      Subroutine collision_walls(l,ib)                !
!     Calculates collisions with walls and boundaries                  !
!######################################################################!
      use multidata
      use mpi
      use vars   
      use vars_pt

      implicit none

      integer l,ib
      double precision damp,stiffness,collision
      double precision lambda_p,lambda_wb,lambda_ww,lambda_wt
      double precision lambda_we,lambda_ws,lambda_wn
      double precision :: e_k,e_yita,ea_k,ea_yita,e_col



!     compute collision parameters
      e_col=1.d0
      e_k = 0.5*rhop_loc(l)*4/3*3.1416*(0.5*dp_loc(l))**3 /                          ! equivelant mass             
     &            (15*dt)**2*(3.1416**2+alog(e_col**2))     
      e_yita = -2*alog(e_col)*sqrt(0.5*rhop_loc(l)*4/3*
     &        3.1416*(0.5*dp_loc(l))**3*e_k) / 
     &            (3.1416**2+alog(e_col**2))

      ea_k = e_k / (rhop_loc(l)*4/3*3.1416*(0.5*dp_loc(l))**3)
      ea_yita = e_yita / (rhop_loc(l)*4/3*3.1416*(0.5*dp_loc(l))**3)

!! ====================> collision loop for particle-wall soft!!!!!!!
!     radius of influence lambda

      lambda_wb = 0.75*0.1*dp_loc(l)+dp_loc(l)*0.5
      lambda_wt = zen-dp_loc(l)*0.5-0.75*0.1*dp_loc(l)
      lambda_ww = 0.75*0.1*dp_loc(l)+dp_loc(l)*0.5
      lambda_we = xen-dp_loc(l)*0.5-0.75*0.1*dp_loc(l)
      lambda_ws = 0.75*0.1*dp_loc(l)+dp_loc(l)*0.5
      lambda_wn = yen-dp_loc(l)*0.5-0.75*0.1*dp_loc(l)      
! ----------------------- collisions with bottom wall ----------------------------------                          
      if ((zp_loc(l).lt.lambda_wb).and.(wp_pt(l).lt.0)) then
            damp = 2*ea_yita * wp_pt(l)
            stiffness = -2*ea_k * MAX((dp_loc(l)*0.5-zp_loc(l)),0.d0)
            collision = - stiffness - damp

            write (6,*)l,wp_pt(l),stiffness,damp

            wp_pt(l) = wp_pt(l) + dt*collision 

                        write (6,*)l,wp_pt(l),collision

      endif 
! -------------------------------------------------------------------------------------
! ----------------------- collisions with top wall -----------------------------------
        if ((zp_loc(l).gt.lambda_wt).and.(wp_pt(l).gt.0)) then
            damp = 2*ea_yita * wp_pt(l)
            stiffness = 2*ea_k * MAX((dp_loc(l)*0.5+zp_loc(l)-zen),0.d0)

            if (dom(ib)%bc_top.eq.3) stiffness = 0

            collision = - stiffness - damp
            wp_pt(l) = wp_pt(l) + dt*collision
        endif
! -------------------------------------------------------------------------------------
! ----------------------- collisions with west wall -----------------------------------
      If (.not.PERIODIC) then
      if ((xp_loc(l).lt.lambda_ww).and.(up_pt(l).lt.0)) then
            damp = 2*ea_yita * up_pt(l)
            stiffness = -2*ea_k * MAX((dp_loc(l)*0.5-xp_loc(l)),0.d0)
            collision = - stiffness - damp
            up_pt(l) = up_pt(l) + dt*collision 
      endif 
      ENDIF
! -------------------------------------------------------------------------------------
! ----------------------- collisions with east wall -----------------------------------
      IF (.not.PERIODIC) then
      if ((xp_loc(l).gt.lambda_we).and.(up_pt(l).gt.0)) then
            damp = 2*ea_yita * up_pt(l)
            stiffness = 2*ea_k * MAX((dp_loc(l)*0.5+xp_loc(l)-xen),0.d0)
            collision = - stiffness - damp
            up_pt(l) = up_pt(l) + dt*collision 
      endif 
      ENDIF
! -------------------------------------------------------------------------------------
! ----------------------- collisions with south wall ----------------------------------
      if ((yp_loc(l).lt.lambda_ws).and.(vp_pt(l).lt.0)) then
            damp = 2*ea_yita * vp_pt(l)
            stiffness = -2*ea_k * MAX((dp_loc(l)*0.5-yp_loc(l)),0.d0)
            collision = - stiffness - damp
            vp_pt(l) = vp_pt(l) + dt*collision 
      endif 
! -------------------------------------------------------------------------------------
! ----------------------- collisions with north wall ----------------------------------
      if ((yp_loc(l).gt.lambda_wn).and.(vp_pt(l).gt.0)) then
            damp = 2*ea_yita * vp_pt(l)
            stiffness = 2*ea_k * MAX((dp_loc(l)*0.5+yp_loc(l)-yen),0.d0)
            collision = - stiffness - damp
            vp_pt(l) = vp_pt(l) + dt*collision 
      endif 
!=============================================================================

      
!     Actualizar velocidad paso previo
            uop_loc(l) = up_pt(l)
            vop_loc(l) = vp_pt(l)
            wop_loc(l) = wp_pt(l)
!     Actualizar posicion de particula
            xp_loc(l)=xp_loc(l)+up_pt(l)*dt
            yp_loc(l)=yp_loc(l)+vp_pt(l)*dt
            zp_loc(l)=zp_loc(l)+wp_pt(l)*dt

      return
      end