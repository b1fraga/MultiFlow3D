!!=======================================================================
!                 LAGRANGIAN PARTICLE TRACKING
!                 Bruño Fraga Bugallo
!                 Cardiff 2013-2017
!                 Stanford 2018
!                 Birmingham 2019-2023
!=======================================================================
!##########################################################################
        subroutine alloc_pt
!##########################################################################
      use vars
      use mpi
      use multidata
      use vars_pt
      implicit none

      integer l,out_cnt   
      logical,allocatable,dimension(:):: out_pt             

!     wtime_refresh = MPI_WTIME ( ) 

      allocate(out_pt(np))

      out_pt=.false.

      out_cnt = 0
      
!     do ib=1,nbp
!           is = dom(ib)%isp
!           ie = dom(ib)%iep
!           js = dom(ib)%jsp
!           je = dom(ib)%jep
!           ks = dom(ib)%ksp
!           ke = dom(ib)%kep

      do l=1,np

      IF (Lcolwall) then

 !================================== collision at wall ========================================
            if (yp_pt(l).lt.yst) then
                  yp_pt(l)=yst+0.75*dp_pt(l)
            elseif (yp_pt(l).gt.yen) then
                  yp_pt(l)=yen-0.75*dp_pt(l)
            endif
            if (zp_pt(l).lt.zst) then
                  zp_pt(l)=zst+0.75*dp_pt(l)
            elseif (zp_pt(l).gt.zen) then
                  zp_pt(l)=zen-0.75*dp_pt(l)
            endif
! =============================================================================================

      ENDIF

!     Comprobar si permanece en el dominio
            if ((xp_pt(l).le.xst).or.(xp_pt(l)
     &.ge.xen)) then
                  out_pt(l) = .TRUE.
                  out_cnt = out_cnt + 1
                  goto 50
            elseif ((yp_pt(l).le.yst).or.(yp_pt(l)
     &.ge.yen)) then
                  out_pt(l) = .TRUE.
                  out_cnt = out_cnt + 1
                  goto 50
            elseif ((zp_pt(l).lt.zst).or.(zp_pt(l)
     &.ge.zen)) then
                  out_pt(l) = .TRUE.
                  out_cnt = out_cnt + 1
                  goto 50
            end if

50    continue


            if (.not.out_pt(l)) then
            xpold(l-out_cnt)=xp_pt(l)
            ypold(l-out_cnt)=yp_pt(l)
            zpold(l-out_cnt)=zp_pt(l)
            uopold(l-out_cnt)=uop_pt(l)
            vopold(l-out_cnt)=vop_pt(l)
            wopold(l-out_cnt)=wop_pt(l)
            dp_old(l-out_cnt)=dp_pt(l)
            rhop_old(l-out_cnt)=rho_pt(l)
            endif

      enddo

      np = np - out_cnt

      if (out_cnt.gt.0) then
            write(202,*) ntime,'Removing',out_cnt,'particles'
     &,'. Total remaining:',np
      endif

      deallocate (out_pt)
      deallocate (xp_pt,yp_pt,zp_pt)
      deallocate (uop_pt,vop_pt,wop_pt)
      deallocate (dp_pt,rho_pt)
      deallocate (Fu,Fv,Fw)

      allocate (xp_pt(np),yp_pt(np),zp_pt(np))
      allocate (uop_pt(np),vop_pt(np),wop_pt(np))
      allocate (dp_pt(np),rho_pt(np))
      allocate (Fu(np),Fv(np),Fw(np))

      do l=1,np
                  xp_pt(l) = xpold(l)
                  yp_pt(l) = ypold(l)
                  zp_pt(l) = zpold(l)
                  uop_pt(l)= uopold(l)
                  vop_pt(l)= vopold(l)
                  wop_pt(l)= wopold(l)
                  dp_pt(l) = dp_old(l)
                  rho_pt(l) = rhop_old(l)
      enddo

      deallocate (xpold,ypold,zpold,uopold,vopold,wopold)
      deallocate (rhop_old,dp_old)

      allocate(xpold(np),ypold(np),zpold(np))
      allocate(uopold(np),vopold(np),wopold(np)) 
      allocate(dp_old(np),rhop_old(np))

!------------------------------------ 
        xpold=xp_pt 
        ypold=yp_pt
        zpold=zp_pt

        uopold=uop_pt
        vopold=vop_pt
        wopold=wop_pt

        dp_old=dp_pt
        rhop_old=rho_pt
!---------------------------------------------------

!       wtime_refresh = MPI_WTIME ( ) - wtime_refresh
!       write(203,*)'reallocate parts: ',wtime_refresh      

      return
      end subroutine


!##########################################################################
        subroutine release_pt
!##########################################################################
      use vars
      use mpi
      use multidata
      use vars_pt
      implicit none

      integer l,np_old,ptnr,tsnr,nfrac,f,m,frac1,frac_end
      double precision random_number_normal,random_number_uniform
      double precision :: xp,yp,zp,uop,vop,wop,Dp,sigma,rho_p
      double precision :: Wx,Wy,Wz,sigma_rho
      double precision :: mindis,dist,distance
      logical :: random

!       wtime_release = MPI_WTIME ( )
!------------------------------------ 
        xpold=xp_pt 
        ypold=yp_pt
        zpold=zp_pt

        uopold=uop_pt
        vopold=vop_pt
        wopold=wop_pt

        dp_old=dp_pt
        rhop_old=rho_pt
!---------------------------------------------------
       

!1 calculating new np

      open(15,file='LPT.cin')
      read(15,*)
      read(15,*) 
      read(15,*) nfrac              !how many Lag fractions you want to calculate BF2023
      read(15,*) 
      read(15,*) 
      
      np_old=np
  
      do f=1,nfrac
            read(15,*) 
            read(15,*) tsnr
            read(15,*) ptnr            
            if ((tsnr.gt.0).and.(mod(itime,tsnr).eq.0)) then
                  np=np+ptnr   
                  write(202,*) ntime,'Releasing',ptnr,'new particles ', 
     &      'within fraction',f
            endif

            read(15,*)                    !dp
            read(15,*)                    !rhop
            read(15,*)                    !release volume
            read(15,*)random
            if (random) then 
                  read(15,*)
            else
                  do l=1,ptnr
                        read(15,*)
                  enddo
            endif
      enddo

      close(15)

      deallocate (xp_pt,yp_pt,zp_pt)
      deallocate (uop_pt,vop_pt,wop_pt)
      deallocate (dp_pt,rho_pt)
      deallocate (Fu,Fv,Fw)

      allocate (xp_pt(np),yp_pt(np),zp_pt(np))
      allocate (uop_pt(np),vop_pt(np),wop_pt(np))
      allocate (dp_pt(np),rho_pt(np))
      allocate (Fu(np),Fv(np),Fw(np))

!2 allocating and deallocating for 1,np_old

      do l=1,np_old
            xp_pt(l)=xpold(l)
            yp_pt(l)=ypold(l)
            zp_pt(l)=zpold(l)

            uop_pt(l)=uopold(l)
            vop_pt(l)=vopold(l)
            wop_pt(l)=wopold(l)    

            dp_pt(l)=dp_old(l)
            rho_pt(l)=rhop_old(l)
      end do
         deallocate (xpold,ypold,zpold)
         allocate(xpold(np),ypold(np),zpold(np))
         deallocate (uopold,vopold,wopold)
         allocate (uopold(np),vopold(np),wopold(np))
         deallocate (dp_old,rhop_old)
         allocate (dp_old(np),rhop_old(np))

!3 initiallising variables for np_old,np

      open(35,file='LPT.cin')       !reopen the file to read the details of every fraction
      read(35,*)   
      read(35,*)                    !PSIcell
      read(35,*)                    !nfrac
      read(35,*)                    !Lcol
      read(35,*)                    !Lcolwall

      frac1=np_old+1
      do f=1,nfrac
            read(35,*)                          !fraction header
            read(35,*) tsnr                     !if 0 then no continuous release, only released on init
            read(35,*) ptnr 
            frac_end=frac1+ptnr-1
            read(35,*) Dp,sigma
            read(35,*) rho_p,sigma_rho
            read(35,*) Wx,Wy,Wz
            read(35,*) random
            if (random) read(35,*)xp,yp,zp,uop,vop,wop

            do l=frac1,frac_end
                  if ((tsnr.gt.0).and.(mod(itime,tsnr).eq.0)) then

                  if (random) then                                      !location

                  mindis=1.1*Dp
                  dist=0

                  do while (dist.lt.mindis)                             !avoiding overlap
                  dist=mindis
                  xp_pt(l)=random_number_uniform(xp-0.5*Wx,xp+0.5*Wx)
                  yp_pt(l)=random_number_uniform(yp-0.5*Wy,yp+0.5*Wy)
                  zp_pt(l)=random_number_uniform(zp-0.5*Wz,zp+0.5*Wz)                      
                  do m=frac1,l
                  if (l.ne.m) then 
                        distance=sqrt((xp_pt(l)-xp_pt(m))**2+(yp_pt(l)
     &                  -yp_pt(m))**2+(zp_pt(l)-zp_pt(m))**2)
                  else
                        distance=1.d9
                  endif
                        dist=min(dist,distance)
                  enddo
                  enddo
                  uop_pt(l)=uop
                  vop_pt(l)=vop
                  wop_pt(l)=wop
                  else                                                  !not random, read from file
                        read(35,*)xp_pt(l),yp_pt(l),zp_pt(l),
     &                  uop_pt(l),vop_pt(l),wop_pt(l)
                  endif                                                 !random location
                  dp_pt(l)= random_number_normal(Dp,sigma)
                  rho_pt(l)= random_number_normal(rho_p,sigma_rho)

                  xpold(l)=xp_pt(l)
                  ypold(l)=yp_pt(l)
                  zpold(l)=zp_pt(l)
                  uopold(l)=uop_pt(l)
                  vopold(l)=vop_pt(l)
                  wopold(l)=wop_pt(l)
                  dp_old(l)=dp_pt(l)
                  Fu(l)=0 ; Fv(l)=0 ; Fw(l)=0
                  else                                                  !tsnr
                        if (.not.random) read(35,*)                           
                  endif                                                 !tsnr
            enddo                                                       !particles
                  
            if ((tsnr.gt.0).and.(mod(itime,tsnr).eq.0)) frac1=frac1+ptnr
            
      enddo                                                             !loop in fracs

      close(35)

      return
      end subroutine
!##########################################################################
        subroutine periodic_pt
!##########################################################################
        use vars
        use vars_pt
        use mpi
        use multidata
        implicit none

        integer l

      do l=1,np

!     Comprobar si permanece en el dominio
      if (xp_pt(l).lt.xst) then
            xp_pt(l)=xp_pt(l)+(xen-xst)                     !bubble comes back at the top (unlikely)
      elseif(xp_pt(l).gt.xen) then
            xp_pt(l)=xp_pt(l)-(xen-xst)                     !bubble comes back at the bottom (likely)
      endif
C       if (yp_pt(l).lt.yst) then
C             yp_pt(l)=yp_pt(l)+(yen-yst)                     !bubble comes back
C       elseif (yp_pt(l).gt.yen) then
C             yp_pt(l)=yp_pt(l)-(yen-yst)                     !bubble comes back
C       endif
C       if (zp_pt(l).lt.zst) then
C             zp_pt(l)=zp_pt(l)+(zen-zst)                     !bubble comes back
C       elseif (zp_pt(l).gt.zen) then
C             zp_pt(l)=zp_pt(l)-(zen-zst)                     !bubble comes back
C       end if

      enddo


      return
      end subroutine
