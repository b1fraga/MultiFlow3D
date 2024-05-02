!======================================================================!
!                 LAGRANGIAN PARTICLE TRACKING                         !
!----------------------------------------------------------------------!
!                           Bruño Fraga                                ! 
!                           Boyang Chen                                !     
!                      Cardiff Uni 2013-2017                           !
!                        Stanford Uni 2018                             !
!                    Uni of Birmingham 2019-2023                       !
!======================================================================!
!######################################################################
      SUBROUTINE MPI_pt
!     Distributes the particles among the processors
!     Creates ghost particles for collision calculations
!######################################################################
      use vars
      use mpi
      use multidata
      use vars_pt

      implicit none

      integer :: l,o,ii,nxdom,nydom,nzdom,ntot,lp
C       integer :: nxdom_iprev,nydom_iprev,nzdom_iprev
C       integer :: nxdom_jprev,nydom_jprev,nzdom_jprev
C       integer :: nxdom_kprev,nydom_kprev,nzdom_kprev
C       integer :: nxdom_inext,nydom_inext,nzdom_inext
C       integer :: nxdom_jnext,nydom_jnext,nzdom_jnext
C       integer :: nxdom_knext,nydom_knext,nzdom_knext

C       integer :: nxdom_corprev1,nydom_corprev1,nzdom_corprev1
C       integer :: nxdom_corprev2,nydom_corprev2,nzdom_corprev2
C       integer :: nxdom_corprev3,nydom_corprev3,nzdom_corprev3
C       integer :: nxdom_corprev4,nydom_corprev4,nzdom_corprev4
C       integer :: nxdom_cornext1,nydom_cornext1,nzdom_cornext1
C       integer :: nxdom_cornext2,nydom_cornext2,nzdom_cornext2
C       integer :: nxdom_cornext3,nydom_cornext3,nzdom_cornext3
C       integer :: nxdom_cornext4,nydom_cornext4,nzdom_cornext4

C       integer :: nxdom_edgprev1,nydom_edgprev1,nzdom_edgprev1
C       integer :: nxdom_edgprev2,nydom_edgprev2,nzdom_edgprev2
C       integer :: nxdom_edgprev3,nydom_edgprev3,nzdom_edgprev3
C       integer :: nxdom_edgprev4,nydom_edgprev4,nzdom_edgprev4
C       integer :: nxdom_edgprev5,nydom_edgprev5,nzdom_edgprev5
C       integer :: nxdom_edgprev6,nydom_edgprev6,nzdom_edgprev6
C       integer :: nxdom_edgnext1,nydom_edgnext1,nzdom_edgnext1
C       integer :: nxdom_edgnext2,nydom_edgnext2,nzdom_edgnext2
C       integer :: nxdom_edgnext3,nydom_edgnext3,nzdom_edgnext3
C       integer :: nxdom_edgnext4,nydom_edgnext4,nzdom_edgnext4
C       integer :: nxdom_edgnext5,nydom_edgnext5,nzdom_edgnext5
C       integer :: nxdom_edgnext6,nydom_edgnext6,nzdom_edgnext6      
      real :: lx,ly,lz,dx,dy,dz
      integer,allocatable,dimension(:)::  lpt_proc,lpt_block
      integer,allocatable,dimension(:)::  id_MPI,id_MPI_loc
      double precision, allocatable, dimension(:):: X_MPI,Y_MPI,Z_MPI
      double precision, allocatable, dimension(:):: X_MPI_loc,Y_MPI_loc
      double precision, allocatable, dimension(:):: Z_MPI_loc
      double precision, allocatable, dimension(:):: U_MPI,V_MPI,W_MPI
      double precision, allocatable, dimension(:):: U_MPI_loc,V_MPI_loc
      double precision, allocatable, dimension(:):: W_MPI_loc
      double precision, allocatable, dimension(:):: dp_MPI,dp_MPI_loc
      double precision, allocatable, dimension(:):: rhop_MPI
      double precision, allocatable, dimension(:):: rhop_MPI_loc

! ---------------------------------------------------
! ghost particles 
! ---------------------------------------------------
C       integer :: ll,oo,iii,ib
C       integer,allocatable,dimension(:)::  idg_MPI,idg_MPI_loc
C       double precision, allocatable, dimension(:):: Xg_MPI,Yg_MPI,Zg_MPI
C       double precision, allocatable, dimension(:):: Xg_MPI_loc
C       double precision, allocatable, dimension(:):: Yg_MPI_loc
C       double precision, allocatable, dimension(:):: Zg_MPI_loc
C       double precision, allocatable, dimension(:):: Ug_MPI,Vg_MPI,Wg_MPI
C       double precision, allocatable, dimension(:):: Ug_MPI_loc
C       double precision, allocatable, dimension(:):: Vg_MPI_loc
C       double precision, allocatable, dimension(:):: Wg_MPI_loc
C       double precision, allocatable, dimension(:):: dpg_MPI,dpg_MPI_loc
C       double precision, allocatable, dimension(:):: rhopg_MPI
C       double precision, allocatable, dimension(:):: rhopg_MPI_loc

C       integer,allocatable,dimension(:):: lpt_proc_iprev,lpt_block_iprev
C       integer,allocatable,dimension(:):: lpt_proc_jprev,lpt_block_jprev 
C       integer,allocatable,dimension(:):: lpt_proc_kprev,lpt_block_kprev
C       integer,allocatable,dimension(:):: lpt_proc_inext,lpt_block_inext
C       integer,allocatable,dimension(:):: lpt_proc_jnext,lpt_block_jnext
C       integer,allocatable,dimension(:):: lpt_proc_knext,lpt_block_knext


      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
      call MPI_BCAST(np,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      ntot=np*nprocs

      allocate(X_MPI(ntot),Y_MPI(ntot),Z_MPI(ntot))
      allocate(U_MPI(ntot),V_MPI(ntot),W_MPI(ntot))
      allocate(dp_MPI(ntot),rhop_MPI(ntot),id_MPI(ntot))
      allocate(X_MPI_loc(np),Y_MPI_loc(np),Z_MPI_loc(np))
      allocate(U_MPI_loc(np),V_MPI_loc(np),W_MPI_loc(np))
      allocate(dp_MPI_loc(np),rhop_MPI_loc(np),id_MPI_loc(np))
      allocate(lpt_block(np),lpt_proc(np))

C       allocate(Xg_MPI(ntot),Yg_MPI(ntot),Zg_MPI(ntot))
C       allocate(Ug_MPI(ntot),Vg_MPI(ntot),Wg_MPI(ntot))
C       allocate(dpg_MPI(ntot),idg_MPI(ntot))
C       allocate(Xg_MPI_loc(np),Yg_MPI_loc(np),Zg_MPI_loc(np))
C       allocate(Ug_MPI_loc(np),Vg_MPI_loc(np),Wg_MPI_loc(np))
C       allocate(dpg_MPI_loc(np),idg_MPI_loc(np))
C ! neighbors
C       allocate(lpt_block_iprev(np),lpt_proc_iprev(np))
C       allocate(lpt_block_jprev(np),lpt_proc_jprev(np))
C       allocate(lpt_block_kprev(np),lpt_proc_kprev(np))
C       allocate(lpt_block_inext(np),lpt_proc_inext(np))
C       allocate(lpt_block_jnext(np),lpt_proc_jnext(np))
C       allocate(lpt_block_knext(np),lpt_proc_knext(np))


      IF (Myrank.eq.0) THEN

      X_MPI = -1.d12
C      Xg_MPI = -1.d12


      lx=xen-xst
      ly=yen-yst
      lz=zen-zst

C       dx = g_dx
C       dy = g_dy
C       dz = g_dz

      ptsinproc=0  ; lpt_block=0  ; lpt_proc=0 !; ptsinproc_g=0

            do l=1,np

                  nxdom=INT((xp_pt(l)-xst-1d-12)/(lx/idom))                   !Domain where point belongs in i

                  nydom=INT((yp_pt(l)-yst-1d-12)/(ly/jdom))                   !Domain where point belongs in j

                  nzdom=INT((zp_pt(l)-zst-1d-12)/(lz/kdom))                   !Domain where point belongs in k


                  lpt_block(l)=idom*jdom*(nzdom)+idom*(nydom)
     &                              +nxdom                              !Domain where point belongs (0 to num_dom)

                  lpt_proc(l)=dom_ad(lpt_block(l))+1                    !Processor where point belongs (1 to num_procs+1)

!                 ptsinblk(lpt_block(l)+1)=ptsinblk(lpt_block(l)+1)+1   !array which tells how many points are in each domain

                  ptsinproc(lpt_proc(l))=ptsinproc(lpt_proc(l))+1       !array which tells how many points are in each processor


!           write(6,*)lpt_block(l)
!                 write(6,*)xp_pt(l)
!                 write(6,*)(lx/idom)*(nxdom)
!           write(6,*)'======================================='

            enddo

C       lpt_block_iprev=0 ; lpt_proc_iprev=0 
C       lpt_block_jprev=0 ; lpt_proc_jprev=0
C       lpt_block_kprev=0 ; lpt_proc_kprev=0
C       lpt_block_inext=0 ; lpt_proc_inext=0
C       lpt_block_jnext=0 ; lpt_proc_jnext=0
C       lpt_block_knext=0 ; lpt_proc_knext=0

C       do lp=1,np
C !..........................................................................
C !=== Previous Neighbor  ===> 
C !..........................................................................
C !---
C        nxdom_iprev = INT((xp_pt(lp)-dx-xst-1d-12)/(lx/idom)) 
C        nydom_iprev = INT((yp_pt(lp)-yst-1d-12)/(ly/jdom))    
C        nzdom_iprev = INT((zp_pt(lp)-zst-1d-12)/(lz/kdom))    
C        lpt_block_iprev(lp) = idom*jdom*(nzdom_iprev)+idom
C      &                            *(nydom_iprev)+nxdom_iprev                                          
C        if (lpt_block_iprev(lp).ne.lpt_block(lp)) then
C        lpt_proc_iprev(lp)=dom_ad(lpt_block_iprev(lp))+1
C        ptsinproc_g(lpt_proc_iprev(lp))=
C      & ptsinproc_g(lpt_proc_iprev(lp))+1 
C        endif
C ! ---
C        nxdom_jprev = INT((xp_pt(lp)-xst-1d-12)/(lx/idom))    
C        nydom_jprev = INT((yp_pt(lp)-dy-yst-1d-12)/(ly/jdom)) 
C        nzdom_jprev = INT((zp_pt(lp)-zst-1d-12)/(lz/kdom))    
C        lpt_block_jprev(lp) = idom*jdom*(nzdom_jprev)+idom*(nydom_jprev)
C      &                              +nxdom_jprev                                    
C        if (lpt_block_jprev(lp).ne.lpt_block(lp)) then
C        lpt_proc_jprev(lp)=dom_ad(lpt_block_jprev(lp))+1
C        ptsinproc_g(lpt_proc_jprev(lp))=
C      & ptsinproc_g(lpt_proc_jprev(lp))+1 
C        endif
C ! ---
C        nxdom_kprev = INT((xp_pt(lp)-xst-1d-12)/(lx/idom))    
C        nydom_kprev = INT((yp_pt(lp)-yst-1d-12)/(ly/jdom))    
C        nzdom_kprev = INT((zp_pt(lp)-dz-zst-1d-12)/(lz/kdom))
C        lpt_block_kprev(lp) = idom*jdom*(nzdom_kprev)
C      &             +idom*(nydom_kprev)+nxdom_kprev                                    
C        if (lpt_block_kprev(lp).ne.lpt_block(lp)) then
C        lpt_proc_kprev(lp)=dom_ad(lpt_block_kprev(lp))+1
C        ptsinproc_g(lpt_proc_kprev(lp))=     
C      & ptsinproc_g(lpt_proc_kprev(lp))+1        
C        endif
C c           endif
C !..........................................................................
C !=== Next Neighbor  ===> 
C !..........................................................................
C ! +--
C        if (abs(xp_pt(lp)-xen).gt.dx) then 
C        nxdom_inext = INT((xp_pt(lp)+dx-xst-1d-12)/(lx/idom)) 
C        nydom_inext = INT((yp_pt(lp)-yst-1d-12)/(ly/jdom))    
C        nzdom_inext = INT((zp_pt(lp)-zst-1d-12)/(lz/kdom))
C        lpt_block_inext(lp) = idom*jdom*(nzdom_inext)+idom*(nydom_inext)
C      &                              +nxdom_inext                                    
C        if (lpt_block_inext(lp).ne.lpt_block(lp)) then
C        lpt_proc_inext(lp)=dom_ad(lpt_block_inext(lp))+1
C        ptsinproc_g(lpt_proc_inext(lp))=
C      & ptsinproc_g(lpt_proc_inext(lp))+1 
C        endif 
C        endif 
C ! -+-
C        if (abs(yp_pt(lp)-yen).gt.dy) then 
C        nxdom_jnext = INT((xp_pt(lp)-xst-1d-12)/(lx/idom))    
C        nydom_jnext = INT((yp_pt(lp)+dy-yst-1d-12)/(ly/jdom)) 
C        nzdom_jnext = INT((zp_pt(lp)-zst-1d-12)/(lz/kdom))
C        lpt_block_jnext(lp) = idom*jdom*(nzdom_jnext)+idom*(nydom_jnext)
C      &                              +nxdom_jnext                                    
C        if (lpt_block_jnext(lp).ne.lpt_block(lp)) then
C        lpt_proc_jnext(lp)=dom_ad(lpt_block_jnext(lp))+1
C        ptsinproc_g(lpt_proc_jnext(lp))=
C      & ptsinproc_g(lpt_proc_jnext(lp))+1 
C        endif
C        endif 
C ! --+
C        if (abs(zp_pt(lp)-zen).gt.dz) then 
C        nxdom_knext = INT((xp_pt(lp)-xst-1d-12)/(lx/idom))    
C        nydom_knext = INT((yp_pt(lp)-yst-1d-12)/(ly/jdom))    
C        nzdom_knext = INT((zp_pt(lp)+dz-zst-1d-12)/(lz/kdom)) 
C        lpt_block_knext(lp) = idom*jdom*(nzdom_knext)+idom*(nydom_knext)
C      &                             +nxdom_knext                                    
C        if (lpt_block_knext(lp).ne.lpt_block(lp)) then
C        lpt_proc_knext(lp)=dom_ad(lpt_block_knext(lp))+1 
C        ptsinproc_g(lpt_proc_knext(lp))=
C      & ptsinproc_g(lpt_proc_knext(lp))+1
C        endif 
C        endif

C             enddo

! --------------------------------
! real particles
! -------------------------------- 
            do o=1,nprocs
!           write(6,*)'1',o,ptsinproc(o)
               ii=1
               DO l=1,np
                  if(lpt_proc(l).eq.(o)) then
                        X_MPI(ii+np*(o-1))=xp_pt(l)                     !superarrays
                        Y_MPI(ii+np*(o-1))=yp_pt(l)
                        Z_MPI(ii+np*(o-1))=zp_pt(l)

                        U_MPI(ii+np*(o-1))=uop_pt(l)
                        V_MPI(ii+np*(o-1))=vop_pt(l)
                        W_MPI(ii+np*(o-1))=wop_pt(l)

                        dp_MPI(ii+np*(o-1))=dp_pt(l)

                        rhop_MPI(ii+np*(o-1))=rho_pt(l)

                        id_MPI(ii+np*(o-1))=lpt_block(l)

                        ii=ii+1
                  endif
               ENDDO
            enddo

C ! --------------------------------
C ! ghost particles
C ! -------------------------------- 

C             do oo=1,nprocs
C !           write(6,*)'1',o,ptsinproc(o)
C                iii=1
C                DO ll=1,np
C c              write(myrank+2000,*) oo,ll,
C c     &           lpt_proc_edgprev1(ll),lpt_proc_edgprev2(ll),
C c     &           lpt_proc_edgprev3(ll),lpt_proc_edgprev3(ll),
C c     &           lpt_proc_edgprev5(ll),lpt_proc_edgprev6(ll),
C c     &           lpt_proc_edgnext1(ll),lpt_proc_edgnext2(ll),
C c     &           lpt_proc_edgnext3(ll),lpt_proc_edgnext4(ll),
C c     &           lpt_proc_edgnext5(ll),lpt_proc_edgnext6(ll)

C clpt_proc_corprev1(ll),
C c     &           lpt_proc_corprev2(ll),lpt_proc_corprev3(ll),
C c     &           lpt_proc_corprev4(ll),lpt_proc_cornext1(ll),
C c     &           lpt_proc_cornext2(ll),
C c     &             lpt_proc_cornext3(ll),lpt_proc_cornext4(ll)

C                   if (                                      ! ghost particles 
C      &  (lpt_proc_iprev(ll).eq.(oo)).or.(lpt_proc_jprev(ll).eq.(oo)).or.
C      &  (lpt_proc_kprev(ll).eq.(oo)).or.(lpt_proc_inext(ll).eq.(oo)).or.
C      &  (lpt_proc_jnext(ll).eq.(oo)).or.
C      &  (lpt_proc_knext(ll).eq.(oo))) then                    !.or.       
C c     &           (lpt_proc_corprev1(ll).eq.(oo)).or.
C c     &           (lpt_proc_corprev2(ll).eq.(oo)).or. 
C c     &           (lpt_proc_corprev3(ll).eq.(oo)).or.
C c     &           (lpt_proc_corprev4(ll).eq.(oo)).or.
C c     &           (lpt_proc_cornext1(ll).eq.(oo)).or.
C c     &           (lpt_proc_cornext2(ll).eq.(oo)).or.
C c     &           (lpt_proc_cornext3(ll).eq.(oo)).or.
C c     &           (lpt_proc_cornext4(ll).eq.(oo)).or.
C c     &           (lpt_proc_edgprev1(ll).eq.(oo)).or.
C c     &           (lpt_proc_edgprev2(ll).eq.(oo)).or.
C c     &           (lpt_proc_edgprev3(ll).eq.(oo)).or.
C c     &           (lpt_proc_edgprev4(ll).eq.(oo)).or.
C c     &           (lpt_proc_edgprev5(ll).eq.(oo)).or.
C c     &           (lpt_proc_edgprev6(ll).eq.(oo)).or.
C c     &           (lpt_proc_edgnext1(ll).eq.(oo)).or.
C c     &           (lpt_proc_edgnext2(ll).eq.(oo)).or.
C c     &           (lpt_proc_edgnext3(ll).eq.(oo)).or.
C c     &           (lpt_proc_edgnext4(ll).eq.(oo)).or.
C c     &           (lpt_proc_edgnext5(ll).eq.(oo)).or.
C c     &           (lpt_proc_edgnext6(ll).eq.(oo))) then 

C c                       write(myrank+1000,*) iii+npt*(oo-1),ll,zp_pt(ll)
C                         Xg_MPI(iii+np*(oo-1))=xp_pt(ll)                      !superarrays
C                         Yg_MPI(iii+np*(oo-1))=yp_pt(ll)
C                         Zg_MPI(iii+np*(oo-1))=zp_pt(ll)

C                         Ug_MPI(iii+np*(oo-1))=uop_pt(ll)
C                         Vg_MPI(iii+np*(oo-1))=vop_pt(ll)
C                         Wg_MPI(iii+np*(oo-1))=wop_pt(ll)

C                         dpg_MPI(iii+np*(oo-1))=dp_pt(ll)
C                         idg_MPI(iii+np*(oo-1))=lpt_block(ll)
C c                       write(myrank+1000,*) Xg_MPI(iii+npt*(oo-1)),xp_pt(ll)
C                         iii=iii+1

C                   endif
C                ENDDO
C             enddo

      ENDIF

      call MPI_BARRIER (MPI_COMM_WORLD,ierr)

      call MPI_BCAST(ptsinproc,nprocs,MPI_INTEGER,0,
     &                  MPI_COMM_WORLD,ierr)

      call MPI_SCATTER(ptsinproc,1,MPI_INTEGER,np_loc,1,
     &            MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!     if (np_loc.eq.0) then
!           np_loc=1
!     endif


! --------------------------------
! real particles
! -------------------------------- 
!      call MPI_BARRIER (MPI_COMM_WORLD,ierr)

        call MPI_SCATTER(X_MPI,np,MPI_DOUBLE_PRECISION,X_MPI_loc,       !location
     &            np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_SCATTER(Y_MPI,np,MPI_DOUBLE_PRECISION,Y_MPI_loc,
     &            np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_SCATTER(Z_MPI,np,MPI_DOUBLE_PRECISION,Z_MPI_loc,
     &            np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        call MPI_SCATTER(U_MPI,np,MPI_DOUBLE_PRECISION,U_MPI_loc,       !velocities (of prev tstep)
     &            np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_SCATTER(V_MPI,np,MPI_DOUBLE_PRECISION,V_MPI_loc,
     &            np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_SCATTER(W_MPI,np,MPI_DOUBLE_PRECISION,W_MPI_loc,
     &            np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        call MPI_SCATTER(dp_MPI,np,MPI_DOUBLE_PRECISION,dp_MPI_loc,     !diameter
     &            np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        call MPI_SCATTER(rhop_MPI,np,MPI_DOUBLE_PRECISION,rhop_MPI_loc, !density
     &            np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        call MPI_SCATTER(id_MPI,np,MPI_INTEGER,id_MPI_loc,              !block to which they belong
     &            np,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

C ! --------------------------------
C ! ghost particles
C ! -------------------------------- 

C         call MPI_SCATTER(Xg_MPI,np,MPI_DOUBLE_PRECISION,Xg_MPI_loc,       !location
C      &            np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
C         call MPI_SCATTER(Yg_MPI,np,MPI_DOUBLE_PRECISION,Yg_MPI_loc,
C      &            np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
C         call MPI_SCATTER(Zg_MPI,np,MPI_DOUBLE_PRECISION,Zg_MPI_loc,
C      &            np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

C         call MPI_SCATTER(Ug_MPI,np,MPI_DOUBLE_PRECISION,Ug_MPI_loc,       !velocities (of prev tstep)
C      &            np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
C         call MPI_SCATTER(Vg_MPI,np,MPI_DOUBLE_PRECISION,Vg_MPI_loc,
C      &            np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
C         call MPI_SCATTER(Wg_MPI,np,MPI_DOUBLE_PRECISION,Wg_MPI_loc,
C      &            np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

C         call MPI_SCATTER(dpg_MPI,np,MPI_DOUBLE_PRECISION,dpg_MPI_loc,     !diameter
C      &            np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

C         call MPI_SCATTER(rhopg_MPI,np,MPI_DOUBLE_PRECISION,              !density
C      &  rhopg_MPI_loc,np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

C         call MPI_SCATTER(idg_MPI,np,MPI_INTEGER,idg_MPI_loc,              !block to which they belong
C      &            np,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


! --------------------------------
! real particles
! -------------------------------- 

      if (np_loc.gt.0) then

      allocate (xp_loc(np_loc),yp_loc(np_loc),zp_loc(np_loc))
      allocate (uop_loc(np_loc),vop_loc(np_loc),wop_loc(np_loc))
      allocate (dp_loc(np_loc),rhop_loc(np_loc))
!     allocate (Fsu(np_loc),Fau(np_loc),Fdu(np_loc),Flu(np_loc))
!     allocate (Fsv(np_loc),Fav(np_loc),Fdv(np_loc),Flv(np_loc))
!     allocate (Faw(np_loc),Fdw(np_loc),Flw(np_loc))
!     allocate (Fgw(np_loc),Fsw(np_loc))
!     allocate (ipux(np_loc),jpvy(np_loc),kpwz(np_loc))
      allocate (id(np_loc))
!     allocate (dh_acumu(np_loc),dh_acumv(np_loc),dh_acumw(np_loc))


            ii=0
            do l=1,np
                  if (X_MPI_loc(l).ne.-1.d12) then    
                        ii=ii+1           
                        xp_loc(l)=X_MPI_loc(l)
                        yp_loc(l)=Y_MPI_loc(l)
                        zp_loc(l)=Z_MPI_loc(l)

                        uop_loc(l)=U_MPI_loc(l)
                        vop_loc(l)=V_MPI_loc(l)
                        wop_loc(l)=W_MPI_loc(l)

                        dp_loc(l)=dp_MPI_loc(l)
                        rhop_loc(l)=rhop_MPI_loc(l)

                        id(l)=id_MPI_loc(l)
                  endif
            enddo
            if (ii.ne.np_loc) then
                  write(6,*)'MPI ERROR in proc:',myrank
                  write(6,*)'np_loc=',np_loc,'=/=',ii
                  stop
            endif

      endif

C ! --------------------------------
C ! ghost particles
C ! -------------------------------- 
C       if (npg_loc.gt.0) then

C       allocate (xpg_loc(npg_loc),ypg_loc(npg_loc),zpg_loc(npg_loc))
C       allocate (uopg_loc(npg_loc),vopg_loc(npg_loc),wopg_loc(npg_loc))
C       allocate (dpg_loc(npg_loc))
C       allocate (rhopg_loc(np_loc))

C             iii=0
C             do ll=1,np
C c                 write(myrank+2000,*) ll,Xg_MPI_loc(ll)
C                   if (Xg_MPI_loc(ll).ne.-1.d12) then  
C                         iii=iii+1         
C                         xpg_loc(ll)=Xg_MPI_loc(ll)
C                         ypg_loc(ll)=Yg_MPI_loc(ll)
C                         zpg_loc(ll)=Zg_MPI_loc(ll)

C                         uopg_loc(ll)=Ug_MPI_loc(ll)
C                         vopg_loc(ll)=Vg_MPI_loc(ll)
C                         wopg_loc(ll)=Wg_MPI_loc(ll)

C                         dpg_loc(ll)=dpg_MPI_loc(ll)
C                         rhopg_loc(l)=rhopg_MPI_loc(l)
C !                       write(myrank+1000,*) xp_loc(l), dp_loc(l), rho_loc(l)
C !                       idg(ll)=idg_MPI_loc(ll)
C                   endif
C             enddo
C             if (iii.ne.npg_loc) then
C               write(6,*)'MPI ERROR in proc (ghost particles):',myrank
C               write(6,*)'npg_loc=',npg_loc,'=/=',iii
C               stop
C             endif

C       endif

C             deallocate(lpt_proc,lpt_block)
C             deallocate(X_MPI,Y_MPI,Z_MPI)
C             deallocate(U_MPI,V_MPI,W_MPI)
C             deallocate(X_MPI_loc,Y_MPI_loc,Z_MPI_loc)
C             deallocate(U_MPI_loc,V_MPI_loc,W_MPI_loc)
C             deallocate(dp_MPI,dp_MPI_loc)
C             deallocate(rhop_MPI,rhop_MPI_loc)
C             deallocate(id_MPI,id_MPI_loc)

C             deallocate(Xg_MPI,Yg_MPI,Zg_MPI)
C             deallocate(Ug_MPI,Vg_MPI,Wg_MPI)
C             deallocate(Xg_MPI_loc,Yg_MPI_loc,Zg_MPI_loc)
C             deallocate(Ug_MPI_loc,Vg_MPI_loc,Wg_MPI_loc)
C             deallocate(dpg_MPI,dpg_MPI_loc)
C             deallocate(idg_MPI,idg_MPI_loc)
C             deallocate(rhopg_MPI,rhopg_MPI_Loc)
C ! neighbors 
C             deallocate(lpt_proc_iprev,lpt_block_iprev)
C             deallocate(lpt_proc_jprev,lpt_block_jprev)
C             deallocate(lpt_proc_kprev,lpt_block_kprev)
C             deallocate(lpt_proc_inext,lpt_block_inext)
C             deallocate(lpt_proc_jnext,lpt_block_jnext)
C             deallocate(lpt_proc_knext,lpt_block_knext)

      RETURN
      END
