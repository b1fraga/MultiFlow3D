module multiflow3d_MPI_pt
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use mpi_f08, only: MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_INTEGER, &
                     MPI_BARRIER, MPI_BCAST, MPI_SCATTER
  use multidata, only: idom, jdom, kdom, dom_ad
  use multiflow3d_mpi, only : ierr, Myrank, nprocs
  use vars, only: g_dx, g_dy, g_dz, np, xen, xst, yen, yst, zen, zst
  use vars_pt, only: Lcol, np_loc, npg_loc, ptsinproc, ptsinproc_g, xp_pt, &
                     yp_pt, zp_pt, uop_pt, vop_pt, wop_pt, dp_pt, rho_pt, &
                     xp_loc, yp_loc, zp_loc, uop_loc, vop_loc, wop_loc, dp_loc, &
                     rhop_loc, id, xpg_loc, ypg_loc, zpg_loc, uopg_loc, &
                     vopg_loc, wopg_loc, dpg_loc, rhopg_loc

  implicit none
  private

  public :: MPI_pt

contains

  !======================================================================!
  !                 LAGRANGIAN PARTICLE TRACKING                         !
  !----------------------------------------------------------------------!
  !                           Bruño Fraga                                !
  !                           Boyang Chen                                !
  !                      Cardiff Uni 2013-2017                           !
  !                        Stanford Uni 2018                             !
  !                    Uni of Birmingham 2019-2024                       !
  !======================================================================!
  !######################################################################
  subroutine MPI_pt
  !     Distributes the particles among the processors
  !     Creates ghost particles for collision calculations
  !######################################################################
  integer :: l,o,ii,nxdom,nydom,nzdom,ntot,lp
  integer :: nxdom_iprev,nydom_iprev,nzdom_iprev
  integer :: nxdom_jprev,nydom_jprev,nzdom_jprev
  integer :: nxdom_kprev,nydom_kprev,nzdom_kprev
  integer :: nxdom_inext,nydom_inext,nzdom_inext
  integer :: nxdom_jnext,nydom_jnext,nzdom_jnext
  integer :: nxdom_knext,nydom_knext,nzdom_knext

  real(dp) :: lx,ly,lz,dx,dy,dz
  integer,allocatable,dimension(:)::  lpt_proc,lpt_block
  integer,allocatable,dimension(:)::  id_MPI,id_MPI_loc
  real(dp), allocatable, dimension(:):: X_MPI,Y_MPI,Z_MPI
  real(dp), allocatable, dimension(:):: X_MPI_loc,Y_MPI_loc
  real(dp), allocatable, dimension(:):: Z_MPI_loc
  real(dp), allocatable, dimension(:):: U_MPI,V_MPI,W_MPI
  real(dp), allocatable, dimension(:):: U_MPI_loc,V_MPI_loc
  real(dp), allocatable, dimension(:):: W_MPI_loc
  real(dp), allocatable, dimension(:):: dp_MPI,dp_MPI_loc
  real(dp), allocatable, dimension(:):: rhop_MPI
  real(dp), allocatable, dimension(:):: rhop_MPI_loc

  ! ---------------------------------------------------
  ! ghost particles
  ! ---------------------------------------------------
  integer :: ll,oo,iii
  integer,allocatable,dimension(:)::  idg_MPI,idg_MPI_loc
  real(dp), allocatable, dimension(:):: Xg_MPI,Yg_MPI
  real(dp), allocatable, dimension(:):: Xg_MPI_loc,Zg_MPI
  real(dp), allocatable, dimension(:):: Yg_MPI_loc
  real(dp), allocatable, dimension(:):: Zg_MPI_loc
  real(dp), allocatable, dimension(:):: Ug_MPI,Vg_MPI
  real(dp), allocatable, dimension(:):: Ug_MPI_loc,Wg_MPI
  real(dp), allocatable, dimension(:):: Vg_MPI_loc
  real(dp), allocatable, dimension(:):: Wg_MPI_loc
  real(dp), allocatable, dimension(:):: dpg_MPI
  real(dp), allocatable, dimension(:):: dpg_MPI_loc
  real(dp), allocatable, dimension(:):: rhopg_MPI
  real(dp), allocatable, dimension(:):: rhopg_MPI_loc

  integer,allocatable,dimension(:):: lpt_proc_iprev,lpt_block_iprev
  integer,allocatable,dimension(:):: lpt_proc_jprev,lpt_block_jprev
  integer,allocatable,dimension(:):: lpt_proc_kprev,lpt_block_kprev
  integer,allocatable,dimension(:):: lpt_proc_inext,lpt_block_inext
  integer,allocatable,dimension(:):: lpt_proc_jnext,lpt_block_jnext
  integer,allocatable,dimension(:):: lpt_proc_knext,lpt_block_knext


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

  if (Lcol) then
     !  ghost parts
     allocate(Xg_MPI(ntot),Yg_MPI(ntot),Zg_MPI(ntot))
     allocate(Ug_MPI(ntot),Vg_MPI(ntot),Wg_MPI(ntot))
     allocate(dpg_MPI(ntot),idg_MPI(ntot),rhopg_MPI(ntot))
     allocate(Xg_MPI_loc(np),Yg_MPI_loc(np),Zg_MPI_loc(np))
     allocate(Ug_MPI_loc(np),Vg_MPI_loc(np),Wg_MPI_loc(np))
     allocate(dpg_MPI_loc(np),rhopg_MPI_loc(np),idg_MPI_loc(np))
     !  neighbours
     allocate(lpt_block_iprev(np),lpt_proc_iprev(np))
     allocate(lpt_block_jprev(np),lpt_proc_jprev(np))
     allocate(lpt_block_kprev(np),lpt_proc_kprev(np))
     allocate(lpt_block_inext(np),lpt_proc_inext(np))
     allocate(lpt_block_jnext(np),lpt_proc_jnext(np))
     allocate(lpt_block_knext(np),lpt_proc_knext(np))

  end if  !Lcol

  if (Myrank==0) then

     X_MPI = -1.0d12
     if (Lcol) Xg_MPI = -1.0d12

     lx=xen-xst
     ly=yen-yst
     lz=zen-zst

     dx = g_dx
     dy = g_dy
     dz = g_dz

     ptsinproc=0
     lpt_block=0
     lpt_proc=0

     do l=1,np

        nxdom=INT((xp_pt(l)-xst-1d-12)/(lx/idom))  !Domain where point belongs in i

        nydom=INT((yp_pt(l)-yst-1d-12)/(ly/jdom))  !Domain where point belongs in j

        nzdom=INT((zp_pt(l)-zst-1d-12)/(lz/kdom))  !Domain where point belongs in k


        lpt_block(l)=idom*jdom*(nzdom)+idom*(nydom) &
                     +nxdom                        !Domain where point belongs (0 to num_dom)

        lpt_proc(l)=dom_ad(lpt_block(l))+1         !Processor where point belongs (1 to num_procs+1)

        ptsinproc(lpt_proc(l))=ptsinproc(lpt_proc(l))+1
        !array which tells how many points are in each processor


     end do

     if (Lcol) then
        ptsinproc_g=0
        lpt_block_iprev=0
        lpt_proc_iprev=0
        lpt_block_jprev=0
        lpt_proc_jprev=0
        lpt_block_kprev=0
        lpt_proc_kprev=0
        lpt_block_inext=0
        lpt_proc_inext=0
        lpt_block_jnext=0
        lpt_proc_jnext=0
        lpt_block_knext=0
        lpt_proc_knext=0

        !..........................................................................
        !=== Ghost particles     ===> location of neighbours
        !..........................................................................
        do lp=1,np
           !=== Previous Neighbour  <===
           nxdom_iprev = INT((xp_pt(lp)-dx-xst-1d-12)/(lx/idom))
           nydom_iprev = INT((yp_pt(lp)-yst-1d-12)/(ly/jdom))
           nzdom_iprev = INT((zp_pt(lp)-zst-1d-12)/(lz/kdom))
           lpt_block_iprev(lp) = idom*jdom*(nzdom_iprev)+idom &
                                 *(nydom_iprev)+nxdom_iprev
           if (lpt_block_iprev(lp)/=lpt_block(lp)) then
              lpt_proc_iprev(lp)=dom_ad(lpt_block_iprev(lp))+1
              ptsinproc_g(lpt_proc_iprev(lp))= &
                                              ptsinproc_g(lpt_proc_iprev(lp))+1
           end if

           ! ---
           nxdom_jprev = INT((xp_pt(lp)-xst-1d-12)/(lx/idom))
           nydom_jprev = INT((yp_pt(lp)-dy-yst-1d-12)/(ly/jdom))
           nzdom_jprev = INT((zp_pt(lp)-zst-1d-12)/(lz/kdom))
           lpt_block_jprev(lp) = idom*jdom*(nzdom_jprev)+idom*(nydom_jprev) &
                                    +nxdom_jprev

           if (lpt_block_jprev(lp)/=lpt_block(lp)) then
              lpt_proc_jprev(lp)=dom_ad(lpt_block_jprev(lp))+1
              ptsinproc_g(lpt_proc_jprev(lp))= &
                                              ptsinproc_g(lpt_proc_jprev(lp))+1
           end if

           ! ---
           nxdom_kprev = INT((xp_pt(lp)-xst-1d-12)/(lx/idom))
           nydom_kprev = INT((yp_pt(lp)-yst-1d-12)/(ly/jdom))
           nzdom_kprev = INT((zp_pt(lp)-dz-zst-1d-12)/(lz/kdom))
           lpt_block_kprev(lp) = idom*jdom*(nzdom_kprev) &
                                 +idom*(nydom_kprev)+nxdom_kprev

           if (lpt_block_kprev(lp)/=lpt_block(lp)) then
              lpt_proc_kprev(lp)=dom_ad(lpt_block_kprev(lp))+1
              ptsinproc_g(lpt_proc_kprev(lp))= &
                                              ptsinproc_g(lpt_proc_kprev(lp))+1
           end if

           !=== Next Neighbor  ===>
           if (abs(xp_pt(lp)-xen)>dx) then
              nxdom_inext = INT((xp_pt(lp)+dx-xst-1d-12)/(lx/idom))
              nydom_inext = INT((yp_pt(lp)-yst-1d-12)/(ly/jdom))
              nzdom_inext = INT((zp_pt(lp)-zst-1d-12)/(lz/kdom))
              lpt_block_inext(lp) = idom*jdom*(nzdom_inext)+idom*(nydom_inext) &
                                    +nxdom_inext
              if (lpt_block_inext(lp)/=lpt_block(lp)) then
                 lpt_proc_inext(lp)=dom_ad(lpt_block_inext(lp))+1
                 ptsinproc_g(lpt_proc_inext(lp))= &
                                                 ptsinproc_g(lpt_proc_inext(lp))+1
              end if
           end if

           ! -+-
           if (abs(yp_pt(lp)-yen)>dy) then
              nxdom_jnext = INT((xp_pt(lp)-xst-1d-12)/(lx/idom))
              nydom_jnext = INT((yp_pt(lp)+dy-yst-1d-12)/(ly/jdom))
              nzdom_jnext = INT((zp_pt(lp)-zst-1d-12)/(lz/kdom))
              lpt_block_jnext(lp) = idom*jdom*(nzdom_jnext)+idom*(nydom_jnext) &
                                    +nxdom_jnext
              if (lpt_block_jnext(lp)/=lpt_block(lp)) then
                 lpt_proc_jnext(lp)=dom_ad(lpt_block_jnext(lp))+1
                 ptsinproc_g(lpt_proc_jnext(lp))= &
                                                 ptsinproc_g(lpt_proc_jnext(lp))+1
              end if
           end if

           ! --+
           if (abs(zp_pt(lp)-zen)>dz) then
              nxdom_knext = INT((xp_pt(lp)-xst-1d-12)/(lx/idom))
              nydom_knext = INT((yp_pt(lp)-yst-1d-12)/(ly/jdom))
              nzdom_knext = INT((zp_pt(lp)+dz-zst-1d-12)/(lz/kdom))
              lpt_block_knext(lp) = idom*jdom*(nzdom_knext)+idom*(nydom_knext) &
                                   +nxdom_knext
              if (lpt_block_knext(lp)/=lpt_block(lp)) then
                 lpt_proc_knext(lp)=dom_ad(lpt_block_knext(lp))+1
                 ptsinproc_g(lpt_proc_knext(lp))= &
                                                 ptsinproc_g(lpt_proc_knext(lp))+1
              end if
           end if
        end do
     end if  !Lcol

     ! --------------------------------
     ! real particles
     ! --------------------------------
     do o=1,nprocs
        ii=1
        do l=1,np
           if(lpt_proc(l)==(o)) then
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
           end if
        end do
     end do

     ! --------------------------------
     ! ghost particles
     ! --------------------------------
     if (Lcol) then
        do oo=1,nprocs
           iii=1
           do ll=1,np

              if (&                             ! ghost particles
                   (lpt_proc_iprev(ll)==(oo)).or.(lpt_proc_jprev(ll)==(oo)).or. &
                   (lpt_proc_kprev(ll)==(oo)).or.(lpt_proc_inext(ll)==(oo)).or. &
                   (lpt_proc_jnext(ll)==(oo)).or. &
                   (lpt_proc_knext(ll)==(oo))) then

                 Xg_MPI(iii+np*(oo-1))=xp_pt(ll)                      !superarrays
                 Yg_MPI(iii+np*(oo-1))=yp_pt(ll)
                 Zg_MPI(iii+np*(oo-1))=zp_pt(ll)

                 Ug_MPI(iii+np*(oo-1))=uop_pt(ll)
                 Vg_MPI(iii+np*(oo-1))=vop_pt(ll)
                 Wg_MPI(iii+np*(oo-1))=wop_pt(ll)

                 rhopg_MPI(iii+np*(oo-1))=rho_pt(ll)
                 dpg_MPI(iii+np*(oo-1))=dp_pt(ll)
                 idg_MPI(iii+np*(oo-1))=lpt_block(ll)
                 iii=iii+1

              end if
           end do
        end do

     end if  !Lcol

  end if  !master proc

  call MPI_BARRIER (MPI_COMM_WORLD,ierr)

  !  number of real particles per block
  call MPI_BCAST(ptsinproc,nprocs,MPI_INTEGER,0, &
                 MPI_COMM_WORLD,ierr)

  call MPI_SCATTER(ptsinproc,1,MPI_INTEGER,np_loc,1, &
                   MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  if (Lcol) then

     !  number of ghost particles per block
     call MPI_BCAST(ptsinproc_g,nprocs,MPI_INTEGER,0, &
                    MPI_COMM_WORLD,ierr)

     call MPI_SCATTER(ptsinproc_g,1,MPI_INTEGER,npg_loc,1, &
                      MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


  end if  !Lcol

  ! --------------------------------
  ! real particles
  ! --------------------------------

  call MPI_SCATTER(X_MPI,np,MPI_DOUBLE_PRECISION,X_MPI_loc,&       !location
                   np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_SCATTER(Y_MPI,np,MPI_DOUBLE_PRECISION,Y_MPI_loc, &
                   np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_SCATTER(Z_MPI,np,MPI_DOUBLE_PRECISION,Z_MPI_loc, &
                   np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  call MPI_SCATTER(U_MPI,np,MPI_DOUBLE_PRECISION,U_MPI_loc,&       !velocities (of prev tstep)
                   np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_SCATTER(V_MPI,np,MPI_DOUBLE_PRECISION,V_MPI_loc, &
                   np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_SCATTER(W_MPI,np,MPI_DOUBLE_PRECISION,W_MPI_loc, &
                   np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  call MPI_SCATTER(dp_MPI,np,MPI_DOUBLE_PRECISION,dp_MPI_loc,&     !diameter
                   np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  call MPI_SCATTER(rhop_MPI,np,MPI_DOUBLE_PRECISION,rhop_MPI_loc,&  !density
                   np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  call MPI_SCATTER(id_MPI,np,MPI_INTEGER,id_MPI_loc,&              !block to which they belong
                   np,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  ! --------------------------------
  ! ghost particles
  ! --------------------------------
  if (Lcol) then

     call MPI_SCATTER(Xg_MPI,np,MPI_DOUBLE_PRECISION,Xg_MPI_loc,&       !location
                      np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     call MPI_SCATTER(Yg_MPI,np,MPI_DOUBLE_PRECISION,Yg_MPI_loc, &
                      np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     call MPI_SCATTER(Zg_MPI,np,MPI_DOUBLE_PRECISION,Zg_MPI_loc, &
                      np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


     call MPI_SCATTER(Ug_MPI,np,MPI_DOUBLE_PRECISION,Ug_MPI_loc,&       !velocities (of prev tstep)
                      np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     call MPI_SCATTER(Vg_MPI,np,MPI_DOUBLE_PRECISION,Vg_MPI_loc, &
                      np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
     call MPI_SCATTER(Wg_MPI,np,MPI_DOUBLE_PRECISION,Wg_MPI_loc, &
                      np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

     call MPI_SCATTER(dpg_MPI,np,MPI_DOUBLE_PRECISION,dpg_MPI_loc,&     !diameter
                      np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

     call MPI_SCATTER(rhopg_MPI,np,MPI_DOUBLE_PRECISION,&              !density
                      rhopg_MPI_loc,np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

     call MPI_SCATTER(idg_MPI,np,MPI_INTEGER,idg_MPI_loc,&              !block to which they belong
                      np,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  end if  !Lcol

  ! --------------------------------
  ! real particles
  ! --------------------------------

  if (np_loc>0) then

     allocate (xp_loc(np_loc),yp_loc(np_loc),zp_loc(np_loc))
     allocate (uop_loc(np_loc),vop_loc(np_loc),wop_loc(np_loc))
     allocate (dp_loc(np_loc),rhop_loc(np_loc))
     allocate (id(np_loc))

     ii=0
     do l=1,np
        if (X_MPI_loc(l)/=-1.0d12) then
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
        end if
     end do

     if (ii/=np_loc) then
        write(6,*)"MPI ERROR in proc:",myrank
        write(6,*)"np_loc=",np_loc,"=/=",ii
        stop
     end if

  end if

  ! --------------------------------
  ! ghost particles
  ! --------------------------------
  if (Lcol) then

     if (npg_loc>0) then

        allocate (xpg_loc(npg_loc),ypg_loc(npg_loc),zpg_loc(npg_loc))
        allocate (uopg_loc(npg_loc),vopg_loc(npg_loc),wopg_loc(npg_loc))
        allocate (dpg_loc(npg_loc),rhopg_loc(npg_loc))

        iii=0
        do ll=1,np
           if (Xg_MPI_loc(ll)/=-1.0d12) then
              iii=iii+1
              xpg_loc(ll)=Xg_MPI_loc(ll)
              ypg_loc(ll)=Yg_MPI_loc(ll)
              zpg_loc(ll)=Zg_MPI_loc(ll)

              uopg_loc(ll)=Ug_MPI_loc(ll)
              vopg_loc(ll)=Vg_MPI_loc(ll)
              wopg_loc(ll)=Wg_MPI_loc(ll)

              dpg_loc(ll)=dpg_MPI_loc(ll)
              rhopg_loc(ll)=rhopg_MPI_loc(ll)

           end if
        end do

        if (iii/=npg_loc) then
           write(6,*)"MPI ERROR in proc (ghost particles):",myrank
           write(6,*)"npg_loc=",npg_loc,"=/=",iii
           stop
        end if

     end if  !npg_loc

  end if  !Lcol

  deallocate(lpt_proc,lpt_block)
  deallocate(X_MPI,Y_MPI,Z_MPI)
  deallocate(U_MPI,V_MPI,W_MPI)
  deallocate(X_MPI_loc,Y_MPI_loc,Z_MPI_loc)
  deallocate(U_MPI_loc,V_MPI_loc,W_MPI_loc)
  deallocate(dp_MPI,dp_MPI_loc)
  deallocate(rhop_MPI,rhop_MPI_loc)
  deallocate(id_MPI,id_MPI_loc)

  if (Lcol) then
     deallocate(Xg_MPI,Yg_MPI,Zg_MPI)
     deallocate(Ug_MPI,Vg_MPI,Wg_MPI)
     deallocate(Xg_MPI_loc,Yg_MPI_loc,Zg_MPI_loc)
     deallocate(Ug_MPI_loc,Vg_MPI_loc,Wg_MPI_loc)
     deallocate(dpg_MPI,dpg_MPI_loc)
     deallocate(idg_MPI,idg_MPI_loc)
     deallocate(rhopg_MPI,rhopg_MPI_Loc)
     ! neighbors
     deallocate(lpt_proc_iprev,lpt_block_iprev)
     deallocate(lpt_proc_jprev,lpt_block_jprev)
     deallocate(lpt_proc_kprev,lpt_block_kprev)
     deallocate(lpt_proc_inext,lpt_block_inext)
     deallocate(lpt_proc_jnext,lpt_block_jnext)
     deallocate(lpt_proc_knext,lpt_block_knext)
  end if  !Lcol


  end subroutine MPI_pt

end module multiflow3d_MPI_pt
