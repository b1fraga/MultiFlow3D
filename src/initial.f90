!##########################################################################
      subroutine initial
!##########################################################################
          use, intrinsic :: iso_fortran_env, only: dp => real64
          use vars
          use multiflow3d_mpi
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          integer :: i,ib,tti,ttj,ttk
          integer :: glevel,gl,mgc_i,mgc_j,mgc_k,is,ie,js,je,ks,ke
          real(dp)    :: ndx,ndy,ndz,nwxend,nwyend,nwzend

          do ib=1,nbp
              tti=dom(ib)%ttc_i
              ttj=dom(ib)%ttc_j
              ttk=dom(ib)%ttc_k

              do i=1,26
                  dom(ib)%tg(i)=i*10**5+dom_id(ib)
              end do

              allocate(dom(ib)%x(tti),dom(ib)%y(ttj),dom(ib)%z(ttk))
              allocate(dom(ib)%xc(tti),dom(ib)%yc(ttj),dom(ib)%zc(ttk))

              allocate(dom(ib)%u(tti,ttj,ttk))
              allocate(dom(ib)%v(tti,ttj,ttk))
              allocate(dom(ib)%w(tti,ttj,ttk))
              allocate(dom(ib)%p(tti,ttj,ttk))
              allocate(dom(ib)%pp(tti,ttj,ttk))
              allocate(dom(ib)%sup(tti,ttj,ttk))
              allocate(dom(ib)%ustar(tti,ttj,ttk))
              allocate(dom(ib)%vstar(tti,ttj,ttk))
              allocate(dom(ib)%wstar(tti,ttj,ttk))

              allocate(dom(ib)%uo(tti,ttj,ttk),dom(ib)%uoo(tti,ttj,ttk))
              allocate(dom(ib)%vo(tti,ttj,ttk),dom(ib)%voo(tti,ttj,ttk))
              allocate(dom(ib)%wo(tti,ttj,ttk),dom(ib)%woo(tti,ttj,ttk))

              allocate(dom(ib)%ap(tti,ttj,ttk),dom(ib)%su(tti,ttj,ttk))
              allocate(dom(ib)%ae(tti,ttj,ttk),dom(ib)%aw(tti,ttj,ttk))
              allocate(dom(ib)%an(tti,ttj,ttk),dom(ib)%as(tti,ttj,ttk))
              allocate(dom(ib)%at(tti,ttj,ttk),dom(ib)%ab(tti,ttj,ttk))

              allocate(dom(ib)%um(tti,ttj,ttk),dom(ib)%vm(tti,ttj,ttk))
              allocate(dom(ib)%wm(tti,ttj,ttk),dom(ib)%pm(tti,ttj,ttk))
              allocate(dom(ib)%uum(tti,ttj,ttk),dom(ib)%vvm(tti,ttj,ttk))
              allocate(dom(ib)%wwm(tti,ttj,ttk),dom(ib)%uvm(tti,ttj,ttk))
              allocate(dom(ib)%uwm(tti,ttj,ttk),dom(ib)%vwm(tti,ttj,ttk))
              allocate(dom(ib)%ppm(tti,ttj,ttk))

              allocate(dom(ib)%ntav1(tti,ttj,ttk))
              allocate(dom(ib)%ntav2(tti,ttj,ttk))
              allocate(dom(ib)%facp1(tti,ttj,ttk))
              allocate(dom(ib)%facp2(tti,ttj,ttk))
              allocate(dom(ib)%facm1(tti,ttj,ttk))
              allocate(dom(ib)%facm2(tti,ttj,ttk))

              dom(ib)%ntav1=0
              dom(ib)%ntav2=0
              dom(ib)%facp1=0
              dom(ib)%facp2=0
              dom(ib)%facm1=1
              dom(ib)%facm2=1
              ntav1_count = 0  !Aleks 04/24
              ntav2_count = 0  !Aleks 04/24

              allocate(dom(ib)%dens(tti,ttj,ttk))
              allocate (dom(ib)%vis(tti,ttj,ttk))

              if (sgs_model==4) then
              allocate(dom(ib)%ksgs(tti,ttj,ttk))
              allocate(dom(ib)%ksgso(tti,ttj,ttk))
              allocate(dom(ib)%eps(tti,ttj,ttk))
              allocate(dom(ib)%epso(tti,ttj,ttk))
              end if
              allocate (dom(ib)%stfcinf(6,pl,ngg))
              if (LENERGY) then
              allocate(dom(ib)%T(tti,ttj,ttk),dom(ib)%To(tti,ttj,ttk))
              allocate(dom(ib)%Tm(tti,ttj,ttk),dom(ib)%Ttm(tti,ttj,ttk))
              allocate(dom(ib)%mu(tti,ttj,ttk))
              end if
              if (LSCALAR) then
              allocate(dom(ib)%S(tti,ttj,ttk),dom(ib)%Sm(tti,ttj,ttk))
              allocate(dom(ib)%So(tti,ttj,ttk),dom(ib)%Stm(tti,ttj,ttk))
              allocate(dom(ib)%sfactor(tti,ttj,ttk))
              end if
              if (L_LSM) then
              allocate(dom(ib)%dens_mg(dom(ib)%tot))
              allocate(dom(ib)%mu(tti,ttj,ttk))
              end if
              if (LAS)   allocate(dom(ib)%mu(tti,ttj,ttk))
              if (differencing==3) allocate(dom(ib)%d1(tti,ttj,ttk), &
        dom(ib)%dphi_dxplus(tti,ttj,ttk), &
        dom(ib)%dphi_dxminus(tti,ttj,ttk), &
        dom(ib)%dphi_dyplus(tti,ttj,ttk), &
        dom(ib)%dphi_dyminus(tti,ttj,ttk), &
        dom(ib)%dphi_dzplus(tti,ttj,ttk), &
        dom(ib)%dphi_dzminus(tti,ttj,ttk))

              allocate (dom(ib)%tauwe(ttj,ttk))
              allocate (dom(ib)%tauww(ttj,ttk))
              allocate (dom(ib)%tauwn(tti,ttk))
              allocate (dom(ib)%tauws(tti,ttk))
              allocate (dom(ib)%tauwt(tti,ttj))
              allocate (dom(ib)%tauwb(tti,ttj))
              allocate (dom(ib)%tauwe2(ttj,ttk))
              allocate (dom(ib)%tauww2(ttj,ttk))
              allocate (dom(ib)%tauwn2(tti,ttk))
              allocate (dom(ib)%tauws2(tti,ttk))
              allocate (dom(ib)%tauwt2(tti,ttj))
              allocate (dom(ib)%tauwb2(tti,ttj))

              if(solver==2) then
              allocate (dom(ib)%faz(ngrd_gl))
              end if

              allocate (dom(ib)%sendb_m1(ngg))
              allocate (dom(ib)%sendb_p1(ngg))
              allocate (dom(ib)%recvb_m1(ngg))
              allocate (dom(ib)%recvb_p1(ngg))
              allocate (dom(ib)%sendb_m2(ngg))
              allocate (dom(ib)%sendb_p2(ngg))
              allocate (dom(ib)%recvb_m2(ngg))
              allocate (dom(ib)%recvb_p2(ngg))
              allocate (dom(ib)%sendb_m3(ngg))
              allocate (dom(ib)%sendb_p3(ngg))
              allocate (dom(ib)%recvb_m3(ngg))
              allocate (dom(ib)%recvb_p3(ngg))

              allocate (dom(ib)%sc1m(ngc),dom(ib)%sc1p(ngc))
              allocate (dom(ib)%rc1m(ngc),dom(ib)%rc1p(ngc))
              allocate (dom(ib)%sc2m(ngc),dom(ib)%sc2p(ngc))
              allocate (dom(ib)%rc2m(ngc),dom(ib)%rc2p(ngc))
              allocate (dom(ib)%sc3m(ngc),dom(ib)%sc3p(ngc))
              allocate (dom(ib)%rc3m(ngc),dom(ib)%rc3p(ngc))
              allocate (dom(ib)%sc4m(ngc),dom(ib)%sc4p(ngc))
              allocate (dom(ib)%rc4m(ngc),dom(ib)%rc4p(ngc))

              allocate (dom(ib)%se1m(nge),dom(ib)%se1p(nge))
              allocate (dom(ib)%re1m(nge),dom(ib)%re1p(nge))
              allocate (dom(ib)%se2m(nge),dom(ib)%se2p(nge))
              allocate (dom(ib)%re2m(nge),dom(ib)%re2p(nge))
              allocate (dom(ib)%se3m(nge),dom(ib)%se3p(nge))
              allocate (dom(ib)%re3m(nge),dom(ib)%re3p(nge))
              allocate (dom(ib)%se4m(nge),dom(ib)%se4p(nge))
              allocate (dom(ib)%re4m(nge),dom(ib)%re4p(nge))
              allocate (dom(ib)%se5m(nge),dom(ib)%se5p(nge))
              allocate (dom(ib)%re5m(nge),dom(ib)%re5p(nge))
              allocate (dom(ib)%se6m(nge),dom(ib)%se6p(nge))
              allocate (dom(ib)%re6m(nge),dom(ib)%re6p(nge))


              dom(ib)%x(1)=dom(ib)%xsl +(-pl+1)*dom(ib)%dx
              dom(ib)%y(1)=dom(ib)%ysl +(-pl+1)*dom(ib)%dy
              dom(ib)%z(1)=dom(ib)%zsl +(-pl+1)*dom(ib)%dz
              dom(ib)%xc(1)=dom(ib)%x(1)-0.5_dp*dom(ib)%dx
              dom(ib)%yc(1)=dom(ib)%y(1)-0.5_dp*dom(ib)%dy
              dom(ib)%zc(1)=dom(ib)%z(1)-0.5_dp*dom(ib)%dz

              do i=2,dom(ib)%ttc_i
                  dom(ib)%x(i)=dom(ib)%x(i-1)+dom(ib)%dx
              end do
              do i=2,dom(ib)%ttc_j
                  dom(ib)%y(i)=dom(ib)%y(i-1)+dom(ib)%dy
              end do
              do i=2,dom(ib)%ttc_k
                  dom(ib)%z(i)=dom(ib)%z(i-1)+dom(ib)%dz
              end do

              do i=1,dom(ib)%ttc_i
                  dom(ib)%xc(i)=dom(ib)%x(i)-0.5_dp*dom(ib)%dx
              end do
              do i=1,dom(ib)%ttc_j
                  dom(ib)%yc(i)=dom(ib)%y(i)-0.5_dp*dom(ib)%dy
              end do
              do i=1,dom(ib)%ttc_k
                  dom(ib)%zc(i)=dom(ib)%z(i)-0.5_dp*dom(ib)%dz
              end do

              if (dom(ib)%inext>=0) then
              if(abs(dom(ib)%x(dom(ib)%iep)-dom(ib)%xel)>1e-5_dp) then
              print*,"mycpu#:",myrank," error-11"
              stop
              end if
              end if

              if (dom(ib)%jnext>=0) then
              if(abs(dom(ib)%y(dom(ib)%jep)-dom(ib)%yel)>1e-5_dp) then
              print*,"mycpu#:",myrank," error-12"
              stop
              end if
              end if

              if (dom(ib)%knext>=0) then
              if(abs(dom(ib)%z(dom(ib)%kep)-dom(ib)%zel)>1e-5_dp) then
              print*,"mycpu#:",myrank," error-13"
              stop
              end if
              end if

              if(solver==2 .and. ngrd_gl>=2) then

              is=dom(ib)%isp
              ie=dom(ib)%iep
              js=dom(ib)%jsp
              je=dom(ib)%jep
              ks=dom(ib)%ksp
              ke=dom(ib)%kep

              do glevel=2,ngrd_gl
                  if(glevel>dom(ib)%ngrid) then
                  gl=dom(ib)%ngrid
                  else
                  gl=glevel
                  end if
                  mgc_i=(ie-is+1)/2**(gl-1)+2
                  mgc_j=(je-js+1)/2**(gl-1)+2
                  mgc_k=(ke-ks+1)/2**(gl-1)+2

                  ndx=dom(ib)%dx*2**(gl-1)
                  ndy=dom(ib)%dy*2**(gl-1)
                  ndz=dom(ib)%dz*2**(gl-1)

                  nwxend=dom(ib)%x(dom(ib)%isp-1)+ndx*(mgc_i-2)
                  nwyend=dom(ib)%y(dom(ib)%jsp-1)+ndy*(mgc_j-2)
                  nwzend=dom(ib)%z(dom(ib)%ksp-1)+ndz*(mgc_k-2)

                  if((abs(dom(ib)%x(dom(ib)%iep)-nwxend)>1E-8_dp) &
            .or.(abs(dom(ib)%y(dom(ib)%jep)-nwyend)>1E-8_dp) &
            .or.(abs(dom(ib)%z(dom(ib)%kep)-nwzend)>1E-8_dp)) then
                  print*,"==ERROR==> in multigrid: max ngrid value"
                  stop
                  end if
              end do
              end if

          end do


      end subroutine initial
!##########################################################################
      subroutine iniflux
!##########################################################################
          use vars
          use multidata
          use multiflow3d_mpi
          use, intrinsic :: iso_fortran_env, only: dp => real64

          implicit none
          integer :: i,j,k,ib,ispr,iepr,jspr,jepr,kspr,kepr
          real(dp) :: buffer_flomas

          MPI_FLT = MPI_DOUBLE_PRECISION

          flomas=0.0_dp

          do ib=1,nbp
              ispr=pl+1
              iepr=dom(ib)%ttc_i-pl
              jspr=pl+1
              jepr=dom(ib)%ttc_j-pl
              kspr=pl+1
              kepr=dom(ib)%ttc_k-pl

              if(dom(ib)%iprev<0) then
              if (dom(ib)%bc_west<61 .and. dom(ib)%bc_west/=4) then
              do j=jspr,jepr
                  do k=kspr,kepr
                      if (L_LSM) then
                      if (dom(ib)%phi(dom(ib)%isu-1,j,k) >= 0.0_dp) then
                      flomas=flomas+dom(ib)%u(dom(ib)%isu-1,j,k)* &
                dom(ib)%dy*dom(ib)%dz
                      end if
                      else
                      flomas=flomas+dom(ib)%u(dom(ib)%isu-1,j,k)* &
                dom(ib)%dy*dom(ib)%dz
                      end if
                  end do
              end do
              end if
              end if

              if(dom(ib)%jprev<0) then
              if (dom(ib)%bc_south<61 .and. dom(ib)%bc_south/=4) then
              do k=kspr,kepr
                  do i=ispr,iepr
                      if (L_LSM) then
                      if (dom(ib)%phi(i,dom(ib)%jsv-1,k) >= 0.0_dp) then
                      flomas=flomas+dom(ib)%v(i,dom(ib)%jsv-1,k)* &
                dom(ib)%dx*dom(ib)%dz
                      end if
                      else
                      flomas=flomas+dom(ib)%v(i,dom(ib)%jsv-1,k)* &
                dom(ib)%dx*dom(ib)%dz
                      end if
                  end do
              end do
              end if
              end if

              if(dom(ib)%kprev<0) then
              if (dom(ib)%bc_bottom<61.and.dom(ib)%bc_bottom/=4) then
              do j=jspr,jepr
                  do i=ispr,iepr
                      if (L_LSM) then
                      if (dom(ib)%phi(i,j,dom(ib)%ksw-1) >= 0.0_dp) then
                      flomas=flomas+dom(ib)%w(i,j,dom(ib)%ksw-1)* &
                dom(ib)%dx*dom(ib)%dy
                      end if
                      else
                      flomas=flomas+dom(ib)%w(i,j,dom(ib)%ksw-1)* &
                dom(ib)%dx*dom(ib)%dy
                      end if
                  end do
              end do
              end if
              end if
          end do

          buffer_flomas = flomas
          call MPI_ALLREDUCE(buffer_flomas,flomas,1,MPI_FLT,MPI_SUM, &
    MPI_COMM_WORLD,ierr)

          return
      end subroutine iniflux
!##########################################################################
      subroutine correctoutflux
!##########################################################################
          use, intrinsic :: iso_fortran_env, only: dp => real64
          use vars
          use multidata
          use multiflow3d_mpi
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          integer :: i,j,k,ib,ispr,iepr,jspr,jepr,kspr,kepr
          real(dp) :: fmout,fct,buffer_fmout

          MPI_FLT = MPI_DOUBLE_PRECISION

          fmout=0.0_dp

          do ib=1,nbp
              ispr=pl+1
              iepr=dom(ib)%ttc_i-pl
              jspr=pl+1
              jepr=dom(ib)%ttc_j-pl
              kspr=pl+1
              kepr=dom(ib)%ttc_k-pl

              if(dom(ib)%inext<0) then
              if (dom(ib)%bc_east<61 .and. dom(ib)%bc_east/=4) then
              do j=jspr,jepr
                  do k=kspr,kepr
                      if (L_LSM) then
                      if (dom(ib)%phi(dom(ib)%ieu+1,j,k) >= 0.0_dp) then
                      fmout=fmout+dom(ib)%u(dom(ib)%ieu+1,j,k)* &
                dom(ib)%dy*dom(ib)%dz
                      end if
                      else
                      fmout=fmout+dom(ib)%u(dom(ib)%ieu+1,j,k)* &
                dom(ib)%dy*dom(ib)%dz
                      end if
                  end do
              end do
              end if
              end if

              if(dom(ib)%jnext<0) then
              if (dom(ib)%bc_north<61 .and. dom(ib)%bc_north/=4) then
              do k=kspr,kepr
                  do i=ispr,iepr
                      if (L_LSM) then
                      if (dom(ib)%phi(i,dom(ib)%jev+1,k) >= 0.0_dp) then
                      fmout=fmout+dom(ib)%v(i,dom(ib)%jev+1,k)* &
                dom(ib)%dx*dom(ib)%dz
                      end if
                      else
                      fmout=fmout+dom(ib)%v(i,dom(ib)%jev+1,k)* &
                dom(ib)%dx*dom(ib)%dz
                      end if
                  end do
              end do
              end if
              end if

              if(dom(ib)%knext<0) then
              if (dom(ib)%bc_top<61 .and. dom(ib)%bc_top/=4) then
              do j=jspr,jepr
                  do i=ispr,iepr
                      if (L_LSM) then
                      if (dom(ib)%phi(i,j,dom(ib)%kew+1) >= 0.0_dp) then
                      fmout=fmout+dom(ib)%w(i,j,dom(ib)%kew+1)* &
                dom(ib)%dx*dom(ib)%dy
                      end if
                      else
                      fmout=fmout+dom(ib)%w(i,j,dom(ib)%kew+1)* &
                dom(ib)%dx*dom(ib)%dy
                      end if
                  end do
              end do
              end if
              end if
          end do

          buffer_fmout = fmout
          call MPI_ALLREDUCE(buffer_fmout,fmout,1,MPI_FLT,MPI_SUM, &
    MPI_COMM_WORLD,ierr)

          fct=flomas/(fmout+1.E-30_dp)

          Mdef=flomas-fmout

          do ib=1,nbp
              if(dom(ib)%inext<0) then
              if (dom(ib)%bc_east<61 .and. dom(ib)%bc_east/=4) then
              do j=dom(ib)%jsu,dom(ib)%jeu
                  do k=dom(ib)%ksu,dom(ib)%keu
                      dom(ib)%u(dom(ib)%ieu+1,j,k)= &
                dom(ib)%u(dom(ib)%ieu+1,j,k)*fct
                  end do
              end do

              do j=dom(ib)%jsv,dom(ib)%jev
                  do k=dom(ib)%ksv,dom(ib)%kev
                      dom(ib)%v(dom(ib)%iev+1,j,k)= &
                dom(ib)%v(dom(ib)%iev+1,j,k)*fct
                  end do
              end do

              do j=dom(ib)%jsw,dom(ib)%jew
                  do k=dom(ib)%ksw,dom(ib)%kew
                      dom(ib)%w(dom(ib)%iew+1,j,k)= &
                dom(ib)%w(dom(ib)%iew+1,j,k)*fct
                  end do
              end do
              end if
              end if

              if(dom(ib)%jnext<0) then
              if (dom(ib)%bc_north<61 .and. dom(ib)%bc_north/=4) then
              do i=dom(ib)%isu,dom(ib)%ieu
                  do k=dom(ib)%ksu,dom(ib)%keu
                      dom(ib)%u(i,dom(ib)%jeu+1,k)= &
                dom(ib)%u(i,dom(ib)%jeu+1,k)*fct
                  end do
              end do

              do i=dom(ib)%isv,dom(ib)%iev
                  do k=dom(ib)%ksv,dom(ib)%kev
                      dom(ib)%v(i,dom(ib)%jev+1,k)= &
                dom(ib)%v(i,dom(ib)%jev+1,k)*fct
                  end do
              end do

              do i=dom(ib)%isw,dom(ib)%iew
                  do k=dom(ib)%ksw,dom(ib)%kew
                      dom(ib)%w(i,dom(ib)%jew+1,k)= &
                dom(ib)%w(i,dom(ib)%jew+1,k)*fct
                  end do
              end do
              end if
              end if

              if(dom(ib)%knext<0) then
              if (dom(ib)%bc_top<61 .and. dom(ib)%bc_top/=4) then
              do i=dom(ib)%isu,dom(ib)%ieu
                  do j=dom(ib)%jsu,dom(ib)%jeu
                      dom(ib)%u(i,j,dom(ib)%keu+1)= &
                dom(ib)%u(i,j,dom(ib)%keu+1)*fct
                  end do
              end do

              do i=dom(ib)%isv,dom(ib)%iev
                  do j=dom(ib)%jsv,dom(ib)%jev
                      dom(ib)%v(i,j,dom(ib)%kev+1)= &
                dom(ib)%v(i,j,dom(ib)%kev+1)*fct
                  end do
              end do

              do i=dom(ib)%isw,dom(ib)%iew
                  do j=dom(ib)%jsw,dom(ib)%jew
                      dom(ib)%w(i,j,dom(ib)%kew+1)= &
                dom(ib)%w(i,j,dom(ib)%kew+1)*fct
                  end do
              end do
              end if
              end if
          end do

          return
      end subroutine correctoutflux
!##########################################################################
      subroutine initflowfield
!##########################################################################
          use, intrinsic :: iso_fortran_env, only: dp => real64
          use vars
          use multiflow3d_mpi
          use multidata
          use module_lsm
          use multiflow3d_sem, only: sem_initial
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          integer :: i,j,k,ib,tti,ttj,ttk,pll
          integer :: sn,sn2
          integer :: inind,jnind,knind
          real(dp) :: dum,ubw,ube,ubs,ubn,ubt,ubb,vb,wb,lz,dummy
          real(dp), dimension(21) :: dm
          character(len=8)   :: chb1
          character(len=25)  :: gf

          dm=0.d0

          do ib=1,nbp

              tti=dom(ib)%ttc_i
              ttj=dom(ib)%ttc_j
              ttk=dom(ib)%ttc_k

              if (LRESTART) then

              qzero=ubulk  !brunho2014
              open (unit=700, file="final_ctime.dat")
              read (700,"(i8,3F15.6)") ntime,ctime,forcn,qstpn,count &
        ,ntav1_count,ntav2_count
              close (700)
              if (.not.reinitmean) then  !Aleks 04/24
              dom(ib)%ntav1=ntav1_count
              dom(ib)%ntav2=ntav2_count
              ntav_restart=ntav2_count
              end if
!===============================================================

              write(chb1,"(i8)") dom_id(ib)
              sn=len(trim(adjustl(chb1)))
              chb1=repeat("0",(4-sn))//trim(adjustl(chb1))
              gf="tecbin"//trim(adjustl(chb1))//".bin"
              open (unit=700, file=gf, form="unformatted",status="old")

              read (700) tti,ttj,ttk
              read (700) pll
              read (700) inind,jnind,knind

              if(pll/=pl .and. myrank==0) then
              print*,"&*&* different number of overlapping layers!!",pl,pll
              write(numfile,*) "&*&* different number of overlapping layers!!"
              stop
              end if

              do k=1,ttk
                  do j=1,ttj
                      do i=1,tti

                          read (700) dummy,dummy,dummy,dm(1),dm(2),dm(3), &
                    dm(4),dm(5),dm(6),dm(7),dm(8),dm(9),dm(10),dm(11),dm(12), &
                    dm(13),dm(14),dm(15),dm(16),dm(17),dm(18),dm(19),dm(20)
!     & dm(21)

                          dom(ib)%p  (i,j,k)=dm(1)
                          dom(ib)%pm (i,j,k)=dm(2)
                          dom(ib)%ppm(i,j,k)=dm(3)
                          dom(ib)%u  (i,j,k)=dm(4)
                          dom(ib)%um (i,j,k)=dm(5)
                          dom(ib)%uum(i,j,k)=dm(6)
                          dom(ib)%v  (i,j,k)=dm(7)
                          dom(ib)%vm (i,j,k)=dm(8)
                          dom(ib)%vvm(i,j,k)=dm(9)
                          dom(ib)%w  (i,j,k)=dm(10)
                          dom(ib)%wm (i,j,k)=dm(11)
                          dom(ib)%wwm(i,j,k)=dm(12)
                          dom(ib)%uvm(i,j,k)=dm(13)
                          dom(ib)%uwm(i,j,k)=dm(14)
                          dom(ib)%vwm(i,j,k)=dm(15)
                          dom(ib)%vis(i,j,k)=dm(16)
                          if (LSCALAR) dom(ib)%S(i,j,k)=dm(17)
                          if (LSCALAR) dom(ib)%Sm(i,j,k) = dm(18)
                          if (LENERGY) dom(ib)%T(i,j,k)=dm(19)
                          if (LENERGY) dom(ib)%Tm(i,j,k)=dm(20)
                          !if (LENERGY) dom(ib)%Ttm(i,j,k)=dm(21)

                      end do
                  end do
              end do
              close (700)
!===============================================================

              if (reinitmean) then
              dom(ib)%um   = 0.0_dp
              dom(ib)%vm   = 0.0_dp
              dom(ib)%wm   = 0.0_dp
              dom(ib)%pm   = 0.0_dp
              dom(ib)%uum  = 0.0_dp
              dom(ib)%vvm  = 0.0_dp
              dom(ib)%wwm  = 0.0_dp
              dom(ib)%uvm  = 0.0_dp
              dom(ib)%uwm  = 0.0_dp
              dom(ib)%vwm  = 0.0_dp
              dom(ib)%ppm  = 0.0_dp
              dom(ib)%Tm   = 0.0_dp
              dom(ib)%Ttm  = 0.0_dp
              ctime=0.0_dp
              ntime=0
              if (L_LSM) dom(ib)%phim  = 0.0_dp
              end if

              else  !no restart !cold initialisation


              qzero=ubulk                               !brunho2014
              qstpn=qzero
              forcn=2.0_dp/(Re*qzero)
              ctime=0.0_dp
              ntime=0
              dom(ib)%u=Ubulk
              dom(ib)%uo=Ubulk
              dom(ib)%uoo=Ubulk

              lz=zen-zst

              if (L_LSM) call init_lsm

              dom(ib)%p=0.0_dp

!======================STRATIFICATION CONDITIONS========================
!     if (LAS) then
!                do k=1,ttk
!                  do j=1,ttj
!                    do i=1,tti
!       if (dom(ib)%z(k).gt.0.66) then
!              dom(ib)%T(i,j,k)=5.;  dom(ib)%To(i,j,k)=5.
!       elseif (dom(ib)%z(k).le.0.66.and.dom(ib)%z(k).gt.0.33) then
!              dom(ib)%T(i,j,k)=0.;  dom(ib)%To(i,j,k)=0.
!       elseif (dom(ib)%z(k).le.0.33) then
!              dom(ib)%T(i,j,k)=-5.;  dom(ib)%To(i,j,k)=-5.
!             end do
!             end do
!           end do
!     endif
!=======================================================================

              dom(ib)%v=0.0_dp
              dom(ib)%w=0.0_dp
              dom(ib)%vo=0.0_dp
              dom(ib)%voo=0.0_dp
              dom(ib)%wo=0.0_dp
              dom(ib)%woo=0.0_dp

              dom(ib)%dens=dens
              dom(ib)%vis=rrey

              if (LENERGY) then
              dom(ib)%T=Tinit
              dom(ib)%To=Tinit
              dom(ib)%Tm=0.0_dp
              dom(ib)%Ttm=0.0_dp
              dom(ib)%mu=rrey*dens
              call energy_init
              end if
              if (LSCALAR) then
              dom(ib)%S=0.0_dp
              dom(ib)%So=0.0_dp
              dom(ib)%Sm=0.0_dp
              dom(ib)%Stm=0.0_dp
              call sediment_init
              end if
              if (L_LSM)          dom(ib)%mu=rrey*dens
              if (LAS) then
              dom(ib)%mu=rrey*dens
              call Active_scalar
              end if
              if (LNonNewt) call NonNewtonian

              dom(ib)%um   = 0.0_dp
              dom(ib)%vm   = 0.0_dp
              dom(ib)%wm   = 0.0_dp
              dom(ib)%pm   = 0.0_dp
              dom(ib)%uum  = 0.0_dp
              dom(ib)%vvm  = 0.0_dp
              dom(ib)%wwm  = 0.0_dp
              dom(ib)%uvm  = 0.0_dp
              dom(ib)%uwm  = 0.0_dp
              dom(ib)%vwm  = 0.0_dp
              dom(ib)%ppm  = 0.0_dp

              dom(ib)%tauww  = 0.0_dp
              dom(ib)%tauww2  = 0.0_dp
              dom(ib)%tauwe  = 0.0_dp
              dom(ib)%tauwe2  = 0.0_dp
              dom(ib)%tauws  = 0.0_dp
              dom(ib)%tauws2  = 0.0_dp
              dom(ib)%tauwn  = 0.0_dp
              dom(ib)%tauwn2  = 0.0_dp
              dom(ib)%tauwb  = 0.0_dp
              dom(ib)%tauwb2  = 0.0_dp
              dom(ib)%tauwt  = 0.0_dp
              dom(ib)%tauwt2  = 0.0_dp

              if (sgs_model>2) then
              dom(ib)%ksgs = (3.d0/2.d0)*(ubulk*0.1_dp)**2.0_dp
              dom(ib)%eps  = 0.09_dp**0.75_dp*dom(ib)%ksgs**1.5_dp/(0.07_dp*lz)
              dom(ib)%ksgso = (3.d0/2.d0)*(ubulk*0.1_dp)**2.0_dp
              dom(ib)%epso  = 0.09_dp**0.75_dp*dom(ib)%ksgs**1.5_dp/(0.07_dp*lz)
              end if

              if (trim(keyword)=="channel") then
              if (.not.L_LSM) dom(ib)%u=ubulk
              ubw=ubulk
              ube=ubulk
              ubs=ubulk               !brunho2014
              ubn=ubulk
              ubt=ubulk
              ubb=ubulk
              vb=0.0_dp
              wb=0.0_dp
              else if (trim(keyword)=="cavity") then
              dom(ib)%u=0.0_dp
              ubw=0.0_dp
              ube=0.0_dp
              ubs=0.0_dp
              ubn=2.0_dp
              ubt=0.0_dp
              ubb=0.0_dp
              vb=0.0_dp
              wb=0.0_dp
              else if (trim(keyword)=="column") then
              dom(ib)%u=0.0_dp
              ubw=0.0_dp
              ube=0.0_dp
              ubs=0.0_dp
              ubn=0.0_dp
              ubt=0.0_dp
              ubb=0.0_dp
              vb=0.0_dp
              wb=0.0_dp
              else
              write (6,*) " wrong keyword "
              end if

!..............U=> West and East ...............
              if (dom(ib)%iprev<0) then
              do k=1,ttk
                  do j=1,ttj
                      if (L_LSM) then                                        !I deleted 'or LSM_BASE'
                      if (dom(ib)%zc(k)>length) then
                      dom(ib)%u(dom(ib)%isu-1,j,k) = 0.0_dp
                      end if
                      else
                      dom(ib)%u(dom(ib)%isu-1,j,k) = ubw
                      end if
                  end do
              end do
              end if
              if (dom(ib)%inext<0) then
              do k=1,ttk
                  do j=1,ttj
                      if (L_LSM) then
                      if (dom(ib)%zc(k)>length) then
                      dom(ib)%u(dom(ib)%ieu+1,j,k) = 0.0_dp
                      end if
                      else
                      dom(ib)%u(dom(ib)%ieu+1,j,k) = ube
                      end if
                  end do
              end do
              end if
!.............U=> South and North .................
              if (dom(ib)%jprev<0) then
              do k=1,ttk
                  do i=1,tti
                      if (L_LSM) then
                      if (dom(ib)%zc(k)>length) then
                      dom(ib)%u(i,dom(ib)%jsu-1,k) = 0.0_dp
                      end if
                      else
                      dom(ib)%u(i,dom(ib)%jsu-1,k) = ubs
                      end if
                  end do
              end do
              end if
              if (dom(ib)%jnext<0) then
              do k=1,ttk
                  do i=1,tti
                      if (L_LSM) then
                      if (dom(ib)%zc(k)>length) then
                      dom(ib)%u(i,dom(ib)%jeu+1,k) = 0.0_dp
                      end if
                      else
                      dom(ib)%u(i,dom(ib)%jeu+1,k) = ubn
                      end if
                  end do
              end do
              end if
!.............U=> Bottom and Top .................
              if (dom(ib)%kprev<0) then
              do j=1,ttj
                  do i=1,tti
                      dom(ib)%u(i,j,dom(ib)%ksu-1) = ubb
                  end do
              end do
              end if
              if (dom(ib)%knext<0) then
              do j=1,ttj
                  do i=1,tti
                      if (L_LSM) then
                      dom(ib)%u(i,j,dom(ib)%keu+1) = 0.0_dp
                      else
                      dom(ib)%u(i,j,dom(ib)%keu+1) = ubt
                      end if
                  end do
              end do
              end if
!........... V=> West and East ....................
              if (dom(ib)%iprev<0) then
              do k=1,ttk
                  do j=1,ttj
                      dom(ib)%v(dom(ib)%isv-1,j,k)    = vb
                  end do
              end do
              end if
              if (dom(ib)%inext<0) then
              do k=1,ttk
                  do j=1,ttj
                      dom(ib)%v(dom(ib)%iev+1,j,k)  = vb
                  end do
              end do
              end if
!............V=> South and North ....................
              if (dom(ib)%jprev<0) then
              do k=1,ttk
                  do i=1,tti
                      dom(ib)%v(i,dom(ib)%jsv-1,k)    = vb
                  end do
              end do
              end if
              if (dom(ib)%jnext<0) then
              do k=1,ttk
                  do i=1,tti
                      dom(ib)%v(i,dom(ib)%jev+1,k)  = vb
                  end do
              end do
              end if
!.............V=> Bottom and Top .................
              if (dom(ib)%kprev<0) then
              do j=1,ttj
                  do i=1,tti
                      dom(ib)%v(i,j,dom(ib)%ksv-1)    = vb
                  end do
              end do
              end if
              if (dom(ib)%knext<0) then
              do j=1,ttj
                  do i=1,tti
                      dom(ib)%v(i,j,dom(ib)%kev+1)  = vb
                  end do
              end do
              end if
!........... W=> West and East ....................
              if (dom(ib)%iprev<0) then
              do k=1,ttk
                  do j=1,ttj
                      dom(ib)%w(dom(ib)%isw-1,j,k)    = wb
                  end do
              end do
              end if
              if (dom(ib)%inext<0) then
              do k=1,ttk
                  do j=1,ttj
                      dom(ib)%w(dom(ib)%iew+1,j,k)  = wb
                  end do
              end do
              end if
!............W=> South and North ....................
              if (dom(ib)%jprev<0) then
              do k=1,ttk
                  do i=1,tti
                      dom(ib)%w(i,dom(ib)%jsw-1,k)    = wb
                  end do
              end do
              end if
              if (dom(ib)%jnext<0) then
              do k=1,ttk
                  do i=1,tti
                      dom(ib)%w(i,dom(ib)%jew+1,k)  = wb
                  end do
              end do
              end if
!.............W=> Bottom and Top .................
              if (dom(ib)%kprev<0) then
              do j=1,ttj
                  do i=1,tti
                      dom(ib)%w(i,j,dom(ib)%ksw-1)    = wb
                  end do
              end do
              end if
              if (dom(ib)%knext<0) then
              do j=1,ttj
                  do i=1,tti
                      dom(ib)%w(i,j,dom(ib)%kew+1)  = wb
                  end do
              end do
              end if


!.######### U=> When power law inlet condition, 7 Dic 2015 .##########
!          IF (dom(ib)%bc_west.eq.12 .or. UPROF_SEM.eq.12) THEN
              IF (dom(ib)%bc_west==12) THEN
              do i = dom(ib)%isu-1,dom(ib)%ieu+1
                  do j = dom(ib)%jsu-1,dom(ib)%jeu+1
                      do k = dom(ib)%ksu-1,dom(ib)%keu+1
                          if (dom(ib)%yc(j)<((yen-yst)/2)) then
                          dom(ib)%u(i,j,k) = ubulk*(1.0d0+1.0d0/7.0d0) &
                    *(DABS(2*dom(ib)%yc(j)/(yen-yst)))**(1.d0/7.d0)
                          else
                          dom(ib)%u(i,j,k) = ubulk*(1.0d0+1.0d0/7.0d0) &
                    *(DABS(2*((yen-yst)-dom(ib)%yc(j))/(yen-yst)))**(1.d0/7.d0)
                          end if
                          dom(ib)%u(i,j,k) = dom(ib)%u(i,j,k)*(1.0d0+1.0d0/7.0d0) &
                    *(DABS(dom(ib)%zc(k)/(zen-zst)))**(1.d0/7.d0)
                      end do
                      end do
                      end do
              END IF
!.######### U=> When power law inlet condition, 7 Dic 2015 .##########
              IF (dom(ib)%bc_west==13) THEN
              do i = dom(ib)%isu-1,dom(ib)%ieu+1
                  do j = dom(ib)%jsu-1,dom(ib)%jeu+1
                      do k = dom(ib)%ksu-1,dom(ib)%keu+1
                          if (dom(ib)%yc(j)<((yen-yst)/2)) then
                          dom(ib)%u(i,j,k) = ubulk*(1.0d0+1.0d0/7.0d0) &
                    *(DABS(2*dom(ib)%yc(j)/(yen-yst)))**(1.d0/7.d0)
                          else
                          dom(ib)%u(i,j,k) = ubulk*(1.0d0+1.0d0/7.0d0) &
                    *(DABS(2*((yen-yst)-dom(ib)%yc(j))/(yen-yst)))**(1.d0/7.d0)
                          end if
                      end do
                      end do
                      end do
              END IF

              end if    !No restart

!        Allocate time series

              jtime=itime_end-ntime

              if (ntime*dt<t_start_averaging2) then
              jtime=itime_end-INT(t_start_averaging2/dt)+1
              end if

              allocate(dom(ib)%u_unst(n_unstpt,jtime))
              allocate(dom(ib)%v_unst(n_unstpt,jtime))
              allocate(dom(ib)%w_unst(n_unstpt,jtime))
              allocate(dom(ib)%um_unst(n_unstpt,jtime))
              allocate(dom(ib)%vm_unst(n_unstpt,jtime))
              allocate(dom(ib)%wm_unst(n_unstpt,jtime))
              allocate(dom(ib)%p_unst(n_unstpt,jtime))
              allocate(dom(ib)%pm_unst(n_unstpt,jtime))
              allocate(dom(ib)%ksgs_unst(n_unstpt,jtime))
              allocate(dom(ib)%eps_unst(n_unstpt,jtime))
              allocate(dom(ib)%T_unst(n_unstpt,jtime))  !Aleks 04/24
              allocate(dom(ib)%Tm_unst(n_unstpt,jtime))  !Aleks 04/24

          end do

!.################################################
!.###########  Synthetic Eddy Method, 14 Dic 2015    ##########
          IF (bc_w==8 .and. myrank==0) then
          print*,"Writing the SEM inlet"
          call sem_initial()  !Generate the files for the inlet turbulent field
          print*,"Finish the SEM inlet"
          END IF

   70     format (10e25.8)
   71     format (3F15.6)

      end subroutine initflowfield
!##########################################################################
