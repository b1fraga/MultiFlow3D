!######################################################################
      subroutine init_lsm
!######################################################################
          use vars
          use module_LSM
          use multidata
          use multiflow3d_mpi
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          integer :: i,j,k,ib,tti,ttj,ttk
          integer :: glevel,gl,mgc_i,mgc_j,mgc_k

!READING
          open (unit=13, file="input/lsm.cin")
          read (12,*)
          read (12,*) reinit,ntime_reinit,reldif_LSM,length,accuracy &
    ,cfl_lsm
          read (12,*) LENDS
          read (12,*) L_LSMinit
          read (12,*) L_anim_phi,L_anim_grd
          read (12,*) densl,densg,nul,nug
          read (12,*) slope

          close(13)
!WARNINGS
          if (L_LSMinit .and. (L_anim_phi .or. L_anim_grd)) then
          if (myrank==0) then
          print*,"Error: not possible to output animation files", &
    "  for LSM_init run!"
          end if
          stop
          end if

          if (L_LSMbase .and. L_LSMinit) then
          if (myrank==0) then
          print*,"Error: L_LSMbase and L_LSMinit cannot both be true!"
          end if
          stop
          end if

          if (L_LSMbase .and. L_LSM) then
          if (myrank==0) then
          print*,"Error: L_LSMbase and L_LSM cannot both be true!"
          end if
          stop
          end if

          if (L_LSMinit .and. (.not.L_LSM)) then
          if (myrank==0) then
          print*,"Error: L_LSMinit cannot be true if L_LSM is false!"
          end if
          stop
          end if

          if (L_anim_phi .and. (.not.L_LSM)) then
          if (myrank==0) then
          print*,"Error: L_anim_phi cannot be true if L_LSM is false!"
          end if
          stop
          end if
!ALLOCATIONS
          allocate(dom(ib)%ijkp_lsm(0:dom(ib)%ngrid))
          allocate (dom(ib)%dens_mg(dom(ib)%tot))

          dom(ib)%ijkp_lsm = 0
          dom(ib)%ijkp_lsm(1)=(dom(ib)%ttc_i-2*pl)* &
    (dom(ib)%ttc_j-2*pl)*(dom(ib)%ttc_k-2*pl)
          do glevel=2,dom(ib)%ngrid
              mgc_i=(dom(ib)%iep-dom(ib)%isp+1)/2**(glevel-1)+2
              mgc_j=(dom(ib)%jep-dom(ib)%jsp+1)/2**(glevel-1)+2
              mgc_k=(dom(ib)%kep-dom(ib)%ksp+1)/2**(glevel-1)+2
              dom(ib)%ijkp_lsm(glevel)=dom(ib)%ijkp_lsm(glevel-1)+ &
        (mgc_i-2)*(mgc_j-2)*(mgc_k-2)
          end do
          dom(ib)%tot=dom(ib)%ijkp_lsm(dom(ib)%ngrid)
!INITIALISATIONS
          if (L_LSMbase) then
          do k=2,ttk
              do j=1,ttj
                  do i=1,tti
                      if (dom(ib)%z(k-1)<=length) then
                      dom(ib)%u(i,j,k)=Ubulk
                      dom(ib)%uo(i,j,k)=Ubulk
                      dom(ib)%uoo(i,j,k)=Ubulk
                      else
                      dom(ib)%u(i,j,k)=0.0_dp
                      dom(ib)%uo(i,j,k)=0.0_dp
                      dom(ib)%uoo(i,j,k)=0.0_dp
                      end if
                  end do
              end do
          end do
          else if (L_LSM) then
          do k=2,ttk
              do j=1,ttj
                  do i=1,tti
                      if (dom(ib)%phi(i,j,k)>=0.0_dp) then
                      dom(ib)%u(i,j,k)=Ubulk
                      dom(ib)%uo(i,j,k)=Ubulk
                      dom(ib)%uoo(i,j,k)=Ubulk
                      else
                      dom(ib)%u(i,j,k)=0.0_dp
                      dom(ib)%uo(i,j,k)=0.0_dp
                      dom(ib)%uoo(i,j,k)=0.0_dp
                      end if
                  end do
                  end do
                  end do
          end if

          mul = nul * densl
          mug = nug * densg

      end subroutine init_lsm

! !######################################################################
!       subroutine initial_lsm_3d_channel
! !######################################################################
!       use vars
!       use module_LSM
!       use multidata
!       use multiflow3d_mpi

!       implicit none
!       integer :: i,j,k,ib,tti,ttj,ttk,sn,sn1
!       integer :: is,ie,js,je,ks,ke
!       real(dp) :: b,dummy
!       character*8 :: chb,chb1
!       character*31 :: gridfile
!       character *80 dummyline

!   if (myrank.eq.0) then
!     write(*,'(a)') '**********************************************'
!     write(*,'(a)') '*'
!     write(*,'(a)') '*            TWO-PHASE SIMULATION'
!     if (l_lsmbase) then
!       write(*,'(a,f8.3)') '* Base case with rigid lid at z=',length
!     else
!       write(*,'(a,f8.3)') '* Initial water level: z=',length
!     endif
!     write(*,'(a)') '*'
!     write(*,'(a)') '**********************************************'
!   endif

!       do ib=1,nbp
!         tti=dom(ib)%ttc_i;  ttj=dom(ib)%ttc_j;  ttk=dom(ib)%ttc_k
!         is=dom(ib)%isp; ie=dom(ib)%iep
!         js=dom(ib)%jsp; je=dom(ib)%jep
!         ks=dom(ib)%ksp; ke=dom(ib)%kep

!         allocate (dom(ib)%phi(tti,ttj,ttk),
!      &dom(ib)%phi_init(tti,ttj,ttk),dom(ib)%phi_reinit(tti,ttj,ttk),
!      &dom(ib)%phi_new(tti,ttj,ttk),dom(ib)%dphi_dx(tti,ttj,ttk),
!      &dom(ib)%dphi_dy(tti,ttj,ttk),dom(ib)%dphi_dz(tti,ttj,ttk),
!      &dom(ib)%s_phi0(tti,ttj,ttk),dom(ib)%h_phi(tti,ttj,ttk),
!      &dom(ib)%phim(tti,ttj,ttk))
! !
! ! Read in phi field if restarting from previous solution
! !
!         if (L_LSM .and. lrestart) then
!           write(chb,'(i8)') dom_id(ib)
!           write(chb1,'(i8)') itime_start
!           sn=len(trim(adjustl(chb)))
!           sn1=len(trim(adjustl(chb1)))
!           chb=repeat('0',(4-sn))//trim(adjustl(chb))
!           chb1=repeat('0',(6-sn1))//trim(adjustl(chb1))
!           gridfile='tecout_phi_'//trim(adjustl(chb))//'_'//
!      & 'initial'//'.dat'
!           open (unit=703, file=gridfile)
!           read (703,*) dummyline
!           read (703,*) dummyline
!           read (703,*) dummyline
!           read (703,*) dummyline
!           read (703,*) dummyline
!           do k=ks-1,ke
!             do j=js-1,je
!               do i=is-1,ie
!                 read (703,73) dummy,dummy,dummy,
!      & dom(ib)%phi(i,j,k),dom(ib)%phim(i,j,k),
!      & dom(ib)%dens(i,j,k),dom(ib)%mu(i,j,k)
!               end do
!             end do
!           end do
!  73       format (10e25.8)
!           close(703)
!           dom(ib)%phi_init = 0.0_dp
!           dom(ib)%phi_new = 0.0_dp
!           dom(ib)%phi_reinit = 0.0_dp
!         else
! !
! ! Initialise uniform phi field, if not restarting
! !
!           dom(ib)%phi=1.0_dp
!           dom(ib)%phi_init = 0.0_dp
!           dom(ib)%phi_new = 0.0_dp
!           dom(ib)%phi_reinit = 0.0_dp

!           if (trim(keyword).eq.'channel') then      ! Channel flow case
!             do k=1,ttk
!               do j=1,ttj
!                 do i=1,tti
! !
! ! Set initial free surface profile, if not uniform (this is case dependent)
! ! Initialise length (i.e. distance of free surface from bed)
! !========== Richard cube test =========================================
! !          if ((dom(ib)%xc(i).ge.0).and.(dom(ib)%xc(i).le.0.035)) then
! !            length=0.025
! !          else if ((dom(ib)%xc(i).gt.0.035).and.
! !     &             (dom(ib)%xc(i).le.0.065)) then
! !            length=-0.266666667*dom(ib)%xc(i)+0.0343
! !          else if ((dom(ib)%xc(i).gt.0.065).and.
! !     &             (dom(ib)%xc(i).le.0.25088)) then
! !            length=0.017
! !          end if
! !========== Sibel constriction test ===================================
! !          if ((dom(ib)%xc(i).ge.0).and.(dom(ib)%xc(i).le.0.59)) then
! !            length=0.076
! !          else if ((dom(ib)%xc(i).gt.0.59).and.
! !     &             (dom(ib)%xc(i).le.0.885)) then
! !            length=-0.213559322*dom(ib)%xc(i)+0.202
! !          else if ((dom(ib)%xc(i).gt.0.885).and.
! !     &             (dom(ib)%xc(i).le.1.475)) then
! !            length=0.013
! !          end if
! !======================================================================
! !
! ! Initialise phi (free surface defined by phi=0, phi=-ve above, +ve below)
! !
!                   if (dom(ib)%zc(k).lt.length)   then
!                     dom(ib)%phi(i,j,k) = 1.0*abs(dom(ib)%zc(k)-length)
!                   else if (dom(ib)%zc(k).gt.length)   then
!                     dom(ib)%phi(i,j,k) = -1.0*abs(dom(ib)%zc(k)-length)
!                   else if (dom(ib)%zc(k).eq.length)  then
!                     dom(ib)%phi(i,j,k) = 0.0_dp
!                   end if
!                   dom(ib)%s_phi0(i,j,k)=dom(ib)%phi(i,j,k)
!                   dom(ib)%phi_reinit(i,j,k)=dom(ib)%phi(i,j,k)
!                 end do
!               end do
!             end do

!           else if (trim(keyword).eq.'wave') then   ! Solitary wave case

!             do k=1,ttk
!               do j=1,ttj
!                 do i=1,tti
!                   b=length/(cosh(sqrt(3.0_dp*length)/2.0_dp*(dom(ib)%xc(i))))**2
!                   if (dom(ib)%zc(k).lt.(b+1))   then
!                     dom(ib)%phi(i,j,k) = 1.0_dp*abs(dom(ib)%zc(k)-(b+1))
!                   else if (dom(ib)%zc(k).gt.(b+1))   then
!                     dom(ib)%phi(i,j,k) = -1.0_dp*abs(dom(ib)%zc(k)-(b+1))
!                   else if (dom(ib)%zc(k).eq.(b+1))  then
!                     dom(ib)%phi(i,j,k) = 0.0_dp
!                   end if
!                   dom(ib)%s_phi0(i,j,k)=dom(ib)%phi(i,j,k)
!                 end do
!               end do
!             end do

!           end if

!         end if

!         dt_reinit = cfl_lsm*max(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz)

!       end do
! !
! ! Set level set function phi to a signed distance function
! !
!       if (reinit)  call tvd_rk_reinit
! !
! ! Define density and viscosity above, below and across the surface
! !
!       call heaviside

!       return
!       end subroutine initial_lsm_3d_channel
!######################################################################
      subroutine lsm_3d
!######################################################################
          use vars
          use module_LSM
          use multidata

          implicit none
          integer :: i,j,k,ib

          call tvd_rk_3step

          return
      end subroutine lsm_3d
!######################################################################
      subroutine  tvd_rk_3step
!######################################################################
!**********************************************************************
! 3 step runge kutta routine to return convected phi field
!**********************************************************************
          use vars
          use module_LSM
          use multidata
          use multiflow3d_mpi
          use, intrinsic :: iso_fortran_env, only: dp => real64

          implicit none
          integer :: i,j,k,ib
          real(dp) :: uijk,vijk,wijk,h1,h2,h3
!
! First RK step (phi --> phi_new)
!
          call dphi_a_v_3d(14)

          do ib=1,nbp

              do k=dom(ib)%ksp,dom(ib)%kep
                  do i=dom(ib)%isp,dom(ib)%iep
                      do j=dom(ib)%jsp,dom(ib)%jep

                          if (i==dom(ib)%isp) then
                          uijk=1.0_dp/16.0_dp*(5.0_dp*dom(ib)%u(i-1,j,k)+15.0_dp*dom(ib)%u(i,j,k)- &
                                  5.0_dp*dom(ib)%u(i+1,j,k)+dom(ib)%u(i+2,j,k))
                          else if (i==dom(ib)%iep) then
                          uijk=1.0_dp/16.0_dp*(dom(ib)%u(i-3,j,k)-5.0_dp*dom(ib)%u(i-2,j,k)+ &
                             15.0_dp*dom(ib)%u(i-1,j,k)+5.0_dp*dom(ib)%u(i,j,k))
                          else
                          uijk=1.0_dp/16.0_dp*(-dom(ib)%u(i-2,j,k)+9.0_dp*dom(ib)%u(i-1,j,k)+ &
                               9.0_dp*dom(ib)%u(i,j,k)-dom(ib)%u(i+1,j,k))
                          end if

                          if (j==dom(ib)%jsp) then
                          vijk=1.0_dp/16.0_dp*(5.0_dp*dom(ib)%v(i,j-1,k)+15.0_dp*dom(ib)%v(i,j,k)- &
                                  5.0_dp*dom(ib)%v(i,j+1,k)+dom(ib)%v(i,j+2,k))
                          else if (j==dom(ib)%jep) then
                          vijk=1.0_dp/16.0_dp*(dom(ib)%v(i,j-3,k)-5.0_dp*dom(ib)%v(i,j-2,k)+ &
                             15.0_dp*dom(ib)%v(i,j-1,k)+5.0_dp*dom(ib)%v(i,j,k))
                          else
                          vijk=1.0_dp/16.0_dp*(-dom(ib)%v(i,j-2,k)+9.0_dp*dom(ib)%v(i,j-1,k)+ &
                               9.0_dp*dom(ib)%v(i,j,k)-dom(ib)%v(i,j+1,k))
                          end if

                          if (k==dom(ib)%ksp) then
                          wijk=1.0_dp/16.0_dp*(5.0_dp*dom(ib)%w(i,j,k-1)+15.0_dp*dom(ib)%w(i,j,k)- &
                                  5.0_dp*dom(ib)%w(i,j,k+1)+dom(ib)%w(i,j,k+2))
                          else if (k==dom(ib)%kep) then
                          wijk=1.0_dp/16.0_dp*(dom(ib)%w(i,j,k-3)-5.0_dp*dom(ib)%w(i,j,k-2)+ &
                             15.0_dp*dom(ib)%w(i,j,k-1)+5.0_dp*dom(ib)%w(i,j,k))
                          else
                          wijk=1.0_dp/16.0_dp*(-dom(ib)%w(i,j,k-2)+9.0_dp*dom(ib)%w(i,j,k-1)+ &
                               9.0_dp*dom(ib)%w(i,j,k)-dom(ib)%w(i,j,k+1))
                          end if
!
! Convection
!
                          h1 = uijk*dom(ib)%dphi_dx(i,j,k)
                          h2 = vijk*dom(ib)%dphi_dy(i,j,k)
                          h3 = wijk*dom(ib)%dphi_dz(i,j,k)

                          dom(ib)%phi_new(i,j,k) = dom(ib)%phi(i,j,k)-dt*(h1+h2+h3)

                      end do
                  end do
              end do

          end do

          call bound_lsm(16)    !(phi_new)
!
! 2nd RK step (phi_new --> pih_init)
!
          call dphi_a_v_3d(16)

          do ib=1,nbp

              do k=dom(ib)%ksp,dom(ib)%kep
                  do i=dom(ib)%isp,dom(ib)%iep
                      do j=dom(ib)%jsp,dom(ib)%jep

                          if (i==dom(ib)%isp) then
                          uijk=1.0_dp/16.0_dp*(5.0_dp*dom(ib)%u(i-1,j,k)+15.0_dp*dom(ib)%u(i,j,k)- &
                                  5.0_dp*dom(ib)%u(i+1,j,k)+dom(ib)%u(i+2,j,k))
                          else if (i==dom(ib)%iep) then
                          uijk=1.0_dp/16.0_dp*(dom(ib)%u(i-3,j,k)-5.0_dp*dom(ib)%u(i-2,j,k)+ &
                             15.0_dp*dom(ib)%u(i-1,j,k)+5.0_dp*dom(ib)%u(i,j,k))
                          else
                          uijk=1.0_dp/16.0_dp*(-dom(ib)%u(i-2,j,k)+9.0_dp*dom(ib)%u(i-1,j,k)+ &
                               9.0_dp*dom(ib)%u(i,j,k)-dom(ib)%u(i+1,j,k))
                          end if

                          if (j==dom(ib)%jsp) then
                          vijk=1.0_dp/16.0_dp*(5.0_dp*dom(ib)%v(i,j-1,k)+15.0_dp*dom(ib)%v(i,j,k)- &
                                  5.0_dp*dom(ib)%v(i,j+1,k)+dom(ib)%v(i,j+2,k))
                          else if (j==dom(ib)%jep) then
                          vijk=1.0_dp/16.0_dp*(dom(ib)%v(i,j-3,k)-5.0_dp*dom(ib)%v(i,j-2,k)+ &
                             15.0_dp*dom(ib)%v(i,j-1,k)+5.0_dp*dom(ib)%v(i,j,k))
                          else
                          vijk=1.0_dp/16.0_dp*(-dom(ib)%v(i,j-2,k)+9.0_dp*dom(ib)%v(i,j-1,k)+ &
                               9.0_dp*dom(ib)%v(i,j,k)-dom(ib)%v(i,j+1,k))
                          end if

                          if (k==dom(ib)%ksp) then
                          wijk=1.0_dp/16.0_dp*(5.0_dp*dom(ib)%w(i,j,k-1)+15.0_dp*dom(ib)%w(i,j,k)- &
                                  5.0_dp*dom(ib)%w(i,j,k+1)+dom(ib)%w(i,j,k+2))
                          else if (k==dom(ib)%kep) then
                          wijk=1.0_dp/16.0_dp*(dom(ib)%w(i,j,k-3)-5.0_dp*dom(ib)%w(i,j,k-2)+ &
                             15.0_dp*dom(ib)%w(i,j,k-1)+5.0_dp*dom(ib)%w(i,j,k))
                          else
                          wijk=1.0_dp/16.0_dp*(-dom(ib)%w(i,j,k-2)+9.0_dp*dom(ib)%w(i,j,k-1)+ &
                               9.0_dp*dom(ib)%w(i,j,k)-dom(ib)%w(i,j,k+1))
                          end if
!
! Convection
!
                          h1 = uijk*dom(ib)%dphi_dx(i,j,k)
                          h2 = vijk*dom(ib)%dphi_dy(i,j,k)
                          h3 = wijk*dom(ib)%dphi_dz(i,j,k)

                          dom(ib)%phi_init(i,j,k) = 0.75_dp*dom(ib)%phi(i,j,k) + &
                    0.25_dp*dom(ib)%phi_new(i,j,k) - 0.25_dp*dt*(h1+h2+h3)

                      end do
                  end do
              end do

          end do

          call bound_lsm(17)   !(phi_init)
!
! 3rd RK step (phi_init --> phi)
!
          call dphi_a_v_3d(17)

          do ib=1,nbp

              do k=dom(ib)%ksp,dom(ib)%kep
                  do i=dom(ib)%isp,dom(ib)%iep
                      do j=dom(ib)%jsp,dom(ib)%jep

                          if (i==dom(ib)%isp) then
                          uijk=1.0_dp/16.0_dp*(5.0_dp*dom(ib)%u(i-1,j,k)+15.0_dp*dom(ib)%u(i,j,k)- &
                                  5.0_dp*dom(ib)%u(i+1,j,k)+dom(ib)%u(i+2,j,k))
                          else if (i==dom(ib)%iep) then
                          uijk=1.0_dp/16.0_dp*(dom(ib)%u(i-3,j,k)-5.0_dp*dom(ib)%u(i-2,j,k)+ &
                             15.0_dp*dom(ib)%u(i-1,j,k)+5.0_dp*dom(ib)%u(i,j,k))
                          else
                          uijk=1.0_dp/16.0_dp*(-dom(ib)%u(i-2,j,k)+9.0_dp*dom(ib)%u(i-1,j,k)+ &
                               9.0_dp*dom(ib)%u(i,j,k)-dom(ib)%u(i+1,j,k))
                          end if

                          if (j==dom(ib)%jsp) then
                          vijk=1.0_dp/16.0_dp*(5.0_dp*dom(ib)%v(i,j-1,k)+15.0_dp*dom(ib)%v(i,j,k)- &
                                  5.0_dp*dom(ib)%v(i,j+1,k)+dom(ib)%v(i,j+2,k))
                          else if (j==dom(ib)%jep) then
                          vijk=1.0_dp/16.0_dp*(dom(ib)%v(i,j-3,k)-5.0_dp*dom(ib)%v(i,j-2,k)+ &
                             15.0_dp*dom(ib)%v(i,j-1,k)+5.0_dp*dom(ib)%v(i,j,k))
                          else
                          vijk=1.0_dp/16.0_dp*(-dom(ib)%v(i,j-2,k)+9.0_dp*dom(ib)%v(i,j-1,k)+ &
                               9.0_dp*dom(ib)%v(i,j,k)-dom(ib)%v(i,j+1,k))
                          end if

                          if (k==dom(ib)%ksp) then
                          wijk=1.0_dp/16.0_dp*(5.0_dp*dom(ib)%w(i,j,k-1)+15.0_dp*dom(ib)%w(i,j,k)- &
                                  5.0_dp*dom(ib)%w(i,j,k+1)+dom(ib)%w(i,j,k+2))
                          else if (k==dom(ib)%kep) then
                          wijk=1.0_dp/16.0_dp*(dom(ib)%w(i,j,k-3)-5.0_dp*dom(ib)%w(i,j,k-2)+ &
                             15.0_dp*dom(ib)%w(i,j,k-1)+5.0_dp*dom(ib)%w(i,j,k))
                          else
                          wijk=1.0_dp/16.0_dp*(-dom(ib)%w(i,j,k-2)+9.0_dp*dom(ib)%w(i,j,k-1)+ &
                               9.0_dp*dom(ib)%w(i,j,k)-dom(ib)%w(i,j,k+1))
                          end if
!
! Convection
!
                          h1 = uijk*dom(ib)%dphi_dx(i,j,k)
                          h2 = vijk*dom(ib)%dphi_dy(i,j,k)
                          h3 = wijk*dom(ib)%dphi_dz(i,j,k)

                          dom(ib)%phi(i,j,k) = 1.0_dp/3.0_dp*dom(ib)%phi(i,j,k) + &
                    2.0_dp/3.0_dp*dom(ib)%phi_init(i,j,k) - 2.0_dp/3.0_dp*dt*(h1+h2+h3)   ! We now have updated phi field

!
! Hold level to be held constant at inflow and outflow (if required - may help with stability in inflow-outflow sims)
!
                          if (trim(keyword)=="channel" .and. lends) then
                          if (dom(ib)%iprev<0) then
                          if ((i>=dom(ib)%isu).and.(i<=dom(ib)%isu+5)) then
                          if (dom(ib)%zc(k)<length)   then
                          dom(ib)%phi(i,j,k) = 1.0_dp*abs(dom(ib)%zc(k)-length)
                          else if (dom(ib)%zc(k)>length)   then
                          dom(ib)%phi(i,j,k) = -1.0_dp*abs(dom(ib)%zc(k)-length)
                          else if (dom(ib)%zc(k)==length)  then
                          dom(ib)%phi(i,j,k) = 0.0_dp
                          end if
                          end if
                          end if
!          length=0.0354
                          if  (dom(ib)%inext<0) then
                          if ((i<=dom(ib)%ieu).and.(i>=dom(ib)%ieu-5)) then
                          if (dom(ib)%zc(k)<length)   then
                          dom(ib)%phi(i,j,k) = 1.0_dp*abs(dom(ib)%zc(k)-length)
                          else if (dom(ib)%zc(k)>length)   then
                          dom(ib)%phi(i,j,k) = -1.0_dp*abs(dom(ib)%zc(k)-length)
                          else if (dom(ib)%zc(k)==length)  then
                          dom(ib)%phi(i,j,k) = 0.0_dp
                          end if
                          end if
                          end if
                          end if

                      end do
                  end do
              end do

          end do

          call bound_lsm(14)  !(phi)
!
! Reinitialise phi field to a signed distance function
!
          if (reinit) then
          call tvd_rk_reinit     !(phi)
          end if

          return
      end subroutine tvd_rk_3step
!#######################################################################
      subroutine dphi_a_v_3d(op)     !(fi)
!#######################################################################
          use vars
          use module_LSM
          use multiflow3d_mpi
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none

          integer :: i,j,k,ifi,op,ib
          real(dp) :: uijk,vijk,wijk,phi_minus_ip12j
          real(dp) :: phi_plus_ip12j,phi_minus_ijp12
          real(dp) :: phi_plus_ijp12,superbee

          call hj_weno_dxplus_3d(op)
          call hj_weno_dyplus_3d(op)
          call hj_weno_dzplus_3d(op)
          call hj_weno_dxminus_3d(op)
          call hj_weno_dyminus_3d(op)
          call hj_weno_dzminus_3d(op)

          do ib=1,nbp

              do k=dom(ib)%ksp,dom(ib)%kep
                  do i=dom(ib)%isp,dom(ib)%iep
                      do j=dom(ib)%jsp,dom(ib)%jep

                          if (i==dom(ib)%isp) then
                          uijk=1.0_dp/16.0_dp*(5.0_dp*dom(ib)%u(i-1,j,k)+15.0_dp*dom(ib)%u(i,j,k)- &
                                   5.0_dp*dom(ib)%u(i+1,j,k)+dom(ib)%u(i+2,j,k))
                          else if (i==dom(ib)%iep) then
                          uijk=1.0_dp/16.0_dp*(dom(ib)%u(i-3,j,k)-5.0_dp*dom(ib)%u(i-2,j,k)+ &
                              15.0_dp*dom(ib)%u(i-1,j,k)+5.0_dp*dom(ib)%u(i,j,k))
                          else
                          uijk=1.0_dp/16.0_dp*(-dom(ib)%u(i-2,j,k)+9.0_dp*dom(ib)%u(i-1,j,k)+ &
                                9.0_dp*dom(ib)%u(i,j,k)-dom(ib)%u(i+1,j,k))
                          end if

                          if (j==dom(ib)%jsp) then
                          vijk=1.0_dp/16.0_dp*(5.0_dp*dom(ib)%v(i,j-1,k)+15.0_dp*dom(ib)%v(i,j,k)- &
                                   5.0_dp*dom(ib)%v(i,j+1,k)+dom(ib)%v(i,j+2,k))
                          else if (j==dom(ib)%jep) then
                          vijk=1.0_dp/16.0_dp*(dom(ib)%v(i,j-3,k)-5.0_dp*dom(ib)%v(i,j-2,k)+ &
                              15.0_dp*dom(ib)%v(i,j-1,k)+5.0_dp*dom(ib)%v(i,j,k))
                          else
                          vijk=1.0_dp/16.0_dp*(-dom(ib)%v(i,j-2,k)+9.0_dp*dom(ib)%v(i,j-1,k)+ &
                                9.0_dp*dom(ib)%v(i,j,k)-dom(ib)%v(i,j+1,k))
                          end if

                          if (k==dom(ib)%ksp) then
                          wijk=1.0_dp/16.0_dp*(5.0_dp*dom(ib)%w(i,j,k-1)+15.0_dp*dom(ib)%w(i,j,k)- &
                                   5.0_dp*dom(ib)%w(i,j,k+1)+dom(ib)%w(i,j,k+2))
                          else if (k==dom(ib)%kep) then
                          wijk=1.0_dp/16.0_dp*(dom(ib)%w(i,j,k-3)-5.0_dp*dom(ib)%w(i,j,k-2)+ &
                              15.0_dp*dom(ib)%w(i,j,k-1)+5.0_dp*dom(ib)%w(i,j,k))
                          else
                          wijk=1.0_dp/16.0_dp*(-dom(ib)%w(i,j,k-2)+9.0_dp*dom(ib)%w(i,j,k-1)+ &
                                9.0_dp*dom(ib)%w(i,j,k)-dom(ib)%w(i,j,k+1))
                          end if

                          if (uijk>0.0_dp) then
                            dom(ib)%dphi_dx(i,j,k) = dom(ib)%dphi_dxminus(i,j,k)
                          end if
                          if (uijk<0.0_dp) then
                            dom(ib)%dphi_dx(i,j,k) = dom(ib)%dphi_dxplus(i,j,k)
                          end if
                          if (uijk==0.0_dp) then
                            dom(ib)%dphi_dx(i,j,k) = 0.0_dp
                          end if

                          if (vijk>0.0_dp) then
                            dom(ib)%dphi_dy(i,j,k) = dom(ib)%dphi_dyminus(i,j,k)
                          end if
                          if (vijk<0.0_dp) then
                            dom(ib)%dphi_dy(i,j,k) = dom(ib)%dphi_dyplus(i,j,k)
                          end if
                          if (vijk==0.0_dp) then
                            dom(ib)%dphi_dy(i,j,k) = 0.0_dp
                          end if

                          if (wijk>0.0_dp) then
                            dom(ib)%dphi_dz(i,j,k) = dom(ib)%dphi_dzminus(i,j,k)
                          end if
                          if (wijk<0.0_dp) then
                            dom(ib)%dphi_dz(i,j,k) = dom(ib)%dphi_dzplus(i,j,k)
                          end if
                          if (wijk==0.0_dp) then
                            dom(ib)%dphi_dz(i,j,k) = 0.0_dp
                          end if

                      end do
                  end do
              end do

          end do

          return
      end subroutine dphi_a_v_3d
!######################################################################
      subroutine tvd_rk_reinit
!######################################################################
!**********************************************************************
! 3 step runge kutta routine which returns re-initialised phi field
!**********************************************************************
          use vars
          use module_LSM
          use multidata
          use multiflow3d_mpi
          use, intrinsic :: iso_fortran_env, only: dp => real64

          implicit none
          integer :: i,j,k,it,ib
          logical :: bool
          real(dp) :: max_abs,abs_dphi,abs_phidiff,max_phidiff
          real(dp) :: local_max_abs,local_max_phidiff,dt_reinit

          bool=.false.
          it=0

          do while (.not.bool)
!
! 1st step
!
              call dphi_for_reinit(14)      ! (phi)

              do ib=1,nbp

                  do k=dom(ib)%ksp,dom(ib)%kep
                      do i=dom(ib)%isp,dom(ib)%iep
                          do j=dom(ib)%jsp,dom(ib)%jep

                              abs_dphi = sqrt(dom(ib)%dphi_dx(i,j,k)**2+ &
                        dom(ib)%dphi_dy(i,j,k)**2+dom(ib)%dphi_dz(i,j,k)**2)

!     dom(ib)%abs_dphi_check(i,j,k)=abs_dphi

                              dom(ib)%s_phi0(i,j,k) = dom(ib)%phi(i,j,k)/ &
                            sqrt(dom(ib)%phi(i,j,k)**2+(abs_dphi**2)*dom(ib)%dx**2)

                              dom(ib)%phi_reinit(i,j,k) =dom(ib)%phi(i,j,k) + &
                        dt_reinit*(dom(ib)%s_phi0(i,j,k)- &
                         dom(ib)%s_phi0(i,j,k)*abs_dphi)

                          end do
                      end do
                  end do

              end do

              call bound_lsm(15)  ! (phi_reinit)
!
! 2nd step
!
              call dphi_for_reinit(15)  ! (phi_reinit)

              do ib=1,nbp

                  do k=dom(ib)%ksp,dom(ib)%kep
                      do i=dom(ib)%isp,dom(ib)%iep
                          do j=dom(ib)%jsp,dom(ib)%jep

                              abs_dphi = sqrt(dom(ib)%dphi_dx(i,j,k)**2+ &
                        dom(ib)%dphi_dy(i,j,k)**2+dom(ib)%dphi_dz(i,j,k)**2)

                              dom(ib)%s_phi0(i,j,k) = dom(ib)%phi_reinit(i,j,k)/ &
                        sqrt(dom(ib)%phi_reinit(i,j,k)**2+(abs_dphi**2)*dom(ib)%dx**2)

                              dom(ib)%phi_reinit(i,j,k) = 0.75_dp*dom(ib)%phi(i,j,k)+ &
                        0.25_dp*dom(ib)%phi_reinit(i,j,k)+ &
                        0.25_dp*dt_reinit*(dom(ib)%s_phi0(i,j,k)- &
                        dom(ib)%s_phi0(i,j,k)*abs_dphi)

                          end do
                      end do
                  end do

              end do

              call bound_lsm(15)  ! (phi_reinit)
!
! 3rd step
!
              call dphi_for_reinit(15)  ! (phi_reinit)

              do ib=1,nbp

                  do k=dom(ib)%ksp,dom(ib)%kep
                      do i=dom(ib)%isp,dom(ib)%iep
                          do j=dom(ib)%jsp,dom(ib)%jep

                              abs_dphi = sqrt(dom(ib)%dphi_dx(i,j,k)**2+ &
                        dom(ib)%dphi_dy(i,j,k)**2+dom(ib)%dphi_dz(i,j,k)**2)

                              dom(ib)%s_phi0(i,j,k) = dom(ib)%phi_reinit(i,j,k)/ &
                        sqrt(dom(ib)%phi_reinit(i,j,k)**2+(abs_dphi**2)*dom(ib)%dx**2)

                              dom(ib)%phi_reinit(i,j,k) =1.0_dp/3.0_dp*dom(ib)%phi(i,j,k)+ &
                        2.0_dp/3.0_dp*dom(ib)%phi_reinit(i,j,k)+ &
                        2.0_dp/3.0_dp*dt_reinit*(dom(ib)%s_phi0(i,j,k)- &
                        dom(ib)%s_phi0(i,j,k)*abs_dphi)

                          end do
                      end do
                  end do

              end do

              call bound_lsm(15)  ! (phi_reinit)

!
! Check convergence
!
              call dphi_for_reinit(15)  ! (phi_reinit)

              max_abs = 0.0_dp
              max_phidiff = 0.0_dp

              do ib=1,nbp

                  do k=dom(ib)%ksp,dom(ib)%kep
                      do i=dom(ib)%isp,dom(ib)%iep
                          do j=dom(ib)%jsp,dom(ib)%jep

                              abs_dphi = sqrt(dom(ib)%dphi_dx(i,j,k)**2+ &
                        dom(ib)%dphi_dy(i,j,k)**2+dom(ib)%dphi_dz(i,j,k)**2)

                              max_abs  = max(max_abs,abs_dphi)
                              abs_phidiff=abs(dom(ib)%phi(i,j,k)-dom(ib)%phi_reinit(i,j,k))
                              max_phidiff = max(max_phidiff,abs_phidiff)

                          end do
                      end do
                  end do

              end do

              local_max_abs = max_abs

              call mpi_allreduce(local_max_abs,max_abs,1,mpi_flt,mpi_max, &
        mpi_comm_world,ierr)

              local_max_phidiff = max_phidiff

              call mpi_allreduce(local_max_phidiff,max_phidiff,1,mpi_flt, &
        mpi_max,mpi_comm_world,ierr)

              if ((max_phidiff<reldif_lsm).and.(it>=1)) then
              bool=.true.
              else
              if (it>=ntime_reinit) then
              bool=.true.
              else
              it=it+1
              end if
              end if

              do ib=1,nbp
                  do k=dom(ib)%ksp,dom(ib)%kep
                      do i=dom(ib)%isp,dom(ib)%iep
                          do j=dom(ib)%jsp,dom(ib)%jep
                              dom(ib)%phi(i,j,k) = dom(ib)%phi_reinit(i,j,k)
                          end do
                      end do
                  end do
              end do

              call bound_lsm(14)  ! (phi)

          end do

          if (myrank==0) then
          write(*,*) "norm v (reinit)", max_abs, "needed steps", it
          write(numfile3,"(i8,f18.8,i8,2f18.8)") ntime,max_abs,it,ctime, &
    dt
          end if

          return
      end subroutine tvd_rk_reinit
!######################################################################
      subroutine dphi_for_reinit(op)
!######################################################################
          use vars
          use module_LSM
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64

          implicit none

          integer :: i,j,k,ifi,op,ib
          real(dp) :: s_phi012
          real(dp), pointer, dimension(:,:,:)::fi
          real(dp) :: lsv,lssig,xm,xp,ym,yp,zm,zp

          call hj_weno_dxplus_3d(op)
          call hj_weno_dyplus_3d(op)
          call hj_weno_dzplus_3d(op)
          call hj_weno_dxminus_3d(op)
          call hj_weno_dyminus_3d(op)
          call hj_weno_dzminus_3d(op)

          do ib=1,nbp

              select case (op)
                case (14)  ! (phi)
                  fi => dom(ib)%phi
                case (15)  ! (phi_reinit)
                  fi => dom(ib)%phi_reinit
              end select

              do k=dom(ib)%ksp,dom(ib)%kep
                  do i=dom(ib)%isp,dom(ib)%iep
                      do j=dom(ib)%jsp,dom(ib)%jep

                          lsv = fi(i,j,k)
                          lssig = dom(ib)%s_phi0(i,j,k)  ! s_phi0=phi initially

                          xm=dom(ib)%dphi_dxminus(i,j,k)
                          xp=dom(ib)%dphi_dxplus(i,j,k)
                          ym=dom(ib)%dphi_dyminus(i,j,k)
                          yp=dom(ib)%dphi_dyplus(i,j,k)
                          zm=dom(ib)%dphi_dzminus(i,j,k)
                          zp=dom(ib)%dphi_dzplus(i,j,k)

                          if ((xm*lssig>0.0_dp).and.(xp*lssig>-xm*lssig)) then
                          dom(ib)%dphi_dx(i,j,k)=dom(ib)%dphi_dxminus(i,j,k)
                          end if

                          if((xp*lssig<0.0_dp).and.(xm*lssig<-xp*lssig)) then
                          dom(ib)%dphi_dx(i,j,k)=dom(ib)%dphi_dxplus(i,j,k)
                          end if

                          if ((xp*lssig>0.0_dp).and.(xm*lssig<0.0_dp)) then
                          dom(ib)%dphi_dx(i,j,k)=0.0_dp
                          end if

                          if ((ym*lssig>0.0_dp).and.(yp*lssig>-ym*lssig)) then
                          dom(ib)%dphi_dy(i,j,k)=dom(ib)%dphi_dyminus(i,j,k)
                          end if

                          if ((yp*lssig<0.0_dp).and.(ym*lssig<-yp*lssig)) then
                          dom(ib)%dphi_dy(i,j,k)=dom(ib)%dphi_dyplus(i,j,k)
                          end if

                          if ((yp*lssig>0.0_dp).and.(ym*lssig<0.0_dp)) then
                          dom(ib)%dphi_dy(i,j,k) = 0.0_dp
                          end if

                          if ((zm*lssig>0.0_dp).and.(zp*lssig>-zm*lssig)) then
                          dom(ib)%dphi_dz(i,j,k)=dom(ib)%dphi_dzminus(i,j,k)
                          end if

                          if((zp*lssig<0.0_dp).and.(zm*lssig<-zp*lssig)) then
                          dom(ib)%dphi_dz(i,j,k)=dom(ib)%dphi_dzplus(i,j,k)
                          end if

                          if ((zp*lssig>0.0_dp).and.(zm*lssig<0.0_dp)) then
                          dom(ib)%dphi_dz(i,j,k)=0.0_dp
                          end if

                      end do
                  end do
              end do

          end do

          return
      end subroutine dphi_for_reinit
!######################################################################
      subroutine hj_weno_dxplus_3d(op)
!######################################################################
          use vars
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64

          implicit none

          integer :: i,j,k,l,k_star,npp,nq,nr,op,ib
          real(dp), pointer, dimension(:,:,:)::fi
          real(dp) :: q1,q2,q3,c_star,v1,v2,v3,v4,v5,v11,v22
          real(dp) :: v33,v44,v55,s1,s2,s3,a1,a2,a3,w1,w2,w3
          real(dp) :: e

          do ib=1,nbp

              select case (op)
                case (1)
                  fi => dom(ib)%u
                  npp = dom(ib)%niul
                  nq = dom(ib)%njul
                  nr = dom(ib)%nkul
                case (2)
                  fi => dom(ib)%v
                  npp = dom(ib)%nivl
                  nq = dom(ib)%njvl
                  nr = dom(ib)%nkvl
                case (3)
                  fi => dom(ib)%w
                  npp = dom(ib)%niwl
                  nq = dom(ib)%njwl
                  nr = dom(ib)%nkwl
                case (14)  ! (phi)
                  fi => dom(ib)%phi
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
                case (15)  ! (phi_reinit)
                  fi => dom(ib)%phi_reinit
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
                case (16)  ! (phi_new)
                  fi => dom(ib)%phi_new
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
                case (17)  ! (phi_init)
                  fi => dom(ib)%phi_init
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
              end select

              if (dom(ib)%iprev<0 .and. dom(ib)%bc_west/=5) then
              do k=1,nr
                  do j=1,nq
                      fi(2,j,k)    = fi(3,j,k)
                      fi(1,j,k)    = fi(2,j,k)
                  end do
              end do
              else if (dom(ib)%inext<0 .and. dom(ib)%bc_east/=5) then
              do k=1,nr
                  do j=1,nq
                      fi(npp-1,j,k) = fi(npp-2,j,k)
                      fi(npp,j,k)   = fi(npp-1,j,k)
                  end do
              end do
              end if
!
! Compute 1st order spatial derivative
!
              do k=1,nr
                  do i=1,npp-1
                      do j=1,nq
                          dom(ib)%d1(i,j,k) = (fi(i+1,j,k) - fi(i,j,k))/dom(ib)%dx
                      end do
                  end do
              end do
!
! Compute 5th order derivative (+ve direction)
!
              dom(ib)%dphi_dxplus = 0.0_dp

              do k=1,nr
                  do i=1,npp-6
                      do j=1,nq

                          l=i

                          v1 = dom(ib)%d1(l+5,j,k)
                          v2 = dom(ib)%d1(l+4,j,k)
                          v3 = dom(ib)%d1(l+3,j,k)
                          v4 = dom(ib)%d1(l+2,j,k)
                          v5 = dom(ib)%d1(l+1,j,k)

                          v11 = v1**2
                          v22 = v2**2
                          v33 = v3**2
                          v44 = v4**2
                          v55 = v5**2

                          e = (1.0e-6_dp) * max(v11,v22,v33,v44,v55) + 1E-99_dp

                          s1 = (13.0_dp/12.0_dp)*(v1-2.0_dp*v2+v3)**2 + 0.25_dp*(v1-4.0_dp*v2+3.0_dp*v3)**2
                          s2 = (13.0_dp/12.0_dp)*(v2-2.0_dp*v3+v4)**2 + 0.25_dp*(v2-v4)**2
                          s3 = (13.0_dp/12.0_dp)*(v3-2.0_dp*v4+v5)**2 + 0.25_dp*(3.0_dp*v3-4.0_dp*v4+v5)**2

                          a1 = 0.1_dp/(e+s1)**2
                          a2 = 0.6_dp/(e+s2)**2
                          a3 = 0.3_dp/(e+s3)**2

                          w1 = a1/(a1+a2+a3)
                          w2 = a2/(a1+a2+a3)
                          w3 = a3/(a1+a2+a3)

                          dom(ib)%dphi_dxplus(i+3,j,k) = w1*(v1*(1.0_dp/3.0_dp)-v2*(7.0_dp/6.0_dp)+ &
                    v3*(11.0_dp/6.0_dp)) + w2*(-v2*(1.0_dp/6.0_dp)+v3*(5.0_dp/6.0_dp)+v4*(1.0_dp/3.0_dp))+ &
                    w3*(v3*(1.0_dp/3.0_dp)+v4*(5.0_dp/6.0_dp)-v5*(1.0_dp/6.0_dp))

                      end do
                  end do
              end do

          end do

          return
      end subroutine hj_weno_dxplus_3d
!######################################################################
      subroutine hj_weno_dxminus_3d(op)
!######################################################################
          use vars
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64

          implicit none

          integer :: i,j,k,l,k_star,npp,nq,nr,op,ib
          real(dp), pointer, dimension(:,:,:)::fi
          real(dp) :: q1,q2,q3,c_star,v1,v2,v3,v4,v5,v11,v22
          real(dp) :: v33,v44,v55,s1,s2,s3,a1,a2,a3,w1,w2,w3
          real(dp) :: e

          do ib=1,nbp

              select case (op)
                case (1)
                  fi => dom(ib)%u
                  npp = dom(ib)%niul
                  nq = dom(ib)%njul
                  nr = dom(ib)%nkul
                case (2)
                  fi => dom(ib)%v
                  npp = dom(ib)%nivl
                  nq = dom(ib)%njvl
                  nr = dom(ib)%nkvl
                case (3)
                  fi => dom(ib)%w
                  npp = dom(ib)%niwl
                  nq = dom(ib)%njwl
                  nr = dom(ib)%nkwl
                case (14)  ! (phi)
                  fi => dom(ib)%phi
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
                case (15)  ! (phi_reinit)
                  fi => dom(ib)%phi_reinit
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
                case (16)  ! (phi_new)
                  fi => dom(ib)%phi_new
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
                case (17)  ! (phi_init)
                  fi => dom(ib)%phi_init
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
              end select

              if (dom(ib)%iprev<0 .and. dom(ib)%bc_west/=5) then
              do k=1,nr
                  do j=1,nq
                      fi(2,j,k)    = fi(3,j,k)
                      fi(1,j,k)    = fi(2,j,k)
                  end do
              end do
              else if (dom(ib)%inext<0 .and. dom(ib)%bc_east/=5) then
              do k=1,nr
                  do j=1,nq
                      fi(npp-1,j,k) = fi(npp-2,j,k)
                      fi(npp,j,k)   = fi(npp-1,j,k)
                  end do
              end do
              end if
!
! Compute 1st order spatial derivative
!
              do k=1,nr
                  do i=1,npp-1
                      do j=1,nq
                          dom(ib)%d1(i,j,k) = (fi(i+1,j,k) - fi(i,j,k))/dom(ib)%dx
                      end do
                  end do
              end do
!
! Compute 5th order derivative (-ve direction)
!
              dom(ib)%dphi_dxminus = 0.0_dp

              do k=1,nr
                  do i=1,npp-6
                      do j=1,nq

                          l=i-1

                          v1 = dom(ib)%d1(l+1,j,k)
                          v2 = dom(ib)%d1(l+2,j,k)
                          v3 = dom(ib)%d1(l+3,j,k)
                          v4 = dom(ib)%d1(l+4,j,k)
                          v5 = dom(ib)%d1(l+5,j,k)

                          v11 = v1**2
                          v22 = v2**2
                          v33 = v3**2
                          v44 = v4**2
                          v55 = v5**2

                          e = (1.0e-6_dp) * max(v11,v22,v33,v44,v55) + 1E-99_dp

                          s1 = (13.0_dp/12.0_dp)*(v1-2.0_dp*v2+v3)**2 + 0.25_dp*(v1-4.0_dp*v2+3.0_dp*v3)**2
                          s2 = (13.0_dp/12.0_dp)*(v2-2.0_dp*v3+v4)**2 + 0.25_dp*(v2-v4)**2
                          s3 = (13.0_dp/12.0_dp)*(v3-2.0_dp*v4+v5)**2 + 0.25_dp*(3.0_dp*v3-4.0_dp*v4+v5)**2

                          a1 = 0.1_dp/(e+s1)**2
                          a2 = 0.6_dp/(e+s2)**2
                          a3 = 0.3_dp/(e+s3)**2

                          w1 = a1/(a1+a2+a3)
                          w2 = a2/(a1+a2+a3)
                          w3 = a3/(a1+a2+a3)

                          dom(ib)%dphi_dxminus(i+3,j,k) = w1*(v1*(1.0_dp/3.0_dp)-v2*(7.0_dp/6.0_dp)+ &
                    v3*(11.0_dp/6.0_dp)) + w2*(-v2*(1.0_dp/6.0_dp)+v3*(5.0_dp/6.0_dp)+v4*(1.0_dp/3.0_dp))+ &
                    w3*(v3*(1.0_dp/3.0_dp)+v4*(5.0_dp/6.0_dp)-v5*(1.0_dp/6.0_dp))

                      end do
                  end do
              end do

          end do

          return
      end subroutine hj_weno_dxminus_3d
!######################################################################
      subroutine hj_weno_dyplus_3d(op)
!######################################################################
          use vars
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64

          implicit none

          integer :: i,j,k,l,k_star,npp,nq,nr,op,ib
          real(dp), pointer, dimension(:,:,:)::fi
          real(dp) :: q1,q2,q3,c_star,v1,v2,v3,v4,v5,v11,v22
          real(dp) :: v33,v44,v55,s1,s2,s3,a1,a2,a3,w1,w2,w3
          real(dp) :: e

          do ib=1,nbp

              select case (op)
                case (1)
                  fi => dom(ib)%u
                  npp = dom(ib)%niul
                  nq = dom(ib)%njul
                  nr = dom(ib)%nkul
                case (2)
                  fi => dom(ib)%v
                  npp = dom(ib)%nivl
                  nq = dom(ib)%njvl
                  nr = dom(ib)%nkvl
                case (3)
                  fi => dom(ib)%w
                  npp = dom(ib)%niwl
                  nq = dom(ib)%njwl
                  nr = dom(ib)%nkwl
                case (14)  ! (phi)
                  fi => dom(ib)%phi
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
                case (15)  ! (phi_reinit)
                  fi => dom(ib)%phi_reinit
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
                case (16)  ! (phi_new)
                  fi => dom(ib)%phi_new
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
                case (17)  ! (phi_init)
                  fi => dom(ib)%phi_init
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
              end select

              if (dom(ib)%jprev<0 .and. dom(ib)%bc_south/=5) then
              do k=1,nr
                  do i=1,npp
                      fi(i,2,k)    = fi(i,3,k)
                      fi(i,1,k)    = fi(i,2,k)
                  end do
              end do
              else if (dom(ib)%jnext<0 .and. dom(ib)%bc_north/=5) then
              do k=1,nr
                  do i=1,npp
                      fi(i,nq-1,k) = fi(i,nq-2,k)
                      fi(i,nq,k)   = fi(i,nq-1,k)
                  end do
              end do
              end if
!
! Compute 1st order spatial derivative
!
              do k=1,nr
                  do i=1,npp
                      do j=1,nq-1
                          dom(ib)%d1(i,j,k) = (fi(i,j+1,k) - fi(i,j,k))/dom(ib)%dy
                      end do
                  end do
              end do
!
! Compute 5th order derivative (+ve direction)
!
              dom(ib)%dphi_dyplus = 0.0_dp

              do k=1,nr
                  do i=1,npp
                      do j=1,nq-6
                          l=j

                          v1 = dom(ib)%d1(i,l+5,k)
                          v2 = dom(ib)%d1(i,l+4,k)
                          v3 = dom(ib)%d1(i,l+3,k)
                          v4 = dom(ib)%d1(i,l+2,k)
                          v5 = dom(ib)%d1(i,l+1,k)

                          v11 = v1**2
                          v22 = v2**2
                          v33 = v3**2
                          v44 = v4**2
                          v55 = v5**2

                          e = (1.0e-6_dp) * max(v11,v22,v33,v44,v55) + 1E-99_dp

                          s1 = (13.0_dp/12.0_dp)*(v1-2.0_dp*v2+v3)**2 + 0.25_dp*(v1-4.0_dp*v2+3.0_dp*v3)**2
                          s2 = (13.0_dp/12.0_dp)*(v2-2.0_dp*v3+v4)**2 + 0.25_dp*(v2-v4)**2
                          s3 = (13.0_dp/12.0_dp)*(v3-2.0_dp*v4+v5)**2 + 0.25_dp*(3.0_dp*v3-4.0_dp*v4+v5)**2

                          a1 = 0.1_dp/(e+s1)**2
                          a2 = 0.6_dp/(e+s2)**2
                          a3 = 0.3_dp/(e+s3)**2

                          w1 = a1/(a1+a2+a3)
                          w2 = a2/(a1+a2+a3)
                          w3 = a3/(a1+a2+a3)

                          dom(ib)%dphi_dyplus(i,j+3,k) = w1*(v1*(1.0_dp/3.0_dp)-v2*(7.0_dp/6.0_dp)+ &
                    v3*(11.0_dp/6.0_dp)) + w2*(-v2*(1.0_dp/6.0_dp)+v3*(5.0_dp/6.0_dp)+v4*(1.0_dp/3.0_dp))+ &
                    w3*(v3*(1.0_dp/3.0_dp)+v4*(5.0_dp/6.0_dp)-v5*(1.0_dp/6.0_dp))

                      end do
                  end do
              end do

          end do

          return
      end subroutine hj_weno_dyplus_3d
!######################################################################
      subroutine hj_weno_dyminus_3d(op)
!######################################################################
          use vars
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64

          implicit none

          integer :: i,j,k,l,k_star,npp,nq,nr,op,ib
          real(dp), pointer, dimension(:,:,:)::fi
          real(dp) :: q1,q2,q3,c_star,v1,v2,v3,v4,v5,v11,v22
          real(dp) :: v33,v44,v55,s1,s2,s3,a1,a2,a3,w1,w2,w3
          real(dp) :: e

          do ib=1,nbp

              select case (op)
                case (1)
                  fi => dom(ib)%u
                  npp = dom(ib)%niul
                  nq = dom(ib)%njul
                  nr = dom(ib)%nkul
                case (2)
                  fi => dom(ib)%v
                  npp = dom(ib)%nivl
                  nq = dom(ib)%njvl
                  nr = dom(ib)%nkvl
                case (3)
                  fi => dom(ib)%w
                  npp = dom(ib)%niwl
                  nq = dom(ib)%njwl
                  nr = dom(ib)%nkwl
                case (14)  ! (phi)
                  fi => dom(ib)%phi
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
                case (15)  ! (phi_reinit)
                  fi => dom(ib)%phi_reinit
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
                case (16)  ! (phi_new)
                  fi => dom(ib)%phi_new
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
                case (17)  ! (phi_init)
                  fi => dom(ib)%phi_init
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
              end select

              if (dom(ib)%jprev<0 .and. dom(ib)%bc_south/=5) then
              do k=1,nr
                  do i=1,npp
                      fi(i,2,k)    = fi(i,3,k)
                      fi(i,1,k)    = fi(i,2,k)
                  end do
              end do
              else if (dom(ib)%jnext<0 .and. dom(ib)%bc_north/=5) then
              do k=1,nr
                  do i=1,npp
                      fi(i,nq-1,k) = fi(i,nq-2,k)
                      fi(i,nq,k)   = fi(i,nq-1,k)
                  end do
              end do
              end if
!
! Compute 1st order spatial derivative
!
              do k=1,nr
                  do i=1,npp
                      do j=1,nq-1
                          dom(ib)%d1(i,j,k) = (fi(i,j+1,k) - fi(i,j,k))/dom(ib)%dy
                      end do
                  end do
              end do
!
! Compute 5th order derivative (-ve direction)
!
              dom(ib)%dphi_dyminus = 0.0_dp

              do k=1,nr
                  do i=1,npp
                      do j=1,nq-6

                          l=j-1

                          v1 = dom(ib)%d1(i,l+1,k)
                          v2 = dom(ib)%d1(i,l+2,k)
                          v3 = dom(ib)%d1(i,l+3,k)
                          v4 = dom(ib)%d1(i,l+4,k)
                          v5 = dom(ib)%d1(i,l+5,k)

                          v11 = v1**2
                          v22 = v2**2
                          v33 = v3**2
                          v44 = v4**2
                          v55 = v5**2

                          e = (1.0e-6_dp) * max(v11,v22,v33,v44,v55) + 1E-99_dp

                          s1 = (13.0_dp/12.0_dp)*(v1-2.0_dp*v2+v3)**2 + 0.25_dp*(v1-4.0_dp*v2+3.0_dp*v3)**2
                          s2 = (13.0_dp/12.0_dp)*(v2-2.0_dp*v3+v4)**2 + 0.25_dp*(v2-v4)**2
                          s3 = (13.0_dp/12.0_dp)*(v3-2.0_dp*v4+v5)**2 + 0.25_dp*(3.0_dp*v3-4.0_dp*v4+v5)**2

                          a1 = 0.1_dp/(e+s1)**2
                          a2 = 0.6_dp/(e+s2)**2
                          a3 = 0.3_dp/(e+s3)**2

                          w1 = a1/(a1+a2+a3)
                          w2 = a2/(a1+a2+a3)
                          w3 = a3/(a1+a2+a3)

                          dom(ib)%dphi_dyminus(i,j+3,k) = w1*(v1*(1.0_dp/3.0_dp)-v2*(7.0_dp/6.0_dp)+ &
                    v3*(11.0_dp/6.0_dp)) + w2*(-v2*(1.0_dp/6.0_dp)+v3*(5.0_dp/6.0_dp)+v4*(1.0_dp/3.0_dp))+ &
                    w3*(v3*(1.0_dp/3.0_dp)+v4*(5.0_dp/6.0_dp)-v5*(1.0_dp/6.0_dp))

                      end do
                  end do
              end do

          end do

          return
      end subroutine hj_weno_dyminus_3d
!######################################################################
      subroutine hj_weno_dzplus_3d(op)
!######################################################################
          use vars
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64

          implicit none

          integer :: i,j,k,l,k_star,npp,nq,nr,op,ib
          real(dp), pointer, dimension(:,:,:)::fi
          real(dp) :: e,q1,q2,q3,c_star,v1,v2,v3,v4,v5,v11,v22
          real(dp) :: v33,v44,v55,s1,s2,s3,a1,a2,a3,w1,w2,w3

          do ib=1,nbp

              select case (op)
                case (1)
                  fi => dom(ib)%u
                  npp = dom(ib)%niul
                  nq = dom(ib)%njul
                  nr = dom(ib)%nkul
                case (2)
                  fi => dom(ib)%v
                  npp = dom(ib)%nivl
                  nq = dom(ib)%njvl
                  nr = dom(ib)%nkvl
                case (3)
                  fi => dom(ib)%w
                  npp = dom(ib)%niwl
                  nq = dom(ib)%njwl
                  nr = dom(ib)%nkwl
                case (14)  ! (phi)
                  fi => dom(ib)%phi
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
                case (15)  ! (phi_reinit)
                  fi => dom(ib)%phi_reinit
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
                case (16)  ! (phi_new)
                  fi => dom(ib)%phi_new
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
                case (17)  ! (phi_init)
                  fi => dom(ib)%phi_init
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
              end select

              if (dom(ib)%kprev<0 .and. dom(ib)%bc_bottom/=5) then
              do i=1,npp
                  do j=1,nq
                      fi(i,j,2)    = fi(i,j,3)
                      fi(i,j,1)    = fi(i,j,2)
                  end do
              end do
              else if (dom(ib)%knext<0 .and. dom(ib)%bc_top/=5) then
              do i=1,npp
                  do j=1,nq
                      fi(i,j,nr-1) = fi(i,j,nr-2)
                      fi(i,j,nr)   = fi(i,j,nr-1)
                  end do
              end do
              end if
!
! Compute 1st order spatial derivative
!
              do k=1,nr-1
                  do i=1,npp
                      do j=1,nq
                          dom(ib)%d1(i,j,k) = (fi(i,j,k+1) - fi(i,j,k))/dom(ib)%dz
                      end do
                  end do
              end do
!
! Compute 5th order derivative (-ve direction)
!
              dom(ib)%dphi_dzplus = 0.0_dp

              do k=1,nr-6
                  do i=1,npp
                      do j=1,nq

                          l=k

                          v1 = dom(ib)%d1(i,j,l+5)
                          v2 = dom(ib)%d1(i,j,l+4)
                          v3 = dom(ib)%d1(i,j,l+3)
                          v4 = dom(ib)%d1(i,j,l+2)
                          v5 = dom(ib)%d1(i,j,l+1)

                          v11 = v1**2
                          v22 = v2**2
                          v33 = v3**2
                          v44 = v4**2
                          v55 = v5**2

                          e = (1.0e-6_dp) * max(v11,v22,v33,v44,v55) + 1E-99_dp

                          s1 = (13.0_dp/12.0_dp)*(v1-2.0_dp*v2+v3)**2 + 0.25_dp*(v1-4.0_dp*v2+3.0_dp*v3)**2
                          s2 = (13.0_dp/12.0_dp)*(v2-2.0_dp*v3+v4)**2 + 0.25_dp*(v2-v4)**2
                          s3 = (13.0_dp/12.0_dp)*(v3-2.0_dp*v4+v5)**2 + 0.25_dp*(3.0_dp*v3-4.0_dp*v4+v5)**2

                          a1 = 0.1_dp/(e+s1)**2
                          a2 = 0.6_dp/(e+s2)**2
                          a3 = 0.3_dp/(e+s3)**2

                          w1 = a1/(a1+a2+a3)
                          w2 = a2/(a1+a2+a3)
                          w3 = a3/(a1+a2+a3)

                          dom(ib)%dphi_dzplus(i,j,k+3) = w1*(v1*(1.0_dp/3.0_dp)-v2*(7.0_dp/6.0_dp)+ &
                    v3*(11.0_dp/6.0_dp)) + w2*(-v2*(1.0_dp/6.0_dp)+v3*(5.0_dp/6.0_dp)+v4*(1.0_dp/3.0_dp))+ &
                    w3*(v3*(1.0_dp/3.0_dp)+v4*(5.0_dp/6.0_dp)-v5*(1.0_dp/6.0_dp))

                      end do
                  end do
              end do

          end do

          return
      end subroutine hj_weno_dzplus_3d
!######################################################################
      subroutine hj_weno_dzminus_3d(op)
!######################################################################
          use vars
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64

          implicit none

          integer :: i,j,k,l,k_star,npp,nq,nr,op,ib
          real(dp), pointer, dimension(:,:,:)::fi
          real(dp) :: q1,q2,q3,c_star,v1,v2,v3,v4,v5,v11,v22
          real(dp) :: v33,v44,v55,s1,s2,s3,a1,a2,a3,w1,w2,w3
          real(dp) :: e

          do ib=1,nbp

              select case (op)
                case (1)
                  fi => dom(ib)%u
                  npp = dom(ib)%niul
                  nq = dom(ib)%njul
                  nr = dom(ib)%nkul
                case (2)
                  fi => dom(ib)%v
                  npp = dom(ib)%nivl
                  nq = dom(ib)%njvl
                  nr = dom(ib)%nkvl
                case (3)
                  fi => dom(ib)%w
                  npp = dom(ib)%niwl
                  nq = dom(ib)%njwl
                  nr = dom(ib)%nkwl
                case (14)  ! (phi)
                  fi => dom(ib)%phi
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
                case (15)  ! (phi_reinit)
                  fi => dom(ib)%phi_reinit
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
                case (16)  ! (phi_new)
                  fi => dom(ib)%phi_new
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
                case (17)  ! (phi_init)
                  fi => dom(ib)%phi_init
                  npp = dom(ib)%ttc_i
                  nq = dom(ib)%ttc_j
                  nr = dom(ib)%ttc_k
              end select

              if (dom(ib)%kprev<0 .and. dom(ib)%bc_bottom/=5) then
              do i=1,npp
                  do j=1,nq
                      fi(i,j,2)    = fi(i,j,3)
                      fi(i,j,1)    = fi(i,j,2)
                  end do
              end do
              else if (dom(ib)%knext<0 .and. dom(ib)%bc_top/=5) then
              do i=1,npp
                  do j=1,nq
                      fi(i,j,nr-1) = fi(i,j,nr-2)
                      fi(i,j,nr)   = fi(i,j,nr-1)
                  end do
              end do
              end if
!
! Compute 1st order spatial derivative
!
              do k=1,nr-1
                  do i=1,npp
                      do j=1,nq
                          dom(ib)%d1(i,j,k) = (fi(i,j,k+1) - fi(i,j,k))/dom(ib)%dz
                      end do
                  end do
              end do
!
! Compute 5th order derivative (-ve direction)
!
              dom(ib)%dphi_dzminus = 0.0_dp

              do k=1,nr-6
                  do i=1,npp
                      do j=1,nq

                          l=k-1

                          v1 = dom(ib)%d1(i,j,l+1)
                          v2 = dom(ib)%d1(i,j,l+2)
                          v3 = dom(ib)%d1(i,j,l+3)
                          v4 = dom(ib)%d1(i,j,l+4)
                          v5 = dom(ib)%d1(i,j,l+5)

                          v11 = v1**2
                          v22 = v2**2
                          v33 = v3**2
                          v44 = v4**2
                          v55 = v5**2

                          e = (1.0e-6_dp) * max(v11,v22,v33,v44,v55) + 1E-99_dp

                          s1 = (13.0_dp/12.0_dp)*(v1-2.0_dp*v2+v3)**2 + 0.25_dp*(v1-4.0_dp*v2+3.0_dp*v3)**2
                          s2 = (13.0_dp/12.0_dp)*(v2-2.0_dp*v3+v4)**2 + 0.25_dp*(v2-v4)**2
                          s3 = (13.0_dp/12.0_dp)*(v3-2.0_dp*v4+v5)**2 + 0.25_dp*(3.0_dp*v3-4.0_dp*v4+v5)**2

                          a1 = 0.1_dp/(e+s1)**2
                          a2 = 0.6_dp/(e+s2)**2
                          a3 = 0.3_dp/(e+s3)**2

                          w1 = a1/(a1+a2+a3)
                          w2 = a2/(a1+a2+a3)
                          w3 = a3/(a1+a2+a3)

                          dom(ib)%dphi_dzminus(i,j,k+3) = w1*(v1*(1.0_dp/3.0_dp)-v2*(7.0_dp/6.0_dp)+ &
                    v3*(11.0_dp/6.0_dp)) + w2*(-v2*(1.0_dp/6.0_dp)+v3*(5.0_dp/6.0_dp)+v4*(1.0_dp/3.0_dp))+ &
                    w3*(v3*(1.0_dp/3.0_dp)+v4*(5.0_dp/6.0_dp)-v5*(1.0_dp/6.0_dp))

                      end do
                  end do
              end do

          end do

          return
      end subroutine hj_weno_dzminus_3d
!#####################################################################
      subroutine heaviside
!######################################################################
          use vars
          use module_LSM
          use multiflow3d_mpi
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64

          implicit none
          integer :: i,j,k,n_epsl,ib,tti,ttj,ttk
          real(dp) :: epsl
          real(dp), parameter :: pi = 3.14159265359_dp
!
! Define an infinitely differentiable smoothed heaviside function h_phi
!
          n_epsl = 2

          do ib=1,nbp

              epsl = n_epsl*max(dom(ib)%dx,dom(ib)%dy,dom(ib)%dz)

              do k=dom(ib)%ksp-pl,dom(ib)%kep+pl
                  do i=dom(ib)%isp-pl,dom(ib)%iep+pl
                      do j=dom(ib)%jsp-pl,dom(ib)%jep+pl

                          if (dom(ib)%phi(i,j,k)<(-1.0_dp*epsl)) then
                            dom(ib)%h_phi(i,j,k) = 0.0_dp
                          end if  ! h_phi=0 above free surface

                          if (dom(ib)%phi(i,j,k)>(epsl)) dom(ib)%h_phi(i,j,k) = 1.0_dp  ! h_phi=1.0_dp below free surface

                          if (abs(dom(ib)%phi(i,j,k))<=epsl) then
!
! Transition zone across free surface (2 grid cells width either side)
!
                          dom(ib)%h_phi(i,j,k)=0.5_dp*(1.0_dp+dom(ib)%phi(i,j,k)/ &
                    epsl+1.0_dp/pi*sin(pi*dom(ib)%phi(i,j,k)/epsl))

                          end if

                          dom(ib)%dens(i,j,k)=densg+(densl-densg)*dom(ib)%h_phi(i,j,k)  !dens = densg above free surface, densl below
                          dom(ib)%mu(i,j,k)=mug+(mul-mug)*dom(ib)%h_phi(i,j,k)  !mu = mug above free surface, mul below

                      end do
                  end do

              end do

          end do

          call bound_lsm(18)  !(dens)
          call bound_lsm(19)  !(mu)

          return
      end subroutine heaviside
!#####################################################################
