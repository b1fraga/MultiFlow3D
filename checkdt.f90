!##########################################################################
      subroutine checkdt
!##########################################################################
          use vars
          use multiflow3d_mpi
          use multidata
          use module_LSM
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          integer :: i,j,k,ib
          real(dp) :: dxx,dyy,dzz,umax,vmax,wmax,dtmax
          real(dp) :: dtmax1,dtvisc,dtvisc1,dtthr
          real(dp) :: uc,vc,wc
          real(dp) :: buffer_umax,buffer_vmax,buffer_wmax
          real(dp) :: buffer_dtmax,dt1,small
          real(dp) :: Cu,Cv,Cw

          umax=0.0_dp
          vmax=0.0_dp
          wmax=0.0_dp
          small=1e-30_dp

          MPI_FLT = MPI_DOUBLE_PRECISION

          do ib=1,nbp
              do i=dom(ib)%isu,dom(ib)%ieu
                  do j=dom(ib)%jsu,dom(ib)%jeu
                      do k=dom(ib)%ksu,dom(ib)%keu
                          umax=max(umax,abs(dom(ib)%u(i,j,k)))
                      end do
                  end do
              end do
          end do

          buffer_umax=umax
          call MPI_ALLREDUCE (buffer_umax,umax,1,MPI_FLT,MPI_MAX, &
                     MPI_COMM_WORLD,ierr )

          do ib=1,nbp
              do i=dom(ib)%isv,dom(ib)%iev
                  do j=dom(ib)%jsv,dom(ib)%jev
                      do k=dom(ib)%ksv,dom(ib)%kev
                          vmax=max(vmax,abs(dom(ib)%v(i,j,k)))
                      end do
                  end do
              end do
          end do

          buffer_vmax=vmax
          call MPI_ALLREDUCE (buffer_vmax,vmax,1,MPI_FLT,MPI_MAX, &
                     MPI_COMM_WORLD,ierr )

          do ib=1,nbp
              do i=dom(ib)%isw,dom(ib)%iew
                  do j=dom(ib)%jsw,dom(ib)%jew
                      do k=dom(ib)%ksw,dom(ib)%kew
                          wmax=max(wmax,abs(dom(ib)%w(i,j,k)))
                      end do
                  end do
              end do
          end do

          buffer_wmax=wmax
          call MPI_ALLREDUCE (buffer_wmax,wmax,1,MPI_FLT,MPI_MAX, &
                     MPI_COMM_WORLD,ierr )


          dtvisc=1e10
          if(SGS) then
          do ib=1,nbp
              dxx=dom(ib)%dx*dom(ib)%dx
              dyy=dom(ib)%dy*dom(ib)%dy
              dzz=dom(ib)%dz*dom(ib)%dz
              do i=dom(ib)%isp,dom(ib)%iep
                  do j=dom(ib)%jsp,dom(ib)%jep
                      do k=dom(ib)%ksp,dom(ib)%kep
                          uc=0.5_dp*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                          vc=0.5_dp*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                          wc=0.5_dp*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                          dtvisc1=1.0_dp/( abs(uc/dom(ib)%dx)+ &
                    abs(vc/dom(ib)%dy)+abs(wc/dom(ib)%dz)+ &
                    2.0_dp*dom(ib)%vis(i,j,k)*(1.0_dp/dxx+1.0_dp/dyy+1.0_dp/dzz)+small)
                          dtvisc=min(dtvisc,dtvisc1)
                          if(LENERGY) then
                          dtthr=0.5_dp*Re*Pr/(1.0_dp/dxx + 1.0_dp/dyy + 1.0_dp/dzz)
                          dtvisc=min(dtvisc,dtthr)
                          end if
                      end do
                  end do
              end do
          end do
          else
          do ib=1,nbp
              dxx=dom(ib)%dx*dom(ib)%dx
              dyy=dom(ib)%dy*dom(ib)%dy
              dzz=dom(ib)%dz*dom(ib)%dz
              dtvisc1=1.0_dp/(1.0_dp/(dxx) + 1.0_dp/(dyy)+ 1.0_dp/(dzz))*Re/2.0_dp
              dtvisc=min(dtvisc,dtvisc1)
              if(LENERGY) then
              dtthr=0.5_dp*Re*Pr/(1.0_dp/dxx + 1.0_dp/dyy + 1.0_dp/dzz)
              dtvisc=min(dtvisc,dtthr)
              end if
          end do
          end if


          dtmax=1e10
          do ib=1,nbp
              dtmax1=min(dom(ib)%dx/umax,dom(ib)%dy/vmax,dom(ib)%dz/wmax)
              dtmax = min(dtmax1,dtmax)
          end do

          dtmax = min(dtvisc,dtmax)
          dtmax = safety_factor * dtmax
          dt=min(dt*1.1_dp,dtmax)

          if (L_LSM)  then
          do ib=1,nbp
              dtvisc=max(mul/densl,mug/densg)* &
        (2.0_dp/(dxx)+2.0_dp/(dyy)+2.0_dp/(dzz))
              Cu=1.0_dp/((umax/dom(ib)%dx+dtvisc)+sqrt((umax/dom(ib)%dx+ &
        dtvisc)**2+4.0_dp*abs(gx)/dom(ib)%dx))
              Cv=1.0_dp/((vmax/dom(ib)%dy+dtvisc)+sqrt((vmax/dom(ib)%dy+ &
        dtvisc)**2+4.0_dp*abs(gy)/dom(ib)%dy))
              Cw=1.0_dp/((wmax/dom(ib)%dz+dtvisc)+sqrt((wmax/dom(ib)%dz+ &
        dtvisc)**2+4.0_dp*abs(gz)/dom(ib)%dz))
              dt = 2.0_dp*safety_factor*min(Cu,Cv,Cw)
          end do
          end if

          buffer_dtmax=dt
!        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
          call MPI_ALLREDUCE (buffer_dtmax,dt1,1,MPI_FLT,MPI_MIN, &
                     MPI_COMM_WORLD,ierr )

          dt=dt1

          if(itime/=itime_start) then
          if(dt<dtavg*0.1_dp) then
          print*,"#*#*#*#*#*# dt becomes smaller, check result!!!!"
          if (myrank==0) then
            write(numfile,*) "#*#*#*# dt becomes smaller, check result!!!!"
          end if
          call tecgrid(itime)
          call tecplot_p(itime)
          call tecplot_u(itime)
          call tecplot_v(itime)
          call tecplot_w(itime)
          call tecbin(itime)
          if(myrank==0) then
          open (unit=101, file="final_ctime.dat")
          write (101,"(i8,3F15.6)") &
    ntime,ctime,forcn,qstpn,count
          close(101)
          end if
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
          call MPI_BARRIER (MPI_COMM_WORLD,ierr)
          stop
          end if
          dtsum=dtsum+dt
          else
          dtsum=dt
          end if
          dtavg=dtsum/(itime-itime_start+1)

          return
      end subroutine checkdt
!##########################################################################
