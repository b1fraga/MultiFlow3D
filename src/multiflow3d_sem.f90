module multiflow3d_sem
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use vars, only: dt, g_dy, g_dz, itmax_sem, ne_sem, ti_sem, ubulk, &
                  yen, yst, zen, zst
  use multidata, only: idom, jdom, kdom, ycor, zcor
  use multiflow3d_mpi, only: myrank
  implicit none
  private

  public :: sem_initial

  contains

    subroutine sem_initial
    !this program implements the sem method described in n. jarrin thesis, chp. 4
      real(dp) :: vol,ly,lz,enne,xmin,xmax,ymin
      real(dp) :: ymax,zmin,zmax,pi,u0,hu,riz,rdivz,uave(3)
      real(dp) :: sigma_value,maxvsem,minvsem
      real(dp) :: rand_num
      integer :: divy,divz,iy,iz,ii,n,i,j,m,it,iglobal
      integer,allocatable,dimension(:)::lsy,lsz,ley,lez
      real(dp),dimension(:,:,:), allocatable :: vsem ,usem
      real(dp),dimension(:,:), allocatable :: x_eddy,epsilo,molt
      real(dp),dimension(:,:), allocatable :: sigma
      real(dp),dimension(:), allocatable ::  ksem
      real(dp),dimension(:) :: x_point(3),reynolds(6)
      real(dp),dimension(:) :: temp(3),temp2(3)
      real(dp), dimension(:,:) ::  r(3,3)
      character(len=44) :: fileglobal
      integer,allocatable:: elemyst(:),elemyen(:),elemzst(:),elemzen(:)
      integer,allocatable:: iddom(:),ljdom(:),lkdom(:)
      !the box dimensions are difined as [xlenght] * [ly] * [lz]
      ![divx], [divy] & [divz] are the numbers of spacial points.
      pi = 2*acos(0.0d0)
      ly  = yen-yst
      lz= zen-zst
      u0 = ubulk

      allocate(iddom(jdom*kdom),ljdom(jdom*kdom),lkdom(jdom*kdom))
      allocate(elemyst(jdom*kdom),elemzst(kdom*jdom))
      allocate(elemyen(jdom*kdom),elemzen(kdom*jdom))
      allocate(ley(jdom),lez(kdom),lsy(jdom),lsz(kdom))

      ljdom=0
      lkdom=0
      divy = 0
      divz = 0
      elemyst=0
      elemyen=0
      elemzst=0
      elemzen=0
      i=0
      ley=0
      lez=0
      lsy=0
      lsz=0

      !id of the block
      do n=1,kdom
         do j=1,jdom
            i=i+1
            iddom(i)=idom*(j-1)+(jdom*idom*(n-1))
         end do
      end do

      !division on the y and z directions
      do j=1,jdom
         ley(j)=ley(j-1)+ &
                nint((ycor(iddom(j),2)-ycor(iddom(j),1))/g_dy) +1
         if(j==1) then
           lsy(j)=1
         else
           lsy(j)=ley(j-1)+1
         end if
         divy = divy + nint((ycor(iddom(j),2)-ycor(iddom(j),1))/g_dy) +1
      end do

      do n=1,kdom
         lez(n)=lez(n-1)+ &
                nint((zcor(iddom(n),2)-zcor(iddom(n),1))/g_dz)+1
         if(n==1) then
           lsz(n)=1
         else
           lsz(n)=lez(n-1)+1
         end if
         divz = divz + nint((zcor(iddom(n),2)-zcor(iddom(n),1))/g_dz) +1
      end do

      !the divisions for each of the domains is determined:
      i=0
      do n=1,kdom
         do j=1,jdom
            i=i+1
            elemyst(i) = lsy(j)
            elemzst(i) = lsz(n)
            elemyen(i) = ley(j)
            elemzen(i) = lez(n)
            ljdom(i) = elemyen(i) -  elemyst(i) + 1
            lkdom(i) = elemzen(i) -  elemzst(i) + 1
         end do
      end do

      write(6,*)"divisions :",divy,divz
      write(6,*)"# blocks  :",jdom,kdom
      write(6,*)lsy(1),ley(1),lsy(2),ley(2),ley(2)-lsy(2)+1
      write(6,*)lsz(1),lez(1),lsz(2),lez(2),lez(2)-lsz(2)+1

      !the amplitude of the vortices. this can be changed to have larger structures!!!!!!
      sigma_value=min(16.d0*g_dy,lz/2.d0,ly/2.d0)  !isotropic

      !the number of eddies is the inlet surface divided by the surface of each turbulent spot
      ne_sem = (ly*lz)/sigma_value**2

      ![n](integer) and [enne](real(dp)) represent the number of eddies
      n = int(ne_sem)
      enne = real(n)

      if(myrank==0) write(6,*) "the number of sem eddies is :", n

      !  [reynolds(6)] is a vector with the six elements of reynolds stresses.
      !  |reynolds(1)  reynolds(2)  reynolds(4)|
      !  |reynolds(2)  reynolds(3)  reynolds(5)|
      !  |reynolds(4)  reynolds(5)  reynolds(6)|
      reynolds=[(ti_sem*u0)**2, 0.0d0,(ti_sem*u0)**2, 0.0d0, &
                    0.0d0,(ti_sem*u0)**2]   ![m/s]

      !allocation of the eddies vector.
      ![vsem(divx,divy,divz,3)] is the instantaneous velocity vector in the point with
      !x,y,z components
      ![x_eddy(3,n)] is the n-th eddy location x,y,z
      ![epsilo(3,n)] is the n-th eddy intensity in x,y,z
      ![molt(3,n)]   is the matrix product [r(3,3)] * [epsilo(3,n)]
      ![ksem(n)]     is a vector used to check whitch eddy is outside the box after the convection

      allocate(vsem(divy,divz,3),usem(divy,divz,3),ksem(n))
      allocate(x_eddy(3,n),epsilo(3,n),molt(3,n),sigma(divy,divz))
      print *, "total time interval = ", dt * itmax_sem, " [s]"

      xmin=1.d8
      xmax=1.d-8
      ymin=1.d8
      ymax=1.d-8
      zmin=1.d8
      zmax=1.d-8

      !definition of eddy length scale and initial velocity field
      do iy=1,divy
         do iz=1,divz
            sigma(iy,iz)=sigma_value    !!!  min(8.0*g_dy,0.20d0) !isotropic
            !set the inlet velocity prof.
            usem(iy,iz,:)=[ u0 ,0.d0,0.d0]
            !   if(uprof_sem.eq.1) usem(iy,iz,:)=(/ u0 ,0.d0,0.d0/)      !elli
            !       if(uprof_sem.eq.15) then
            !        riz=iz                         !elli
            !      usem(iy,iz,:)=(/ u0-0.04324*exp(-11.1764*riz), 0.d0,0.d0/)
            !   endif
            uave(:) = uave(:) + usem(iy,iz,:)
            !calculation of the eddy box parameters
            xmin=min(0.0d0-sigma(iy,iz),xmin)
            xmax=max(0.0d0+sigma(iy,iz),xmax)
            ymin=min(0.0d0-sigma(iy,iz),ymin)
            ymax=max(ly   +sigma(iy,iz),ymax)
            zmin=min(0.0d0-sigma(iy,iz),zmin)
            zmax=max(lz   +sigma(iy,iz),zmax)
         end do
      end do

      uave = uave / (divy * divz)
      vol  = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
      !generation of the eddy location inside the box and initialization of the [ksem] vector
      do ii=1,n
         call random_number(rand_num)
         x_eddy(1,ii) = (xmax - xmin) * rand_num + xmin
         call random_number(rand_num)
         x_eddy(2,ii) = (ymax - ymin) * rand_num + ymin
         call random_number(rand_num)
         x_eddy(3,ii) = (zmax - zmin) * rand_num + zmin
         ksem(ii) = 0
         !initialization of the intensities. for every direction the average intensity value is calculated and
         !it is forced to be lower than the [vlim] value
         call random_number(rand_num)
         epsilo(1,ii) = (rand_num*2.0d0 - 1.0d0)
         call random_number(rand_num)
         epsilo(2,ii) = (rand_num*2.0d0 - 1.0d0)
         call random_number(rand_num)

         epsilo(3,ii) = (rand_num*2.0d0 - 1.0d0)
      end do

      !initialization of the [r(3,3)] matrix with the cholensky decomposition
      !of the reynolds stress tensor
      r = 0
      r(1,1) = dsqrt(reynolds(1))
      r(2,1) = reynolds(2) / r(1,1)
      r(2,2) = dsqrt(reynolds(3) - r(2,1)*r(2,1))
      r(3,1) = reynolds(4) / r(1,1)
      r(3,2) = (reynolds(5) - r(2,1)*r(3,1)) / r(2,2)
      r(3,3) = dsqrt(reynolds(6) - r(3,1)*r(3,1) - r(3,2)*r(3,2))
      !beginning of time iterations
      do it=1,itmax_sem         !parallelize this loop
         if(mod(it,50)==0) then
            write(*,*)"iteration ",it,"in progress.time: ",(it-1) * dt,"[s]"
         end if

         !printings of global velocity and of convection velocity
         do i=1,jdom*kdom
            iglobal=500+i
            write(fileglobal,"(a13,i4.4,a1,i6.6,a4)") &
                 "inflow/inlet_",iddom(i),"_",it,".dat"
            open(unit=iglobal,file=fileglobal,status="unknown", &
                 action="write")
            write (iglobal,*)"variables=up,vp,wp"
            write (iglobal,*) &
            "zone "," i=",ljdom(i),","," j=",lkdom(i),", k= ",1," f=point"
         end do

         molt = matmul(r,epsilo)  ! !aij*epsij   matrix multiplication
         maxvsem=0.d0
         minvsem=1000000.d0

         !------beginning of spatial iteration
         do iy = 1,divy
            do iz = 1,divz
               !x_point = grid point coordinates
               x_point=[0.0d0,iy*ly/divy+ymin,iz*lz/divz+zmin]
               vsem(iy,iz,:)=[ 0.d0, 0.d0,  0.d0 ]
               !------------beginning of eddies iterations
               do ii=1,n
                  temp(:) = dabs(x_point(:) - x_eddy(:,ii))
                  if (temp(1)<sigma(iy,iz) .and. temp(2)<sigma(iy,iz) .and. &
                       temp(3)<sigma(iy,iz)) then

                     vsem(iy,iz,:)=vsem(iy,iz,:)+molt(:,ii) &
                          *(dsqrt(1.5d0)**3.0d0)*dsqrt(vol)/dsqrt(sigma(iy,iz)**3) &
                          *(1.0d0- dabs(x_point(:) - x_eddy(:,ii))/sigma(iy,iz)) &
                          *(1.0d0- dabs(x_point(:) - x_eddy(:,ii))/sigma(iy,iz)) &
                          *(1.0d0- dabs(x_point(:) - x_eddy(:,ii))/sigma(iy,iz))
                     !f(x)=sqrt(3/2)*(1-abs(x)) jarrin et al. 2009
                  end if
               end do
               ! instananeous velocity=mean velocity + sem fluctuation velocity
               vsem(iy,iz,:) = vsem(iy,iz,:) / dsqrt(enne)
               if(vsem(iy,iz,1)>maxvsem) maxvsem=vsem(iy,iz,1)
               if(vsem(iy,iz,1)<=minvsem) minvsem=vsem(iy,iz,1)
            end do
         end do

         !    write(6,*) '========='
         write(6,"(i6,2f15.6)")it,maxvsem,minvsem

         !------end of spatial iterations
         do i=1,jdom*kdom
            iglobal=500+i
            do m=elemzst(i),elemzen(i)
               do j=elemyst(i),elemyen(i)
                  write(500+i,"(3e15.6)") vsem(j,m,:)                !turn down precision for large files
               end do
            end do
            close (unit=iglobal)
         end do

         !--------beginning of eddies convection iterations
         do ii=1,n
            !re-calculation of the eddies position. if any eddy goes beyond the box limits
            !it is restarted at the surface facing the exit.
            x_eddy(1,ii) = x_eddy(1,ii) + uave(1) * dt
            x_eddy(2,ii) = x_eddy(2,ii) + uave(2) * dt
            x_eddy(3,ii) = x_eddy(3,ii) + uave(3) * dt
            !after the eddies convections are necessary some tests to check if any eddy is now outside
            !the sem box defined earlier
            !when a eddy is re-generate the ksem(i) factor assume the 1 value. this value is used later to
            !generate a new intensity for the new eddy
            if (x_eddy(1,ii) > xmax) then
               x_eddy(1,ii) = xmin
               call random_number(rand_num)
               x_eddy(2,ii) = (ymax - ymin) * rand_num + ymin
               call random_number(rand_num)
               x_eddy(3,ii) = (zmax - zmin) * rand_num + zmin
               ksem(ii) = 1
            else if (x_eddy(1,ii) < xmin) then
               x_eddy(1,ii) = xmax
               call random_number(rand_num)
               x_eddy(2,ii) = (ymax - ymin) * rand_num + zmin
               call random_number(rand_num)
               x_eddy(3,ii) = (zmax - zmin) * rand_num + zmin
               ksem(ii) = 1
            else if (x_eddy(2,ii) > ymax) then
               call random_number(rand_num)
               x_eddy(1,ii) = (xmax - xmin) * rand_num + xmin
               x_eddy(2,ii) = ymin
               call random_number(rand_num)
               x_eddy(3,ii) = (zmax - zmin) * rand_num + zmin
               ksem(ii) = 1
            else if (x_eddy(2,ii) < ymin) then
               call random_number(rand_num)
               x_eddy(1,ii) = (xmax - xmin) * rand_num + xmin
               x_eddy(2,ii) = ymax
               call random_number(rand_num)
               x_eddy(3,ii) = (zmax - zmin) * rand_num + zmin
               ksem(ii) = 1
            else if  (x_eddy(3,ii) > zmax) then
               call random_number(rand_num)
               x_eddy(1,ii) = (xmax - xmin) * rand_num + xmin
               call random_number(rand_num)
               x_eddy(2,ii) = (ymax - ymin) * rand_num + zmin
               x_eddy(3,ii) = zmin
               ksem(ii) = 1
            else if (x_eddy(3,ii) < zmin) then
               call random_number(rand_num)
               x_eddy(1,ii) = (xmax - xmin) * rand_num + xmin
               call random_number(rand_num)
               x_eddy(2,ii) = (ymax - ymin) * rand_num + zmin
               x_eddy(3,ii) = zmax
               ksem(ii) = 1
            end if
            !intensity generation for the re-created eddies. we are using the ksem factor as explained fator.
            if (ksem(ii)== 1) then
               call random_number(rand_num)
               epsilo(3,ii) = (rand_num*2.0d0 - 1.0d0)
               call random_number(rand_num)
               epsilo(2,ii) = (rand_num*2.0d0 - 1.0d0)
               call random_number(rand_num)
               epsilo(1,ii) = (rand_num*2.0d0 - 1.0d0)
            end if
            ksem(ii) = 0
         end do
         !---------end of eddies convections iterations
      end do
      !----end of time iterations

    end subroutine sem_initial

end module multiflow3d_sem
