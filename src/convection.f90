!##########################################################################
      subroutine convection
!##########################################################################
          use vars
          use multidata
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          integer :: i,j,k,ib,ii,jj,kk
          real(dp) :: du2dx,dv2dy,dw2dz,vel
          real(dp) :: duvdx,duvdy,duwdx,duwdz,dvwdy,dvwdz
          real(dp) :: up12,um12,vp12,vm12,wp12,wm12
          real(dp) :: uijk,vijk,wijk
          real(dp) :: dxx,dyy,dzz,dudx,dudy,dudz
          real(dp) :: dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
          real(dp) :: sp2,sp1,sn,sm1,sm2

          do ib=1,nbp
              dom(ib)%ustar=dom(ib)%u
              dom(ib)%vstar=dom(ib)%v
              dom(ib)%wstar=dom(ib)%w
          end do

          if (differencing==3) then
          call HJ_WENO_dx(1)
          call HJ_WENO_dy(1)
          call HJ_WENO_dz(1)
          end if

          do ib=1,nbp

              do k=dom(ib)%ksu,dom(ib)%keu
                  do i=dom(ib)%isu,dom(ib)%ieu
                      do j=dom(ib)%jsu,dom(ib)%jeu

                          if(differencing==1) then

                          up12=0.5_dp*(dom(ib)%u(i,j,k)+dom(ib)%u(i+1,j,k))
                          um12=0.5_dp*(dom(ib)%u(i,j,k)+dom(ib)%u(i-1,j,k))
                          du2dx=(up12**2-um12**2)/dom(ib)%dx
                          duvdy=0.25_dp*( (dom(ib)%v(i,j,k)+dom(ib)%v(i+1,j,k))* &
                    (dom(ib)%u(i,j,k)+dom(ib)%u(i,j+1,k)) - &
                    (dom(ib)%v(i,j-1,k)+dom(ib)%v(i+1,j-1,k))* &
                    (dom(ib)%u(i,j,k)+dom(ib)%u(i,j-1,k)) )/dom(ib)%dy
                          duwdz=0.25_dp*( (dom(ib)%w(i,j,k)+dom(ib)%w(i+1,j,k))* &
                    (dom(ib)%u(i,j,k)+dom(ib)%u(i,j,k+1)) - &
                    (dom(ib)%w(i,j,k-1)+dom(ib)%w(i+1,j,k-1))* &
                    (dom(ib)%u(i,j,k)+dom(ib)%u(i,j,k-1)) )/dom(ib)%dz

                          else if(differencing==2) then
                          ii=i+2
                          sp2=((-dom(ib)%u(ii+1,j,k)+9.0_dp*dom(ib)%u(ii,j,k)+ &
                    9.0_dp*dom(ib)%u(ii-1,j,k)-dom(ib)%u(ii-2,j,k))/16.0_dp)**2
                          ii=i+1
                          sp1=((-dom(ib)%u(ii+1,j,k)+9.0_dp*dom(ib)%u(ii,j,k)+ &
                    9.0_dp*dom(ib)%u(ii-1,j,k)-dom(ib)%u(ii-2,j,k))/16.0_dp)**2
                          ii=i
                          sn=((-dom(ib)%u(ii+1,j,k)+9.0_dp*dom(ib)%u(ii,j,k)+ &
                    9.0_dp*dom(ib)%u(ii-1,j,k)-dom(ib)%u(ii-2,j,k))/16.0_dp)**2
                          ii=i-1
                          sm1=((-dom(ib)%u(ii+1,j,k)+9.0_dp*dom(ib)%u(ii,j,k)+ &
                    9.0_dp*dom(ib)%u(ii-1,j,k)-dom(ib)%u(ii-2,j,k))/16.0_dp)**2
                          du2dx=(-sp2+27.0_dp*sp1-27.0_dp*sn+sm1)/(24.0_dp*dom(ib)%dx)

                          jj=j+1
                          sp1=((-dom(ib)%u(i,jj+2,k)+9.0_dp*dom(ib)%u(i,jj+1,k)+ &
                    9.0_dp*dom(ib)%u(i,jj,k)-dom(ib)%u(i,jj-1,k))/16.0_dp)* &
                    ((-dom(ib)%v(i+2,jj,k)+9.0_dp*dom(ib)%v(i+1,jj,k)+ &
                    9.0_dp*dom(ib)%v(i,jj,k)-dom(ib)%v(i-1,jj,k))/16.0_dp)
                          jj=j
                          sn=((-dom(ib)%u(i,jj+2,k)+9.0_dp*dom(ib)%u(i,jj+1,k)+ &
                    9.0_dp*dom(ib)%u(i,jj,k)-dom(ib)%u(i,jj-1,k))/16.0_dp)* &
                    ((-dom(ib)%v(i+2,jj,k)+9.0_dp*dom(ib)%v(i+1,jj,k)+ &
                    9.0_dp*dom(ib)%v(i,jj,k)-dom(ib)%v(i-1,jj,k))/16.0_dp)
                          jj=j-1
                          sm1=((-dom(ib)%u(i,jj+2,k)+9.0_dp*dom(ib)%u(i,jj+1,k)+ &
                    9.0_dp*dom(ib)%u(i,jj,k)-dom(ib)%u(i,jj-1,k))/16.0_dp)* &
                    ((-dom(ib)%v(i+2,jj,k)+9.0_dp*dom(ib)%v(i+1,jj,k)+ &
                    9.0_dp*dom(ib)%v(i,jj,k)-dom(ib)%v(i-1,jj,k))/16.0_dp)
                          jj=j-2
                          sm2=((-dom(ib)%u(i,jj+2,k)+9.0_dp*dom(ib)%u(i,jj+1,k)+ &
                    9.0_dp*dom(ib)%u(i,jj,k)-dom(ib)%u(i,jj-1,k))/16.0_dp)* &
                    ((-dom(ib)%v(i+2,jj,k)+9.0_dp*dom(ib)%v(i+1,jj,k)+ &
                    9.0_dp*dom(ib)%v(i,jj,k)-dom(ib)%v(i-1,jj,k))/16.0_dp)
                          duvdy=(-sp1+27.0_dp*sn-27.0_dp*sm1+sm2)/(24.0_dp*dom(ib)%dy)

                          kk=k+1
                          sp1=((-dom(ib)%u(i,j,kk+2)+9.0_dp*dom(ib)%u(i,j,kk+1)+ &
                    9.0_dp*dom(ib)%u(i,j,kk)-dom(ib)%u(i,j,kk-1))/16.0_dp)* &
                    ((-dom(ib)%w(i+2,j,kk)+9.0_dp*dom(ib)%w(i+1,j,kk)+ &
                    9.0_dp*dom(ib)%w(i,j,kk)-dom(ib)%w(i-1,j,kk))/16.0_dp)
                          kk=k
                          sn=((-dom(ib)%u(i,j,kk+2)+9.0_dp*dom(ib)%u(i,j,kk+1)+ &
                    9.0_dp*dom(ib)%u(i,j,kk)-dom(ib)%u(i,j,kk-1))/16.0_dp)* &
                    ((-dom(ib)%w(i+2,j,kk)+9.0_dp*dom(ib)%w(i+1,j,kk)+ &
                    9.0_dp*dom(ib)%w(i,j,kk)-dom(ib)%w(i-1,j,kk))/16.0_dp)
                          kk=k-1
                          sm1=((-dom(ib)%u(i,j,kk+2)+9.0_dp*dom(ib)%u(i,j,kk+1)+ &
                    9.0_dp*dom(ib)%u(i,j,kk)-dom(ib)%u(i,j,kk-1))/16.0_dp)* &
                    ((-dom(ib)%w(i+2,j,kk)+9.0_dp*dom(ib)%w(i+1,j,kk)+ &
                    9.0_dp*dom(ib)%w(i,j,kk)-dom(ib)%w(i-1,j,kk))/16.0_dp)
                          kk=k-2
                          sm2=((-dom(ib)%u(i,j,kk+2)+9.0_dp*dom(ib)%u(i,j,kk+1)+ &
                    9.0_dp*dom(ib)%u(i,j,kk)-dom(ib)%u(i,j,kk-1))/16.0_dp)* &
                    ((-dom(ib)%w(i+2,j,kk)+9.0_dp*dom(ib)%w(i+1,j,kk)+ &
                    9.0_dp*dom(ib)%w(i,j,kk)-dom(ib)%w(i-1,j,kk))/16.0_dp)
                          duwdz=(-sp1+27.0_dp*sn-27.0_dp*sm1+sm2)/(24.0_dp*dom(ib)%dz)

                          else if(differencing==3) then
                          if (dom(ib)%u(i,j,k)>0.0_dp) then
                          dudx = dom(ib)%dphi_dxminus(i,j,k)
                          else if (dom(ib)%u(i,j,k)<0.0_dp) then
                          dudx = dom(ib)%dphi_dxplus(i,j,k)
                          else
                          dudx = 0.0_dp
                          end if

                          vijk=0.25_dp*(dom(ib)%v(i,j,k)+dom(ib)%v(i+1,j,k)+ &
                    dom(ib)%v(i,j-1,k)+dom(ib)%v(i+1,j-1,k))
                          if (vijk>0.0_dp) then
                          dudy = dom(ib)%dphi_dyminus(i,j,k)
                          else if (vijk<0.0_dp) then
                          dudy = dom(ib)%dphi_dyplus(i,j,k)
                          else
                          dudy = 0.0_dp
                          end if

                          wijk=0.25_dp*(dom(ib)%w(i,j,k)+dom(ib)%w(i+1,j,k)+ &
                    dom(ib)%w(i,j,k-1)+dom(ib)%w(i+1,j,k-1))
                          if (wijk>0.0_dp) then
                          dudz = dom(ib)%dphi_dzminus(i,j,k)
                          else if (wijk<0.0_dp) then
                          dudz = dom(ib)%dphi_dzplus(i,j,k)
                          else
                          dudz = 0.0_dp
                          end if

                          du2dx=dom(ib)%u(i,j,k)*dudx
                          duvdy=vijk*dudy
                          duwdz=wijk*dudz

                          end if

                          if(conv_sch==1) then
                          dom(ib)%ustar(i,j,k)=(dom(ib)%u(i,j,k)- &
                    dt*(du2dx+duvdy+duwdz))
                          else if(conv_sch==2) then
                          vel=dom(ib)%uo(i,j,k)
                          dom(ib)%uo(i,j,k)=(du2dx+duvdy+duwdz)
                          dom(ib)%ustar(i,j,k)=(dom(ib)%u(i,j,k)- &
                    dt*(1.5_dp*dom(ib)%uo(i,j,k)-0.5_dp*vel))
                          end if

                      end do
                  end do
              end do
          end do


          if (differencing==3) then
          call HJ_WENO_dx(2)
          call HJ_WENO_dy(2)
          call HJ_WENO_dz(2)
          end if

          do ib=1,nbp

              do k=dom(ib)%ksv,dom(ib)%kev
                  do i=dom(ib)%isv,dom(ib)%iev
                      do j=dom(ib)%jsv,dom(ib)%jev

                          if(differencing==1) then
                          vp12=0.5_dp*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j+1,k))
                          vm12=0.5_dp*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j-1,k))
                          dv2dy=(vp12**2-vm12**2)/dom(ib)%dy
                          duvdx=0.25_dp*( (dom(ib)%u(i,j,k)+dom(ib)%u(i,j+1,k))* &
                    (dom(ib)%v(i,j,k)+dom(ib)%v(i+1,j,k)) - &
                    (dom(ib)%u(i-1,j,k)+dom(ib)%u(i-1,j+1,k))* &
                    (dom(ib)%v(i,j,k)+dom(ib)%v(i-1,j,k)) )/dom(ib)%dx
                          dvwdz=0.25_dp*( (dom(ib)%w(i,j,k)+dom(ib)%w(i,j+1,k))* &
                    (dom(ib)%v(i,j,k)+dom(ib)%v(i,j,k+1)) - &
                    (dom(ib)%w(i,j,k-1)+dom(ib)%w(i,j+1,k-1))* &
                    (dom(ib)%v(i,j,k)+dom(ib)%v(i,j,k-1)) )/dom(ib)%dz

                          else if(differencing==2) then
                          jj=j+2
                          sp2=((-dom(ib)%v(i,jj+1,k)+9.0_dp*dom(ib)%v(i,jj,k)+ &
                    9.0_dp*dom(ib)%v(i,jj-1,k)-dom(ib)%v(i,jj-2,k))/16.0_dp)**2
                          jj=j+1
                          sp1=((-dom(ib)%v(i,jj+1,k)+9.0_dp*dom(ib)%v(i,jj,k)+ &
                    9.0_dp*dom(ib)%v(i,jj-1,k)-dom(ib)%v(i,jj-2,k))/16.0_dp)**2
                          jj=j
                          sn=((-dom(ib)%v(i,jj+1,k)+9.0_dp*dom(ib)%v(i,jj,k)+ &
                    9.0_dp*dom(ib)%v(i,jj-1,k)-dom(ib)%v(i,jj-2,k))/16.0_dp)**2
                          jj=j-1
                          sm1=((-dom(ib)%v(i,jj+1,k)+9.0_dp*dom(ib)%v(i,jj,k)+ &
                    9.0_dp*dom(ib)%v(i,jj-1,k)-dom(ib)%v(i,jj-2,k))/16.0_dp)**2
                          dv2dy=(-sp2+27.0_dp*sp1-27.0_dp*sn+sm1)/(24.0_dp*dom(ib)%dy)

                          ii=i+1
                          sp1=((-dom(ib)%v(ii+2,j,k)+9.0_dp*dom(ib)%v(ii+1,j,k)+ &
                    9.0_dp*dom(ib)%v(ii,j,k)-dom(ib)%v(ii-1,j,k))/16.0_dp)* &
                    ((-dom(ib)%u(ii,j+2,k)+9.0_dp*dom(ib)%u(ii,j+1,k)+ &
                    9.0_dp*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j-1,k))/16.0_dp)
                          ii=i
                          sn=((-dom(ib)%v(ii+2,j,k)+9.0_dp*dom(ib)%v(ii+1,j,k)+ &
                    9.0_dp*dom(ib)%v(ii,j,k)-dom(ib)%v(ii-1,j,k))/16.0_dp)* &
                    ((-dom(ib)%u(ii,j+2,k)+9.0_dp*dom(ib)%u(ii,j+1,k)+ &
                    9.0_dp*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j-1,k))/16.0_dp)
                          ii=i-1
                          sm1=((-dom(ib)%v(ii+2,j,k)+9.0_dp*dom(ib)%v(ii+1,j,k)+ &
                    9.0_dp*dom(ib)%v(ii,j,k)-dom(ib)%v(ii-1,j,k))/16.0_dp)* &
                    ((-dom(ib)%u(ii,j+2,k)+9.0_dp*dom(ib)%u(ii,j+1,k)+ &
                    9.0_dp*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j-1,k))/16.0_dp)
                          ii=i-2
                          sm2=((-dom(ib)%v(ii+2,j,k)+9.0_dp*dom(ib)%v(ii+1,j,k)+ &
                    9.0_dp*dom(ib)%v(ii,j,k)-dom(ib)%v(ii-1,j,k))/16.0_dp)* &
                    ((-dom(ib)%u(ii,j+2,k)+9.0_dp*dom(ib)%u(ii,j+1,k)+ &
                    9.0_dp*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j-1,k))/16.0_dp)
                          duvdx=(-sp1+27.0_dp*sn-27.0_dp*sm1+sm2)/(24.0_dp*dom(ib)%dx)

                          kk=k+1
                          sp1=((-dom(ib)%v(i,j,kk+2)+9.0_dp*dom(ib)%v(i,j,kk+1)+ &
                    9.0_dp*dom(ib)%v(i,j,kk)-dom(ib)%v(i,j,kk-1))/16.0_dp)* &
                    ((-dom(ib)%w(i,j+2,kk)+9.0_dp*dom(ib)%w(i,j+1,kk)+ &
                    9.0_dp*dom(ib)%w(i,j,kk)-dom(ib)%w(i,j-1,kk))/16.0_dp)
                          kk=k
                          sn=((-dom(ib)%v(i,j,kk+2)+9.0_dp*dom(ib)%v(i,j,kk+1)+ &
                    9.0_dp*dom(ib)%v(i,j,kk)-dom(ib)%v(i,j,kk-1))/16.0_dp)* &
                    ((-dom(ib)%w(i,j+2,kk)+9.0_dp*dom(ib)%w(i,j+1,kk)+ &
                    9.0_dp*dom(ib)%w(i,j,kk)-dom(ib)%w(i,j-1,kk))/16.0_dp)
                          kk=k-1
                          sm1=((-dom(ib)%v(i,j,kk+2)+9.0_dp*dom(ib)%v(i,j,kk+1)+ &
                    9.0_dp*dom(ib)%v(i,j,kk)-dom(ib)%v(i,j,kk-1))/16.0_dp)* &
                    ((-dom(ib)%w(i,j+2,kk)+9.0_dp*dom(ib)%w(i,j+1,kk)+ &
                    9.0_dp*dom(ib)%w(i,j,kk)-dom(ib)%w(i,j-1,kk))/16.0_dp)
                          kk=k-2
                          sm2=((-dom(ib)%v(i,j,kk+2)+9.0_dp*dom(ib)%v(i,j,kk+1)+ &
                    9.0_dp*dom(ib)%v(i,j,kk)-dom(ib)%v(i,j,kk-1))/16.0_dp)* &
                    ((-dom(ib)%w(i,j+2,kk)+9.0_dp*dom(ib)%w(i,j+1,kk)+ &
                    9.0_dp*dom(ib)%w(i,j,kk)-dom(ib)%w(i,j-1,kk))/16.0_dp)
                          dvwdz=(-sp1+27.0_dp*sn-27.0_dp*sm1+sm2)/(24.0_dp*dom(ib)%dz)

                          else if(differencing==3) then
                          uijk=0.25_dp*(dom(ib)%u(i,j,k)+dom(ib)%u(i,j+1,k)+ &
                    dom(ib)%u(i-1,j,k)+dom(ib)%u(i-1,j+1,k))
                          if (uijk>0.0_dp) then
                          dvdx = dom(ib)%dphi_dxminus(i,j,k)
                          else if (uijk<0.0_dp) then
                          dvdx = dom(ib)%dphi_dxplus(i,j,k)
                          else
                          dvdx = 0.0_dp
                          end if

                          if (dom(ib)%v(i,j,k)>0.0_dp) then
                          dvdy = dom(ib)%dphi_dyminus(i,j,k)
                          else if (dom(ib)%v(i,j,k)<0.0_dp) then
                          dvdy = dom(ib)%dphi_dyplus(i,j,k)
                          else
                          dvdy = 0.0_dp
                          end if

                          wijk=0.25_dp*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j+1,k)+ &
                    dom(ib)%w(i,j,k-1)+dom(ib)%w(i,j+1,k-1))
                          if (wijk>0.0_dp) then
                          dvdz = dom(ib)%dphi_dzminus(i,j,k)
                          else if (wijk<0.0_dp) then
                          dvdz = dom(ib)%dphi_dzplus(i,j,k)
                          else
                          dvdz = 0.0_dp
                          end if

                          duvdx=uijk*dvdx
                          dv2dy=dom(ib)%v(i,j,k)*dvdy
                          dvwdz=wijk*dvdz

                          end if

                          if(conv_sch==1) then
                          dom(ib)%vstar(i,j,k)=(dom(ib)%v(i,j,k)- &
                    dt*(duvdx+dv2dy+dvwdz))
                          else if(conv_sch==2) then
                          vel=dom(ib)%vo(i,j,k)
                          dom(ib)%vo(i,j,k)=(duvdx+dv2dy+dvwdz)
                          dom(ib)%vstar(i,j,k)=(dom(ib)%v(i,j,k)- &
                    dt*(1.5_dp*dom(ib)%vo(i,j,k)-0.5_dp*vel))
                          end if

                      end do
                  end do
              end do
          end do

          if (differencing==3) then
          call HJ_WENO_dx(3)
          call HJ_WENO_dy(3)
          call HJ_WENO_dz(3)
          end if

          do ib=1,nbp

              do k=dom(ib)%ksw,dom(ib)%kew
                  do i=dom(ib)%isw,dom(ib)%iew
                      do j=dom(ib)%jsw,dom(ib)%jew

                          if(differencing==1) then
                          wp12=0.5_dp*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k+1))
                          wm12=0.5_dp*(dom(ib)%w(i,j,k)+dom(ib)%w(i,j,k-1))
                          dw2dz=(wp12**2-wm12**2)/dom(ib)%dz
                          duwdx=0.25_dp*( (dom(ib)%u(i,j,k)+dom(ib)%u(i,j,k+1))* &
                    (dom(ib)%w(i,j,k)+dom(ib)%w(i+1,j,k)) - &
                    (dom(ib)%u(i-1,j,k)+dom(ib)%u(i-1,j,k+1))* &
                    (dom(ib)%w(i,j,k)+dom(ib)%w(i-1,j,k)) )/dom(ib)%dx
                          dvwdy=0.25_dp*( (dom(ib)%v(i,j,k)+dom(ib)%v(i,j,k+1))* &
                    (dom(ib)%w(i,j,k)+dom(ib)%w(i,j+1,k)) - &
                    (dom(ib)%v(i,j-1,k)+dom(ib)%v(i,j-1,k+1))* &
                    (dom(ib)%w(i,j,k)+dom(ib)%w(i,j-1,k)) )/dom(ib)%dy

                          else if(differencing==2) then
                          kk=k+2
                          sp2=((-dom(ib)%w(i,j,kk+1)+9.0_dp*dom(ib)%w(i,j,kk)+ &
                    9.0_dp*dom(ib)%w(i,j,kk-1)-dom(ib)%w(i,j,kk-2))/16.0_dp)**2
                          kk=k+1
                          sp1=((-dom(ib)%w(i,j,kk+1)+9.0_dp*dom(ib)%w(i,j,kk)+ &
                    9.0_dp*dom(ib)%w(i,j,kk-1)-dom(ib)%w(i,j,kk-2))/16.0_dp)**2
                          kk=k
                          sn=((-dom(ib)%w(i,j,kk+1)+9.0_dp*dom(ib)%w(i,j,kk)+ &
                    9.0_dp*dom(ib)%w(i,j,kk-1)-dom(ib)%w(i,j,kk-2))/16.0_dp)**2
                          kk=k-1
                          sm1=((-dom(ib)%w(i,j,kk+1)+9.0_dp*dom(ib)%w(i,j,kk)+ &
                    9.0_dp*dom(ib)%w(i,j,kk-1)-dom(ib)%w(i,j,kk-2))/16.0_dp)**2
                          dw2dz=(-sp2+27.0_dp*sp1-27.0_dp*sn+sm1)/(24.0_dp*dom(ib)%dz)

                          ii=i+1
                          sp1=((-dom(ib)%w(ii+2,j,k)+9.0_dp*dom(ib)%w(ii+1,j,k)+ &
                    9.0_dp*dom(ib)%w(ii,j,k)-dom(ib)%w(ii-1,j,k))/16.0_dp)* &
                    ((-dom(ib)%u(ii,j,k+2)+9.0_dp*dom(ib)%u(ii,j,k+1)+ &
                    9.0_dp*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j,k-1))/16.0_dp)
                          ii=i
                          sn=((-dom(ib)%w(ii+2,j,k)+9.0_dp*dom(ib)%w(ii+1,j,k)+ &
                    9.0_dp*dom(ib)%w(ii,j,k)-dom(ib)%w(ii-1,j,k))/16.0_dp)* &
                    ((-dom(ib)%u(ii,j,k+2)+9.0_dp*dom(ib)%u(ii,j,k+1)+ &
                    9.0_dp*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j,k-1))/16.0_dp)
                          ii=i-1
                          sm1=((-dom(ib)%w(ii+2,j,k)+9.0_dp*dom(ib)%w(ii+1,j,k)+ &
                    9.0_dp*dom(ib)%w(ii,j,k)-dom(ib)%w(ii-1,j,k))/16.0_dp)* &
                    ((-dom(ib)%u(ii,j,k+2)+9.0_dp*dom(ib)%u(ii,j,k+1)+ &
                    9.0_dp*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j,k-1))/16.0_dp)
                          ii=i-2
                          sm2=((-dom(ib)%w(ii+2,j,k)+9.0_dp*dom(ib)%w(ii+1,j,k)+ &
                    9.0_dp*dom(ib)%w(ii,j,k)-dom(ib)%w(ii-1,j,k))/16.0_dp)* &
                    ((-dom(ib)%u(ii,j,k+2)+9.0_dp*dom(ib)%u(ii,j,k+1)+ &
                    9.0_dp*dom(ib)%u(ii,j,k)-dom(ib)%u(ii,j,k-1))/16.0_dp)
                          duwdx=(-sp1+27.0_dp*sn-27.0_dp*sm1+sm2)/(24.0_dp*dom(ib)%dx)

                          jj=j+1
                          sp1=((-dom(ib)%w(i,jj+2,k)+9.0_dp*dom(ib)%w(i,jj+1,k)+ &
                    9.0_dp*dom(ib)%w(i,jj,k)-dom(ib)%w(i,jj-1,k))/16.0_dp)* &
                    ((-dom(ib)%v(i,jj,k+2)+9.0_dp*dom(ib)%v(i,jj,k+1)+ &
                    9.0_dp*dom(ib)%v(i,jj,k)-dom(ib)%v(i,jj,k-1))/16.0_dp)
                          jj=j
                          sn=((-dom(ib)%w(i,jj+2,k)+9.0_dp*dom(ib)%w(i,jj+1,k)+ &
                    9.0_dp*dom(ib)%w(i,jj,k)-dom(ib)%w(i,jj-1,k))/16.0_dp)* &
                    ((-dom(ib)%v(i,jj,k+2)+9.0_dp*dom(ib)%v(i,jj,k+1)+ &
                    9.0_dp*dom(ib)%v(i,jj,k)-dom(ib)%v(i,jj,k-1))/16.0_dp)
                          jj=j-1
                          sm1=((-dom(ib)%w(i,jj+2,k)+9.0_dp*dom(ib)%w(i,jj+1,k)+ &
                    9.0_dp*dom(ib)%w(i,jj,k)-dom(ib)%w(i,jj-1,k))/16.0_dp)* &
                    ((-dom(ib)%v(i,jj,k+2)+9.0_dp*dom(ib)%v(i,jj,k+1)+ &
                    9.0_dp*dom(ib)%v(i,jj,k)-dom(ib)%v(i,jj,k-1))/16.0_dp)
                          jj=j-2
                          sm2=((-dom(ib)%w(i,jj+2,k)+9.0_dp*dom(ib)%w(i,jj+1,k)+ &
                    9.0_dp*dom(ib)%w(i,jj,k)-dom(ib)%w(i,jj-1,k))/16.0_dp)* &
                    ((-dom(ib)%v(i,jj,k+2)+9.0_dp*dom(ib)%v(i,jj,k+1)+ &
                    9.0_dp*dom(ib)%v(i,jj,k)-dom(ib)%v(i,jj,k-1))/16.0_dp)
                          dvwdy=(-sp1+27.0_dp*sn-27.0_dp*sm1+sm2)/(24.0_dp*dom(ib)%dy)

                          else if(differencing==3) then
                          uijk=0.25_dp*(dom(ib)%u(i,j,k)+dom(ib)%u(i,j,k+1)+ &
                    dom(ib)%u(i-1,j,k)+dom(ib)%u(i-1,j,k+1))
                          if (uijk>0.0_dp) then
                          dwdx = dom(ib)%dphi_dxminus(i,j,k)
                          else if (uijk<0.0_dp) then
                          dwdx = dom(ib)%dphi_dxplus(i,j,k)
                          else
                          dwdx = 0.0_dp
                          end if

                          vijk=0.25_dp*(dom(ib)%v(i,j,k)+dom(ib)%v(i,j,k+1)+ &
                    dom(ib)%v(i,j-1,k)+dom(ib)%v(i,j-1,k+1))
                          if (vijk>0.0_dp) then
                          dwdy = dom(ib)%dphi_dyminus(i,j,k)
                          else if (vijk<0.0_dp) then
                          dwdy = dom(ib)%dphi_dyplus(i,j,k)
                          else
                          dwdy = 0.0_dp
                          end if

                          if (dom(ib)%w(i,j,k)>0.0_dp) then
                          dwdz = dom(ib)%dphi_dzminus(i,j,k)
                          else if (dom(ib)%w(i,j,k)<0.0_dp) then
                          dwdz = dom(ib)%dphi_dzplus(i,j,k)
                          else
                          dwdz = 0.0_dp
                          end if

                          duwdx=uijk*dwdx
                          dvwdy=vijk*dwdy
                          dw2dz=dom(ib)%w(i,j,k)*dwdz

                          end if

                          if(conv_sch==1) then
                          dom(ib)%wstar(i,j,k)=(dom(ib)%w(i,j,k)- &
                    dt*(duwdx+dvwdy+dw2dz))
                          else if(conv_sch==2) then
                          vel=dom(ib)%wo(i,j,k)
                          dom(ib)%wo(i,j,k)=(duwdx+dvwdy+dw2dz)
                          dom(ib)%wstar(i,j,k)=(dom(ib)%w(i,j,k)- &
                    dt*(1.5_dp*dom(ib)%wo(i,j,k)-0.5_dp*vel))
                          end if


                      end do
                  end do
              end do

          end do

          return
      end subroutine convection
!##########################################################################
