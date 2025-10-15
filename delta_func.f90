!=======================================================================
!           Pablo Ouro Barba
!           Cardiff 2013-2014
!=======================================================================
!######################################################################
      real function phi_r1smth(r)
!######################################################################
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          real, intent(in) :: r
          real :: PI,abr
          PI = 4.D0*DATAN(1.D0)
          abr=SQRT(r*r)
          if (abr>=1.5_dp) then
          phi_r1smth = 0.0_dp
          else if ((abr<1.5_dp).and.(abr>=0.5_dp)) then
          phi_r1smth = 9.0_dp/8.0_dp-3.0_dp*abr/2+abr**2/2
          else if ((abr<0.5_dp).and.(abr>=0.0_dp)) then
          phi_r1smth = 3.0_dp/4.0_dp-abr**2
          end if
          return
      end function
!######################################################################
      real function phi_r2smth(r)
!######################################################################
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          real, intent(in) :: r
          real :: PI
          PI = 4.D0*DATAN(1.D0)
          if (r<=-2.5_dp) then
          phi_r2smth = 0.0_dp
          else if ((r>=-2.5_dp).and.(r<=-1.5_dp)) then
          phi_r2smth= -1.0_dp/8.0_dp/PI*(-5.0_dp*PI-2.0_dp*PI*r+4.0_dp*sin(PI/4.0_dp*(-2.0_dp*r-1.0_dp)))
          else if ((r>=-1.5_dp).and.(r<=0.0_dp)) then
          phi_r2smth = 1.0_dp/4.0_dp/PI*(PI+2.0_dp*sin(PI/4.0_dp*(-2.0_dp*r+1.0_dp)) &
                             -2.0_dp*sin(PI/4.0_dp*(-2.0_dp*r-1.0_dp)))
          else if ((r>=0.0_dp).and.(r<=1.5_dp)) then
          phi_r2smth = 1.0_dp/4.0_dp/PI*(PI+2.0_dp*sin(PI/4.0_dp*(2.0_dp*r+1.0_dp)) &
                             -2.0_dp*sin(PI/4.0_dp*(2.0_dp*r-1.0_dp)))
          else if ((r>=1.5_dp).and.(r<=2.5_dp)) then
          phi_r2smth= -1.0_dp/8.0_dp/PI*(-5.0_dp*PI+2.0_dp*PI*r+4.0_dp*sin(PI/4.0_dp*(2.0_dp*r-1.0_dp)))
          else if (r>=2.5_dp) then
          phi_r2smth = 0.0_dp
          end if

          return
      end function
!######################################################################
      real function phi_r3(r)
!######################################################################
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          real, intent(in) :: r
          if (r<=-1.5_dp) then
          phi_r3 = 0.0_dp
          else if ((r>=-1.5_dp).and.(r<=-0.5_dp)) then
          phi_r3 = 1.0_dp/6.0_dp*(5.0_dp+3.0_dp*r-sqrt(-3.0_dp*(1.0_dp+r)**2+1.0_dp))
          else if ((r>=-0.5_dp).and.(r<=0.0_dp)) then
          phi_r3 = 1.0_dp/3.0_dp*(1.0_dp+sqrt(-3.0_dp*r**2+1.0_dp))
          else if ((r>=0.0_dp).and.(r<=0.5_dp)) then
          phi_r3 = 1.0_dp/3.0_dp*(1.0_dp+sqrt(-3.0_dp*r**2+1.0_dp))
          else if ((r>=0.5_dp).and.(r<=1.5_dp)) then
          phi_r3 = 1.0_dp/6.0_dp*(5.0_dp-3.0_dp*r-sqrt(-3.0_dp*(1.0_dp-r)**2+1.0_dp))
          else if (r>=1.5_dp) then
          phi_r3 = 0.0_dp
          end if

          return
      end function
!######################################################################
      real function phi_r3smth2(r)
!######################################################################
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          real, intent(in) :: r
          real :: PI
          PI = 4.D0*DATAN(1.D0)
          if (r<=-2.0_dp) then
          phi_r3smth2 = 0.0_dp
          else if ((r>=-2.0_dp).and.(r<=-1.0_dp)) then
          phi_r3smth2 = 55.0_dp/48.0_dp - sqrt(3.0_dp)*pi/108.0_dp + 13.0_dp*r/12.0_dp &
    + r**2/4.0_dp + (-2.0_dp*r-3.0_dp)/48.0_dp*sqrt(-12.0_dp*r**2-36.0_dp*r-23.0_dp) &
    + sqrt(3.0_dp)/36.0_dp*ASIN(sqrt(3.0_dp)/2.0_dp*(-2.0_dp*r-3.0_dp))
          else if ((r>=-1.0_dp).and.(r<=0.0_dp)) then
          phi_r3smth2 = 17.0_dp/48.0_dp + sqrt(3.0_dp)*pi/108.0_dp - r/4.0_dp &
    - r**2/4.0_dp + (2.0_dp*r+1.0_dp)/16.0_dp*sqrt(-12.0_dp*r**2-12.0_dp*r+1.0_dp) &
    - sqrt(3.0_dp)/12.0_dp*ASIN(sqrt(3.0_dp)/2.0_dp*(-2.0_dp*r-1.0_dp))
          else if ((r>=0.0_dp).and.(r<=1.0_dp)) then
          phi_r3smth2 = 17.0_dp/48.0_dp + sqrt(3.0_dp)*pi/108.0_dp + r/4.0_dp &
    - r**2/4.0_dp + (-2.0_dp*r+1.0_dp)/16.0_dp*sqrt(-12.0_dp*r**2+12.0_dp*r+1.0_dp) &
    - sqrt(3.0_dp)/12.0_dp*ASIN(sqrt(3.0_dp)/2.0_dp*(2.0_dp*r-1.0_dp))
          else if ((r>=1.0_dp).and.(r<=2.0_dp)) then
          phi_r3smth2 = 55.0_dp/48.0_dp - sqrt(3.0_dp)*pi/108.0_dp - 13.0_dp*r/12.0_dp &
    + r**2/4.0_dp + (2.0_dp*r-3.0_dp)/48.0_dp*sqrt(-12.0_dp*r**2+36.0_dp*r-23.0_dp) &
    + sqrt(3.0_dp)/36.0_dp*ASIN(sqrt(3.0_dp)/2.0_dp*(2.0_dp*r-3.0_dp))
          else if (r>=2.0_dp) then
          phi_r3smth2 = 0.0_dp
          end if

          return
      end function
!######################################################################
      real function phi_r3smth(r)
!######################################################################
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          double precision, intent(in) :: r
!      real :: PI
          real :: scal1,scal2,scal3,scal4,scal5,scal6

          scal1 = 1.095450017660160615261_dp  ! 55.0_dp/48.0_dp - sqrt(3.0_dp)*pi/108.0_dp
          scal2 = 1.083333333333333333333_dp  ! 13.0_dp/12.0_dp
          scal3 = 0.4045499823398393847387_dp  ! 17.0_dp/48.0_dp + sqrt(3.0_dp)*pi/108.0_dp
          scal4 = 0.0481125224324688137091_dp  ! sqrt(3.0_dp)/36.0_dp
          scal5 = 0.1443375672974064411273_dp  ! sqrt(3.0_dp)/12.0_dp
          scal6 = 0.8660254037844386467637_dp  ! sqrt(3.0_dp)/2.0_dp

!       PI = 4.D0*DATAN(1.D0)

          if (r<=-2.0_dp) then
          phi_r3smth = 0.0_dp
          else if ((r>=-2.0_dp).and.(r<=-1.0_dp)) then
          phi_r3smth = scal1 + scal2*r + scal4*ASIN(scal6*(-2.0_dp*r-3.0_dp)) &
    + 0.25_dp*r**2 + (-2.0_dp*r-3.0_dp)/48.0_dp*sqrt(-12.0_dp*r**2-36.0_dp*r-23.0_dp)
          else if ((r>=-1.0_dp).and.(r<=0.0_dp)) then
          phi_r3smth = scal3 - r/4.0_dp - scal5*ASIN(scal6*(-2.0_dp*r-1.0_dp)) &
    - 0.25_dp*r**2 + (2.0_dp*r+1.0_dp)/16.0_dp*sqrt(-12.0_dp*r**2-12.0_dp*r+1.0_dp)
          else if ((r>=0.0_dp).and.(r<=1.0_dp)) then

          phi_r3smth = scal3 + r/4.0_dp - scal5*ASIN(scal6*(2.0_dp*r-1.0_dp)) &
    - 0.25_dp*r**2 + (-2.0_dp*r+1.0_dp)/16.0_dp*sqrt(-12.0_dp*r**2+12.0_dp*r+1.0_dp)
          else if ((r>=1.0_dp).and.(r<=2.0_dp)) then
          phi_r3smth = scal1 - scal2*r + scal4*ASIN(scal6*(2.0_dp*r-3.0_dp)) &
    + 0.25_dp*r**2 + (2.0_dp*r-3.0_dp)/48.0_dp*sqrt(-12.0_dp*r**2+36.0_dp*r-23.0_dp)
          else if (r>=2.0_dp) then
          phi_r3smth = 0.0_dp
          end if

          return
      end function
!######################################################################
      real function phi_r4(r)
!######################################################################
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          real, intent(in) :: r
          if (r<=-2.0_dp) then
          phi_r4 = 0.0_dp
          else if ((r>=-2.0_dp).and.(r<=-1.0_dp)) then
          phi_r4 = 1.0_dp/8.0_dp*(5.0_dp+2.0_dp*r-sqrt(-7.0_dp-12.0_dp*r-4.0_dp*r**2))
          else if ((r>=-1.0_dp).and.(r<=0.0_dp)) then
          phi_r4 = 1.0_dp/8.0_dp*(3.0_dp+2.0_dp*r+sqrt(1.0_dp-4.0_dp*r-4.0_dp*r**2))
          else if ((r>=0.0_dp).and.(r<=1.0_dp)) then
          phi_r4 = 1.0_dp/8.0_dp*(3.0_dp-2.0_dp*r+sqrt(1.0_dp+4.0_dp*r-4.0_dp*r**2))
          else if ((r>=1.0_dp).and.(r<=2.0_dp)) then
          phi_r4 = 1.0_dp/8.0_dp*(5.0_dp-2.0_dp*r-sqrt(-7.0_dp+12.0_dp*r-4.0_dp*r**2))
          else if (r>=2.0_dp) then
          phi_r4 = 0.0_dp
          end if

          return
      end function
!######################################################################
      real function phi_r4smth(r)
!######################################################################
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          real, intent(in) :: r
          real :: PI,ar
          PI = 4.D0*DATAN(1.D0)
          ar=abs(r)
          if ((ar>=0.0_dp).and.(ar<=0.5_dp)) then
          phi_r4smth = 3.0_dp/8.0_dp+ PI/32.0_dp - r*r/4.0_dp
          else if ((ar>=0.5_dp).and.(ar<=1.5_dp)) then
          phi_r4smth = 1.0_dp/4.0_dp + (1.0_dp-ar)/8.0_dp*SQRT(-2.0_dp+8.0_dp*ar -4.0_dp*r*r) &
     -1.0_dp/8.0_dp*ASIN(sqrt(2.0_dp)*(ar-1.0_dp))
          else if ((ar>=1.5_dp).and.(ar<=2.5_dp)) then
          phi_r4smth = 17.0_dp/16.0_dp - PI/64.0_dp - 3.0_dp*ar/4.0_dp+ r*r/8.0_dp &
     + (ar-2.0_dp)/16.0_dp * SQRT(-14.0_dp+ 16.0_dp*ar - 4.0_dp*r*r) &
     + 1.0_dp/16.0_dp*ASIN(sqrt(2.0_dp)*(ar-2.0_dp))
          else if (ar>=2.5_dp) then
          phi_r4smth = 0.0_dp
          end if
          return
      end function


!######################################################################
      real function cd2_0(r)
!######################################################################
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          real, intent(in) :: r
          if (r<=-1.5_dp) then
          cd2_0 = 0.0_dp
          else if ((r>=-1.5_dp).and.(r<=-0.5_dp)) then
          cd2_0 = 0.5_dp*r+0.75_dp
          else if ((r>=-0.5_dp).and.(r<=0.5_dp)) then
          cd2_0 = 0.5_dp
          else if ((r>=0.5_dp).and.(r<=1.5_dp)) then
          cd2_0 = -0.5_dp*r+0.75_dp
          else if (r>=1.5_dp) then
          cd2_0 = 0.0_dp
          end if

          return
      end function

!######################################################################
      real function dcd2_0(r)
!######################################################################
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          real, intent(in) :: r
          if (r<=-1.5_dp) then
          dcd2_0 = 0.0_dp
          else if ((r>=-1.5_dp).and.(r<=-0.5_dp)) then
          dcd2_0 = 0.5_dp
          else if ((r>=-0.5_dp).and.(r<=0.5_dp)) then
          dcd2_0 = 0.0_dp
          else if ((r>=0.5_dp).and.(r<=1.5_dp)) then
          dcd2_0 = -0.5_dp
          else if (r>=1.5_dp) then
          dcd2_0 = 0.0_dp
          end if

          return
      end function

!######################################################################
      real function cd2_1(r)
!######################################################################
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          real, intent(in) :: r
          if (r<=-2.0_dp) then
          cd2_1 = 0.0_dp
          else if ((r>=-2.0_dp).and.(r<=-1.0_dp)) then
          cd2_1 = 0.25_dp*r**2+r+1.0_dp
          else if ((r>=-1.0_dp).and.(r<=1.0_dp)) then
          cd2_1 = -0.25_dp*r**2+0.5_dp
          else if ((r>=1.0_dp).and.(r<=2.0_dp)) then
          cd2_1 = 0.25_dp*r**2-r+1.0_dp
          else if (r>=2.0_dp) then
          cd2_1 = 0.0_dp
          end if

          return
      end function


!######################################################################
      real function dcd2_1(r)
!######################################################################
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          real, intent(in) :: r
          if (r<=-2.0_dp) then
          dcd2_1 = 0.0_dp
          else if ((r>=-2.0_dp).and.(r<=-1.0_dp)) then
          dcd2_1 = 0.5_dp*r+1.0_dp
          else if ((r>=-1.0_dp).and.(r<=1.0_dp)) then
          dcd2_1 = -0.5_dp*r
          else if ((r>=1.0_dp).and.(r<=2.0_dp)) then
          dcd2_1 = 0.5_dp*r-1.0_dp
          else if (r>=2.0_dp) then
          dcd2_1 = 0.0_dp
          end if

          return
      end function
!######################################################################
      double precision function dh(dx,dy,dz,xij,yij,zij,Xl,Yl,Zl,order)
!######################################################################
! ..... The salient properties of the kernels DIRAC DELTA,dh are the following:
! * dh is a continuously diﬀerentiable function and therefore yields
!   a smoother transfer than e.g. linear interpolation.
! * Interpolation using the kernels dh is second-order accurate
!   for smooth ﬁelds (Uhlmann(2005) cf. Section 5.1_dp.1_dp).
! * The support of the regularized delta function is small, which makes
!   the evaluation of the sums in Eq. (9) relatively cheap. In particular,
!   we use the expression for dh deﬁned by Roma et al., involving only
!   three grid points in each coordinate direction.
! * (ORDER=3) A. Roma, C. Peskin, M. Berger, An adaptive version of the
!   immersed boundary method, J. Comput. Phys. 153 (1999)
! * (ORDER=4) C. Peskin, The immersed boundary method,
!   Acta Numerica 11 (2002) 1–39
!
          implicit none

          real,    intent(in) :: dx,dy,dz,xij,yij,zij,Xl,Yl,Zl
          integer, intent(in) :: order
          real(kind=8) ::phi_r2smth,phi_r3smth,phi_r4,phi_r3,phi_r1smth,phi_r4smth
          real(kind=8) :: cd2_0,cd2_1
          select case (order)

            case (1)
              dh =  phi_r1smth((xij-Xl)/dx) &
         * phi_r1smth((yij-Yl)/dy) * phi_r1smth((zij-Zl)/dz)
            case (2)
              dh =   phi_r2smth((xij-Xl)/dx) &
         * phi_r2smth((yij-Yl)/dy) * phi_r2smth((zij-Zl)/dz)
            case (3)
              dh =   phi_r3smth((xij-Xl)/dx) &
          * phi_r3smth((yij-Yl)/dy) * phi_r3smth((zij-Zl)/dz)
            case (4)
              dh =   phi_r4smth((xij-Xl)/dx) &
          * phi_r4smth((yij-Yl)/dy) * phi_r4smth((zij-Zl)/dz)
            case (5)
              dh =   phi_r3((xij-Xl)/dx) &
          * phi_r3((yij-Yl)/dy) * phi_r3((zij-Zl)/dz)
            case (6)
              dh =   phi_r4((xij-Xl)/dx) &
           * phi_r4((yij-Yl)/dy) * phi_r4((zij-Zl)/dz)
            case (7)
              dh =   cd2_0((xij-Xl)/dx) &
           * cd2_0((yij-Yl)/dy) * cd2_0((zij-Zl)/dz)
            case (8)
              dh =   cd2_1((xij-Xl)/dx) &
           * cd2_1((yij-Yl)/dy) * cd2_1((zij-Zl)/dz)

            case default

              print*, '===ERROR==='
              print*, ' order of delta function is not selected '
              stop

          end Select

          return
      end function



!######################################################################
      real function dphi_r3(r)
!######################################################################
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none
          real, intent(in) :: r
          if (r<=-1.5_dp) then
          dphi_r3 = 0.0_dp
          else if ((r>=-1.5_dp).and.(r<=-0.5_dp)) then
          dphi_r3 = 1.0_dp/6.0_dp*(3.0_dp-0.5_dp*(-3.0_dp*(1.0_dp+r)**2+1.0_dp)**(-0.5_dp) &
    *(-6.0_dp*(1.0_dp+r)))
          else if ((r>=-0.5_dp).and.(r<=0.5_dp)) then
          dphi_r3 = 1.0_dp/6.0_dp*(-3.0_dp*r**2+1.0_dp)**(-0.5_dp)*(-6.0_dp*r)
          else if ((r>=0.5_dp).and.(r<=1.5_dp)) then
          dphi_r3 = 1.0_dp/6.0_dp*(-3.0_dp-0.5_dp*(-3.0_dp*(1.0_dp-r)**2+1.0_dp)**(-0.5_dp) &
    *(6.0_dp*(1.0_dp-r)))
          else if (r>=1.5_dp) then
          dphi_r3 = 0.0_dp
          end if

          return
      end function




!######################################################################
      double precision &
function ddh(dx,dy,dz,xij,yij,zij,Xl,Yl,Zl,order,dir)
!######################################################################
          use, intrinsic :: iso_fortran_env, only: dp => real64
          implicit none

          real,    intent(in) :: dx,dy,dz,xij,yij,zij,Xl,Yl,Zl
          integer, intent(in) :: order,dir
          real(kind=8) :: phi_r3,dphi_r3,cd2_0,dcd2_0,cd2_1,dcd2_1

          select case (order)

            case (5)
              select case (dir)

                case (1)
                  ddh = dphi_r3((xij-Xl)/dx) &
              * phi_r3((yij-Yl)/dy) * phi_r3((zij-Zl)/dz) &
              * (-1.0_dp/dx)

                case (2)
                  ddh = phi_r3((xij-Xl)/dx) &
              * dphi_r3((yij-Yl)/dy) * phi_r3((zij-Zl)/dz) &
              * (-1.0_dp/dy)

                case (3)
                  ddh = phi_r3((xij-Xl)/dx) &
              * phi_r3((yij-Yl)/dy) * dphi_r3((zij-Zl)/dz) &
              * (-1.0_dp/dz)

              end select




            case (7)
              select case (dir)

                case (1)
                  ddh = dcd2_0((xij-Xl)/dx) &
              * cd2_0((yij-Yl)/dy) * cd2_0((zij-Zl)/dz) &
              * (-1.0_dp/dx)

                case (2)
                  ddh = cd2_0((xij-Xl)/dx) &
              * dcd2_0((yij-Yl)/dy) * cd2_0((zij-Zl)/dz) &
              * (-1.0_dp/dy)

                case (3)
                  ddh = cd2_0((xij-Xl)/dx) &
              * cd2_0((yij-Yl)/dy) * dcd2_0((zij-Zl)/dz) &
              * (-1.0_dp/dz)

              end select


          end select

          return
      end function
