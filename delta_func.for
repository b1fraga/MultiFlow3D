!=======================================================================
!			Pablo Ouro Barba
!			Cardiff 2013-2014
!=======================================================================
!######################################################################
      real function phi_r1smth(r)
!######################################################################
      implicit none
      real, intent(in) :: r
      real :: PI,abr
       PI = 4.D0*DATAN(1.D0)
       abr=SQRT(r*r)
      if (abr.ge.1.5) then
        phi_r1smth = 0.0
      else if ((abr.lt.1.5).and.(abr.ge.0.5)) then
        phi_r1smth = 9./8.-3.*abr/2+abr**2/2
      else if ((abr.lt.0.5).and.(abr.ge.0.0)) then
        phi_r1smth = 3./4.-abr**2
      end if
      return
      end function
!######################################################################
      real function phi_r2smth(r)
!######################################################################
      implicit none
      real, intent(in) :: r
      real :: PI
       PI = 4.D0*DATAN(1.D0)
      if (r.le.-2.5) then
        phi_r2smth = 0.0
      else if ((r.ge.-2.5).and.(r.le.-1.5)) then
        phi_r2smth= -1./8./PI*(-5.*PI-2.*PI*r+4.*sin(PI/4.*(-2.*r-1.)))
      else if ((r.ge.-1.5).and.(r.le.0.0)) then
        phi_r2smth = 1./4./PI*(PI+2.*sin(PI/4.*(-2.*r+1.))
     &                           -2.*sin(PI/4.*(-2.*r-1.)))
      else if ((r.ge.0.0).and.(r.le.1.5)) then
        phi_r2smth = 1./4./PI*(PI+2.*sin(PI/4.*(2.*r+1.))
     &                           -2.*sin(PI/4.*(2.*r-1.)))
      else if ((r.ge.1.5).and.(r.le.2.5)) then
        phi_r2smth= -1./8./PI*(-5.*PI+2.*PI*r+4.*sin(PI/4.*(2.*r-1.)))
      else if (r.ge.2.5) then
        phi_r2smth = 0.0
      end if

      return
      end function
!######################################################################
      real function phi_r3(r)
!######################################################################
      implicit none
      real, intent(in) :: r
      if (r.le.-1.5) then
        phi_r3 = 0.0
      else if ((r.ge.-1.5).and.(r.le.-0.5)) then
        phi_r3 = 1.0/6.0*(5.0+3.0*r-sqrt(-3.0*(1.0+r)**2+1.0))
      else if ((r.ge.-0.5).and.(r.le.0.0)) then
        phi_r3 = 1.0/3.0*(1.0+sqrt(-3.0*r**2+1.0))
      else if ((r.ge.0.0).and.(r.le.0.5)) then
        phi_r3 = 1.0/3.0*(1.0+sqrt(-3.0*r**2+1.0))
      else if ((r.ge.0.5).and.(r.le.1.5)) then
        phi_r3 = 1.0/6.0*(5.0-3.0*r-sqrt(-3.0*(1.0-r)**2+1.0))
      else if (r.ge.1.5) then
        phi_r3 = 0.0
      end if

      return
      end function
!######################################################################
      real function phi_r3smth2(r)
!######################################################################
      implicit none
      real, intent(in) :: r
      real :: PI
       PI = 4.D0*DATAN(1.D0)
      if (r.le.-2.0) then
        phi_r3smth2 = 0.0
      else if ((r.ge.-2.0).and.(r.le.-1.0)) then
        phi_r3smth2 = 55.0/48.0 - sqrt(3.0)*pi/108.0 + 13.0*r/12.0
     &+ r**2/4.0 + (-2.0*r-3.0)/48.0*sqrt(-12.0*r**2-36.0*r-23.0)
     &+ sqrt(3.0)/36.0*ASIN(sqrt(3.0)/2.0*(-2.0*r-3.0))
      else if ((r.ge.-1.0).and.(r.le.0.0)) then
        phi_r3smth2 = 17.0/48.0 + sqrt(3.0)*pi/108.0 - r/4.0
     &- r**2/4.0 + (2.0*r+1.0)/16.0*sqrt(-12.0*r**2-12.0*r+1.0)
     &- sqrt(3.0)/12.0*ASIN(sqrt(3.0)/2.0*(-2.0*r-1.0))
      else if ((r.ge.0.0).and.(r.le.1.0)) then
        phi_r3smth2 = 17.0/48.0 + sqrt(3.0)*pi/108.0 + r/4.0
     &- r**2/4.0 + (-2.0*r+1.0)/16.0*sqrt(-12.0*r**2+12.0*r+1.0)
     &- sqrt(3.0)/12.0*ASIN(sqrt(3.0)/2.0*(2.0*r-1.0))
      else if ((r.ge.1.0).and.(r.le.2.0)) then
        phi_r3smth2 = 55.0/48.0 - sqrt(3.0)*pi/108.0 - 13.0*r/12.0
     &+ r**2/4.0 + (2.0*r-3.0)/48.0*sqrt(-12.0*r**2+36.0*r-23.0)
     &+ sqrt(3.0)/36.0*ASIN(sqrt(3.0)/2.0*(2.0*r-3.0))
      else if (r.ge.2.0) then
        phi_r3smth2 = 0.0
      end if

      return
      end function
!######################################################################
      real function phi_r3smth(r)
!######################################################################
      implicit none
      double precision, intent(in) :: r
!      real :: PI
      real :: scal1,scal2,scal3,scal4,scal5,scal6

      scal1 = 1.095450017660160615261  ! 55.0/48.0 - sqrt(3.0)*pi/108.0
      scal2 = 1.083333333333333333333  ! 13.0/12.0
      scal3 = 0.4045499823398393847387 ! 17.0/48.0 + sqrt(3.0)*pi/108.0
      scal4 = 0.0481125224324688137091 ! sqrt(3.0)/36.0
      scal5 = 0.1443375672974064411273 ! sqrt(3.0)/12.0
      scal6 = 0.8660254037844386467637 ! sqrt(3.0)/2.0

!       PI = 4.D0*DATAN(1.D0)

      if (r.le.-2.0) then
        phi_r3smth = 0.0
      else if ((r.ge.-2.0).and.(r.le.-1.0)) then
        phi_r3smth = scal1 + scal2*r + scal4*ASIN(scal6*(-2.0*r-3.0))
     &+ 0.25*r**2 + (-2.0*r-3.0)/48.0*sqrt(-12.0*r**2-36.0*r-23.0)
      else if ((r.ge.-1.0).and.(r.le.0.0)) then
        phi_r3smth = scal3 - r/4.0 - scal5*ASIN(scal6*(-2.0*r-1.0))
     &- 0.25*r**2 + (2.0*r+1.0)/16.0*sqrt(-12.0*r**2-12.0*r+1.0)
      else if ((r.ge.0.0).and.(r.le.1.0)) then

        phi_r3smth = scal3 + r/4.0 - scal5*ASIN(scal6*(2.0*r-1.0))
     &- 0.25*r**2 + (-2.0*r+1.0)/16.0*sqrt(-12.0*r**2+12.0*r+1.0)
      else if ((r.ge.1.0).and.(r.le.2.0)) then
        phi_r3smth = scal1 - scal2*r + scal4*ASIN(scal6*(2.0*r-3.0))
     &+ 0.25*r**2 + (2.0*r-3.0)/48.0*sqrt(-12.0*r**2+36.0*r-23.0)
      else if (r.ge.2.0) then
        phi_r3smth = 0.0
      end if

      return
      end function
!######################################################################
      real function phi_r4(r)
!######################################################################
      implicit none
      real, intent(in) :: r
      if (r.le.-2.0) then
        phi_r4 = 0.0
      else if ((r.ge.-2.0).and.(r.le.-1.0)) then
        phi_r4 = 1.0/8.0*(5.0+2.0*r-sqrt(-7.0-12.0*r-4.0*r**2))
      else if ((r.ge.-1.0).and.(r.le.0.0)) then
        phi_r4 = 1.0/8.0*(3.0+2.0*r+sqrt(1.0-4.0*r-4.0*r**2))
      else if ((r.ge.0.0).and.(r.le.1.0)) then
        phi_r4 = 1.0/8.0*(3.0-2.0*r+sqrt(1.0+4.0*r-4.0*r**2))
      else if ((r.ge.1.0).and.(r.le.2.0)) then
        phi_r4 = 1.0/8.0*(5.0-2.0*r-sqrt(-7.0+12.0*r-4.0*r**2))
      else if (r.ge.2.0) then
        phi_r4 = 0.0
      end if

      return
      end function
!######################################################################
      real function phi_r4smth(r)
!######################################################################
      implicit none
      real, intent(in) :: r
      real :: PI,ar
       PI = 4.D0*DATAN(1.D0)
	ar=abs(r)
      if ((ar.ge.0.0).and.(ar.le.0.5)) then
	phi_r4smth = 3.0/8.0+ PI/32.0 - r*r/4.0
      else if ((ar.ge.0.5).and.(ar.le.1.5)) then
	phi_r4smth = 1.0/4.0 + (1.0-ar)/8.0*SQRT(-2.0+8.0*ar -4.0*r*r)
     & -1.0/8.0*ASIN(sqrt(2.0)*(ar-1.0))
      else if ((ar.ge.1.5).and.(ar.le.2.5)) then
	phi_r4smth = 17.0/16.0 - PI/64.0 - 3.0*ar/4.0+ r*r/8.0
     & + (ar-2.0)/16.0 * SQRT(-14.0+ 16.0*ar - 4.0*r*r)
     & + 1.0/16.0*ASIN(sqrt(2.0)*(ar-2.0))
      else if (ar.ge.2.5) then
        phi_r4smth = 0.0
      end if
      return
      end function


!######################################################################
      real function cd2_0(r)
!######################################################################
      implicit none
      real, intent(in) :: r
      if (r.le.-1.5) then
        cd2_0 = 0.0
      else if ((r.ge.-1.5).and.(r.le.-0.5)) then
        cd2_0 = 0.5*r+0.75
      else if ((r.ge.-0.5).and.(r.le.0.5)) then
        cd2_0 = 0.5
      else if ((r.ge.0.5).and.(r.le.1.5)) then
        cd2_0 = -0.5*r+0.75
      else if (r.ge.1.5) then
        cd2_0 = 0.0
      end if

      return
      end function

!######################################################################
      real function dcd2_0(r)
!######################################################################
      implicit none
      real, intent(in) :: r
      if (r.le.-1.5) then
        dcd2_0 = 0.0
      else if ((r.ge.-1.5).and.(r.le.-0.5)) then
        dcd2_0 = 0.5
      else if ((r.ge.-0.5).and.(r.le.0.5)) then
        dcd2_0 = 0.0
      else if ((r.ge.0.5).and.(r.le.1.5)) then
        dcd2_0 = -0.5
      else if (r.ge.1.5) then
        dcd2_0 = 0.0
      end if

      return
      end function

!######################################################################
      real function cd2_1(r)
!######################################################################
      implicit none
      real, intent(in) :: r
      if (r.le.-2.0) then
        cd2_1 = 0.0
      else if ((r.ge.-2.0).and.(r.le.-1.0)) then
        cd2_1 = 0.25*r**2+r+1.0
      else if ((r.ge.-1.0).and.(r.le.1.0)) then
        cd2_1 = -0.25*r**2+0.5
      else if ((r.ge.1.0).and.(r.le.2.0)) then
        cd2_1 = 0.25*r**2-r+1.0
      else if (r.ge.2.0) then
        cd2_1 = 0.0
      end if

      return
      end function


!######################################################################
      real function dcd2_1(r)
!######################################################################
      implicit none
      real, intent(in) :: r
      if (r.le.-2.0) then
        dcd2_1 = 0.0
      else if ((r.ge.-2.0).and.(r.le.-1.0)) then
        dcd2_1 = 0.5*r+1.0
      else if ((r.ge.-1.0).and.(r.le.1.0)) then
        dcd2_1 = -0.5*r
      else if ((r.ge.1.0).and.(r.le.2.0)) then
        dcd2_1 = 0.5*r-1.0
      else if (r.ge.2.0) then
        dcd2_1 = 0.0
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
!   for smooth ﬁelds (Uhlmann(2005) cf. Section 5.1.1).
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
      real*8 ::phi_r2smth,phi_r3smth,phi_r4,phi_r3,phi_r1smth,phi_r4smth
      real*8 :: cd2_0,cd2_1
      select case (order)

         case (1)
      dh =  phi_r1smth((xij-Xl)/dx)
     & * phi_r1smth((yij-Yl)/dy) * phi_r1smth((zij-Zl)/dz)
         case (2)
      dh =   phi_r2smth((xij-Xl)/dx)
     & * phi_r2smth((yij-Yl)/dy) * phi_r2smth((zij-Zl)/dz)
         case (3)
       dh =   phi_r3smth((xij-Xl)/dx)
     &   * phi_r3smth((yij-Yl)/dy) * phi_r3smth((zij-Zl)/dz)
         case (4)
       dh =   phi_r4smth((xij-Xl)/dx)
     &   * phi_r4smth((yij-Yl)/dy) * phi_r4smth((zij-Zl)/dz)
         case (5)
       dh =   phi_r3((xij-Xl)/dx)
     &   * phi_r3((yij-Yl)/dy) * phi_r3((zij-Zl)/dz)
         case (6)
      dh =   phi_r4((xij-Xl)/dx)
     &   * phi_r4((yij-Yl)/dy) * phi_r4((zij-Zl)/dz)
         case (7)
      dh =   cd2_0((xij-Xl)/dx)
     &   * cd2_0((yij-Yl)/dy) * cd2_0((zij-Zl)/dz)
         case (8)
      dh =   cd2_1((xij-Xl)/dx)
     &   * cd2_1((yij-Yl)/dy) * cd2_1((zij-Zl)/dz)

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
      implicit none
      real, intent(in) :: r
      if (r.le.-1.5) then
        dphi_r3 = 0.0
      else if ((r.ge.-1.5).and.(r.le.-0.5)) then
        dphi_r3 = 1.0/6.0*(3.0-0.5*(-3.0*(1.0+r)**2+1.0)**(-0.5)
     &*(-6.0*(1.0+r)))
      else if ((r.ge.-0.5).and.(r.le.0.5)) then
        dphi_r3 = 1.0/6.0*(-3.0*r**2+1.0)**(-0.5)*(-6.0*r)
      else if ((r.ge.0.5).and.(r.le.1.5)) then
        dphi_r3 = 1.0/6.0*(-3.0-0.5*(-3.0*(1.0-r)**2+1.0)**(-0.5)
     &*(6.0*(1.0-r)))
      else if (r.ge.1.5) then
        dphi_r3 = 0.0
      end if

      return
      end function




!######################################################################
      double precision
     &function ddh(dx,dy,dz,xij,yij,zij,Xl,Yl,Zl,order,dir)
!######################################################################

      implicit none

      real,    intent(in) :: dx,dy,dz,xij,yij,zij,Xl,Yl,Zl
      integer, intent(in) :: order,dir
      real*8 :: phi_r3,dphi_r3,cd2_0,dcd2_0,cd2_1,dcd2_1

      select case (order)

      case (5)
      select case (dir)

        case (1)
       ddh = dphi_r3((xij-Xl)/dx)
     &   * phi_r3((yij-Yl)/dy) * phi_r3((zij-Zl)/dz)
     &   * (-1.0/dx)

        case (2)
       ddh = phi_r3((xij-Xl)/dx)
     &   * dphi_r3((yij-Yl)/dy) * phi_r3((zij-Zl)/dz)
     &   * (-1.0/dy)

        case (3)
       ddh = phi_r3((xij-Xl)/dx)
     &   * phi_r3((yij-Yl)/dy) * dphi_r3((zij-Zl)/dz)
     &   * (-1.0/dz)

      end select




      case (7)
      select case (dir)

        case (1)
       ddh = dcd2_0((xij-Xl)/dx)
     &   * cd2_0((yij-Yl)/dy) * cd2_0((zij-Zl)/dz)
     &   * (-1.0/dx)

        case (2)
       ddh = cd2_0((xij-Xl)/dx)
     &   * dcd2_0((yij-Yl)/dy) * cd2_0((zij-Zl)/dz)
     &   * (-1.0/dy)

        case (3)
       ddh = cd2_0((xij-Xl)/dx)
     &   * cd2_0((yij-Yl)/dy) * dcd2_0((zij-Zl)/dz)
     &   * (-1.0/dz)

      end select


      end select

      return
      end function
