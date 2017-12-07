      function I1JSTW(al,be)
! function I1JSTW from 1102.4344
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp)::I1JSTW
      real(dp)::al,be,rtal,rtbe,phimax,phicut,dgauss,
     & I1integrand1,I1integrand2,eps
      common/albe/rtal,rtbe
      external I1integrand1,I1integrand2
!$omp threadprivate(/albe/)
      eps=1.e-3_dp
      rtal=sqrt(al)
      rtbe=sqrt(be)
!--setup limits
      if (abs(rtal-rtbe) .ge. one) then
      phicut=zip
      elseif (rtal+rtbe .le. one) then
      phicut=pi
      else
      phicut=acos((al+be-one)/(two*rtal*rtbe))
      endif
!--perform integration
      if (phicut > zip) then
        I1JSTW=dgauss(I1integrand1,zip,phicut,eps)
      else
        I1JSTW=zip
      endif

      if (al-be-one .gt. zip) then
      phimax=asin(one/rtal)
      I1JSTW=I1JSTW
     & +dgauss(I1integrand2,phicut,phimax,eps)
      endif
      I1JSTW=two/pi*I1JSTW
      return
      end


      function I1integrand1(phi)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp)::I1integrand1
      real(dp)::rtal,rtbe,phi,y,yp,G,rt
      complex(dp)::XSPENZ
      common/albe/rtal,rtbe
!$omp threadprivate(/albe/)
!---statement function
      G(y,phi)=-two*real(XSPENZ(y*exp(im*phi)))
!---end statement function
      rt=one/rtal**2-sin(phi)**2
      if (rt > 0) then
        rt=sqrt(rt)
      else
        rt=zip
      endif
      yp=cos(phi)+rt
      I1integrand1=G(yp,phi)-G(rtbe/rtal,phi)
      return
      end

      function I1integrand2(phi)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp)::I1integrand2
      real(dp)::rtal,rtbe,phi,y,yp,ym,G,rt
      complex(dp)::XSPENZ
      common/albe/rtal,rtbe
!$omp threadprivate(/albe/)
!---statement function
      G(y,phi)=-two*real(XSPENZ(y*exp(im*phi)))
!---end statement function
      rt=one/rtal**2-sin(phi)**2
      if (rt > 0) then
        rt=sqrt(rt)
      else
        rt=zip
      endif
      yp=cos(phi)+rt
      ym=cos(phi)-rt
      I1integrand2=G(yp,phi)-G(ym,phi)
      return
      end

