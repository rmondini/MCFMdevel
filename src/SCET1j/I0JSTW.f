      function I0JSTW(al,be)
! function I0 from 1102.4344
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp)::I0JSTW
      real(dp)::al,be,rtal,rtbe,phimax,phicut,dgauss,
     & I0integrand1,I0integrand2,eps
      common/albe/rtal,rtbe
      external I0integrand1,I0integrand2
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
        I0JSTW=dgauss(I0integrand1,zip,phicut,eps)
      else
        I0JSTW=zip
      endif
      if (al-be-one .gt. zip) then
      phimax=asin(one/rtal)
      I0JSTW=I0JSTW
     & +dgauss(I0integrand2,phicut,phimax,eps)
      endif
      I0JSTW=two/pi*I0JSTW
      return
      end

      function I0integrand1(phi)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp)::I0integrand1
      real(dp)::rtal,rtbe,phi,yp,rt
      common/albe/rtal,rtbe
!$omp threadprivate(/albe/)
      rt=one/rtal**2-sin(phi)**2
      if (rt > 0) then
        rt=sqrt(rt)
      else
        rt=zip
      endif
      yp=cos(phi)+rt
      I0integrand1=log(yp/(rtbe/rtal))
      return
      end

      function I0integrand2(phi)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp)::I0integrand2
      real(dp)::rtal,rtbe,phi,yp,ym,rt
      common/albe/rtal,rtbe
!$omp threadprivate(/albe/)
      rt=one/rtal**2-sin(phi)**2
      if (rt > 0) then
        rt=sqrt(rt)
      else
        rt=zip
      endif
      yp=cos(phi)+rt
      ym=cos(phi)-rt
      I0integrand2=log(yp/ym)
      return
      end



