      subroutine applyRe(epsseries)
c--- Given a series of coefficients in 1/epsilon
c--- multiply by cos(pi*eps)
c--- Note: coefficients are returned in place
      implicit none
      include 'types.f'
      include 'constants.f'
      real(kind=dp) epsseries(5)

c--- e0
      epsseries(5)=epsseries(5)-pi**2/2*epsseries(3)+pi**4/24*epsseries(1)
c--- e1
      epsseries(4)=epsseries(4)-pi**2/2*epsseries(2)
c--- e2
      epsseries(3)=epsseries(3)-pi**2/2*epsseries(1)
c--- e3
      epsseries(2)=epsseries(2)
c--- e4
      epsseries(1)=epsseries(1)

      return
      end
      
