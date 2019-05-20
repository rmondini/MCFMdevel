
      subroutine nloxnlofeval(yy,zz,musq,nloxnloeps0) 
      implicit none 
      include 'types.f'
      include 'constants.f'
      include 'nf.f'  
      include 'ewcouple.f' 
      include 'masses.f'
      include 'scet_const.f'
      integer i
      real(dp) xx

      complex(dp) a1finite,a2finite
      complex(dp) lnrat

      real(dp),intent(in)::yy,zz,musq
      real(dp),intent(out)::nloxnloeps0

      double precision Li2
      external Li2

c---- initialize to zero
      nloxnloeps0=0._dp

c---- value for x
      xx = one-yy-zz

!      write(*,*) 'xx=',xx
!      write(*,*) 'musq=',musq
!      write(*,*) 'hmass**2=',hmass**2
!      write(*,*) 'logminusmh2omu2=',lnrat(-hmass**2,musq)

      a1finite=(pi*(yy + zz)*(4 - 4*cf*Nc + pi**2 - 8*Li2(xx) + 
     -      8*cf*Nc*Li2(yy) + 8*cf*Nc*Li2(zz) - 
     -      2*lnrat(-hmass**2,musq)**2 - 
     -      8*cf*Nc*lnrat(-hmass**2,musq)**2 + 
     -      4*lnrat(-hmass**2,musq)*log(xx) - 8*log(1 - xx)*log(xx) + 
     -      2*log(xx)**2 - 4*lnrat(-hmass**2,musq)*log(yy) - 
     -      8*cf*Nc*lnrat(-hmass**2,musq)*log(yy) + 4*log(xx)*log(yy) + 
     -      8*cf*Nc*log(1 - yy)*log(yy) - 2*log(yy)**2 - 
     -      4*cf*Nc*log(yy)**2 - 4*lnrat(-hmass**2,musq)*log(zz) - 
     -      8*cf*Nc*lnrat(-hmass**2,musq)*log(zz) + 4*log(xx)*log(zz) - 
     -      4*log(yy)*log(zz) - 8*cf*Nc*log(yy)*log(zz) + 
     -      8*cf*Nc*log(1 - zz)*log(zz) - 2*log(zz)**2 - 
     -      4*cf*Nc*log(zz)**2))/(2.*Nc*rt2*yy*zz) 

      a2finite=(pi*(8*cf*Nc - pi**2 + 8*Li2(xx) - 8*cf*Nc*Li2(yy) - 
     -      8*cf*Nc*Li2(zz) + 2*lnrat(-hmass**2,musq)**2 + 
     -      8*cf*Nc*lnrat(-hmass**2,musq)**2 - 
     -      4*lnrat(-hmass**2,musq)*log(xx) + 8*log(1 - xx)*log(xx) - 
     -      2*log(xx)**2 + 4*lnrat(-hmass**2,musq)*log(yy) + 
     -      8*cf*Nc*lnrat(-hmass**2,musq)*log(yy) - 4*log(xx)*log(yy) - 
     -      8*cf*Nc*log(1 - yy)*log(yy) + 2*log(yy)**2 + 
     -      4*cf*Nc*log(yy)**2 + 4*lnrat(-hmass**2,musq)*log(zz) + 
     -      8*cf*Nc*lnrat(-hmass**2,musq)*log(zz) - 4*log(xx)*log(zz) + 
     -      4*log(yy)*log(zz) + 8*cf*Nc*log(yy)*log(zz) - 
     -      8*cf*Nc*log(1 - zz)*log(zz) + 2*log(zz)**2 + 
     -      4*cf*Nc*log(zz)**2))/(Nc*rt2*yy)

!      write(*,*) 'a1finite=',a1finite
!      write(*,*) 'a2finite=',a2finite

      nloxnloeps0=
     & + four*yy*zz*(conjg(a1finite)*a1finite)
     & + two*xx**2*yy/zz*(conjg(a2finite)*a2finite)
     & - two*xx*yy*(conjg(a1finite)*a2finite)
     & - two*xx*yy*(conjg(a2finite)*a1finite)

!      write(*,*) 'nloxnloeps0=',nloxnloeps0

!==== multiply by LO color: Nc*cf
!==== divide by (2*pi)**3
!==== overall factor is as**3*yb**2 [mu**(2*eps) exp(eugam*eps)/(4*pi)**eps]
!==== NLO^2 hard (e.g. rescaled by LO msq and (as/2/pi)**2)

      nloxnloeps0=1/(32._dp*pi**2)*yy*zz/(xx**2+one)*nloxnloeps0

      return 
      
   99 format(a35,e20.12)
      
      end

