
      subroutine nnloevalfinite(yy,zz,logmh2omu2,nnloeps0) 
      implicit none 
      include 'types.f'
      include 'constants.f'
      include 'nf.f' 
      include 'ewcouple.f' 
      include 'masses.f'
      include 'scet_const.f' 
      integer i
      real(dp) xx
      real(dp) a1born,a2born

      real(dp) a1finite,a2finite
      real(dp) finitepartsum
      real(dp) finitepartoneonenoeps,finiteparttwotwo,finitepartonetwo,finiteparttwoone

      real(dp),intent(in)::yy,zz,logmh2omu2
      real(dp),intent(out)::nnloeps0

c--- initialize to zero
      nnloeps0=0._dp
      finitepartsum=0._dp

      a1finite=0._dp
      a2finite=0._dp

      finitepartoneonenoeps=0._dp
      finiteparttwotwo=0._dp
      finitepartonetwo=0._dp
      finiteparttwoone=0._dp

!==================================================================

c---- value for x
      xx = one-yy-zz

c---- form factors at LO with overall factor sqrt(as/(2*pi))
      a1born=two*rt2*pi*(1/yy+1/zz)
      a2born=two*rt2*pi*(-2/yy)

!==================================================================

!==== numerical evaluation of the MI part, infrared part, and renormalization part of A1,A2 after eps-expansion up to eps^0 and cancellation of poles

      call acoeffnnloevalfinite(yy,zz,logmh2omu2,a1finite,a2finite)

!==== construction of all required pieces of the interference

      finitepartoneonenoeps=four*yy*zz*a1finite*a1born
      finiteparttwotwo=two*xx**2*yy/zz*a2finite*a2born
      finitepartonetwo=-two*xx*yy*a1finite*a2born
      finiteparttwoone=-two*xx*yy*a2finite*a1born

      finitepartsum=two*(finitepartoneonenoeps+finiteparttwotwo+finitepartonetwo+finiteparttwoone)

      write(*,*) 'finitepartsum = ',finitepartsum

!==== multiply by LO color: Nc*cf
!==== divide by (2*pi)**3
!==== overall factor is as**3*yb**2 [mu**(2*eps) exp(eugam*eps)/(4*pi)**eps]

!==== NNLO hard (e.g. rescaled by LO msq and (as/2/pi)**2)

      nnloeps0=Nc*cf/(128._dp*pi**2)*yy*zz/(xx**2+one)*finitepartsum

      write(6,99) ' eps  0: ',nnloeps0

      return 
      
   99 format(a35,e20.12)
      
      end

