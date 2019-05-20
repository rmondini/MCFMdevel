
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
      real(dp) finitepartsum
      real(dp) hplevalfinite

      real(dp),intent(in)::yy,zz,logmh2omu2
      real(dp),intent(out)::nnloeps0

c--- initialize to zero
      nnloeps0=0._dp
      finitepartsum=0._dp

!==================================================================

c---- value for x
      xx = one-yy-zz

!==================================================================

!==== numerical evaluation of the MI part, infrared part, and renormalization part of A1,A2 after eps-expansion up to eps^0 and cancellation of poles

      finitepartsum = hplevalfinite(yy,zz,logmh2omu2)

!==== multiply by LO color: Nc*cf
!==== divide by (2*pi)**3
!==== overall factor is as**3*yb**2 [mu**(2*eps) exp(eugam*eps)/(4*pi)**eps]

!==== NNLO hard (e.g. rescaled by LO msq and (as/2/pi)**2)

      nnloeps0=Nc*cf/(128._dp*pi**2)*yy*zz/(xx**2+one)*finitepartsum

!      write(6,99) ' eps  0: ',nnloeps0

      return 
      
   99 format(a35,e20.12)
      
      end

