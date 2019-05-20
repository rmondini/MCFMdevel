
      subroutine nnloevalfinite_proc1(uu,vv,musq,nnloeps0) 
      implicit none 
      include 'types.f'
      include 'constants.f'
      include 'nf.f' 
      include 'ewcouple.f' 
      include 'masses.f'
      include 'scet_const.f' 
      integer i,i1,i2,i3,i4
      real(dp) xx,yy,zz
      real(dp) a1born,a2born

      complex(dp) a1finite,a2finite
      complex(dp) finitepartoneonenoeps,finiteparttwotwo,finitepartonetwo,finiteparttwoone

      real(dp) finitepartsum

      real(dp),intent(in)::uu,vv,musq
      real(dp),intent(out)::nnloeps0 

c---- initialize to zero
      nnloeps0=0._dp
      finitepartsum=0._dp

!==================================================================

c---- value for x,y,z
      xx = one/vv
      yy = -uu/vv
      zz = one-xx-yy

c---- form factors at LO with overall factor sqrt(as/(2*pi))
      a1born=two*rt2*pi*(1/yy+1/zz)
      a2born=two*rt2*pi*(-2/yy)

!==================================================================

!==== MSbar-renormalized coefficients for hard function

      call nnloevalfinitetdhpla1_proc1(uu,vv,musq,a1finite)
      call nnloevalfinitetdhpla2_proc1(uu,vv,musq,a2finite)

!      write(*,*) 'a1finite nnlo proc1 = ',a1finite
!      write(*,*) 'a2finite nnlo proc1 = ',a2finite

!==== construction of all required pieces of the interference

      finitepartoneonenoeps=four*yy*zz*a1finite*a1born
      finiteparttwotwo=two*xx**2*yy/zz*a2finite*a2born
      finitepartonetwo=-two*xx*yy*a1finite*a2born
      finiteparttwoone=-two*xx*yy*a2finite*a1born

      finitepartsum=two*real(finitepartoneonenoeps+finiteparttwotwo+finitepartonetwo+finiteparttwoone,dp)

!==================================================================

!==== multiply by LO color: Nc*cf
!==== divide by (2*pi)**3
!==== overall factor is as**3*yb**2
!==== NNLO hard (e.g. rescaled by LO msq and (as/2/pi)**2)

!==== no overall -1:
!==== for proc1 because we are crossing two fermions from final to initial
!==== for proc2 we are crossing one fermion from final to initial so
!==== there would be a -1 in finitepartsum but then it would be canceled when rescaling by LO msq

      nnloeps0=1/(32._dp*pi**2)*yy*zz/(xx**2+one)*finitepartsum

!      write(*,*) 'nnloeps0 proc1 = ',nnloeps0

      return 
      
   99 format(a35,e20.12) 
      end


