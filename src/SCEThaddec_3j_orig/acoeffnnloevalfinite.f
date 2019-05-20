
      subroutine acoeffnnloevalfinite(yy,zz,logmh2omu2,a1finite,a2finite) 
      implicit none 
      include 'types.f'
      include 'constants.f'  
      include 'ewcouple.f' 
      include 'masses.f' 
      include 'nf.f'
      include 'scet_const.f'
      integer i
      real(dp) hplevalfinite
      real(dp),intent(in)::yy,zz,logmh2omu2
      real(dp),intent(out)::a1finite,a2finite

c--- initialize to zero
      a1finite=0._dp
      a2finite=0._dp

!==================================================================

!==== 5  = A1,eps^0
!==== 10 = A2,eps^0

      a1finite = hplevalfinite(yy,zz,logmh2omu2,5)
      a2finite = hplevalfinite(yy,zz,logmh2omu2,10)

!      write(*,*)
!      write(*,*) 'a1finite = ',a1finite
!      write(*,*)
!      write(*,*) 'a2finite = ',a2finite

!===================================================================

      return 
      
   99 format(a35,e15.8)
      
      end

