
      subroutine acoeffnlorenormeval(yy,zz,logmh2omu2,a1renorm,a2renorm) 
      implicit none 
      include 'types.f'
      include 'constants.f' 
      include 'ewcouple.f' 
      include 'masses.f' 
      include 'nf.f'
      include 'scet_const.f'
      real(dp),intent(in)::yy,zz,logmh2omu2
      real(dp),intent(out)::a1renorm(5),a2renorm(5)

      double precision Li2,Li3,Li4
      external Li2,Li3,Li4

c--- initialize to zero
      a1renorm(:)=0._dp
      a2renorm(:)=0._dp

!==================================================================

!==== eps^(-1)
      a1renorm(2) = ((-11*ca)/6. - 3*cf + nf/3.)*pi*rt2*(1/yy + 1/zz)

!==== eps^(-1)
      a2renorm(2) = ((11*ca + 18*cf - 2*nf)*pi*rt2)/(3.*yy) 

!===================================================================

      return 
      
   99 format(a35,e15.8)
      
      end

