!--   bit which assembles just below cut pieces for H=>bb at a given order

      subroutine Hbb_bc_bit(s12,E1,E2,order,xmsq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'scale.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'hbbparams.f'
      real(dp):: s12,xmsq,E1,E2,L12
      real(dp)::  xx(2),soft1(-1:1),soft2(-1:3),hard(2)
      real(dp) :: jeta1(-1:1),jeta2(-1:3),jetb1(-1:1),jetb2(-1:3)
      integer order
      logical higgsrestframe
      real(dp) :: assemble_dec,fac

      fac=xn*gwsq*mb_eff**2/(4._dp*wmass**2)*2._dp*s12
      
      higgsrestframe=.true. 
           
     
      if(higgsrestframe) then
         L12=zip
         call softqqbis(order,soft1,soft2)
         call jetq(order,sqrt(s12),jeta1,jeta2)
         call jetq(order,sqrt(s12),jetb1,jetb2)
      else
!---- constuct log for soft function 
         L12=log(s12/E1/E2/4._dp)
c---  two jets functions for b-bar 
         call jetq(order,two*E1,jeta1,jeta2)
         call jetq(order,two*E2,jetb1,jetb2)
         call soft_qqdec(order,soft1,soft2,L12)
      endif
      
      
!---- hard function 
      call hard_hbb(s12,musq,hard)

    
      xmsq=assemble_dec(order,
     & 1._dp,1._dp,jeta1,jetb1,
     & jeta2,jetb2,soft1,soft2,hard)


      xmsq=xmsq*fac

      
      return
      end
