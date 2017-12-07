!---  CW Sept 17
!---  this is a piece of the NLO x NLO  qqb => Z (H=>bb) process
!---  It corresponds to : The Real qqb => ZH  x the below cut H=>bb (i.e. Real - Virtual) 
      subroutine qqb_zh_bb_nlopxd_g(p,msqv)
      implicit none
      include 'types.f'
      
      integer:: j,k
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      real(dp):: p(mxpart,4),
     & msqv(-nf:nf,-nf:nf),msq(-nf:nf,-nf:nf),
     & dot,virt,xl12

      xl12=log(two*dot(p,1,2)/musq)
      scheme='dred'
      call spinoru(6,p,za,zb)


!------ REAL PRODUCTION 


!------ BELOW CUT FACTORIZATION H=>bb
      s56=s(5,6)
      call Hbb_bc_bit(s56,1,hdecay)
      
!----- COMBINATION 



      
      return
      end
