
      subroutine qqb_wgamj_v(p,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f' 
      include 'nf.f'
    
      real(dp)::  msq(-nf:nf,-nf:nf),p(mxpart,4)
c$$$       
      complex(dp):: test
      complex(dp):: Wgajprp34_Qu_mhv_lc


      msq(:,:)=zip
      p(:,:)=zip
!------KC testing code
      epinv=0._dp
      epinv2=0._dp
      musq=1._dp
      scale=1._dp
      write(6,*) 'heyup'
      include 'kcpoint_Wgaj.f'
      call spinorz(6,p,za,zb)
      call writeout(p)
      write(6,*) za(1,2),za(3,4)
      test=Wgajprp34_Qu_mhv_lc(1,2,3,4,5,6,za,zb)
      stop
      
      return
      end
      
      
