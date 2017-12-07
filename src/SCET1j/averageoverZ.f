      subroutine averageoverZ(fxn,p,msq)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'masses.f'
      integer:: jcount
      real(dp):: msq(-nf:nf,-nf:nf),msq1(-nf:nf,-nf:nf),
     & p(mxpart,4),pout(mxpart,4)
      pout(:,:)=0._dp
      msq(:,:)=0._dp
C---sum over 6 values;
      do jcount=1,6
      call pgen(p,jcount,3,4,pout)
      call fxn(pout,msq1,2)
      msq(:,:)=msq(:,:)+msq1(:,:)
      enddo
      msq(:,:)=msq(:,:)*(zwidth**2/(4._dp*esq*(le**2+re**2)))
 
      return
      end
