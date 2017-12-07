      subroutine gg_Hgaga_Hpi_gs(p,msq)
      implicit none
      include 'types.f'
c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c----for interference between H->ga ga and continuum background 
c     g(-p1)+g(-p2) --> [H -->] ga(p3)+ga(p4) + g(p5)
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'qqgg.f'
      integer:: j,k
      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      real(dp):: msq15_2(-nf:nf,-nf:nf),msq25_1(-nf:nf,-nf:nf),
     & msq15_2v(-nf:nf,-nf:nf),msq25_1v(-nf:nf,-nf:nf),
     & sub15_2(4),sub25_1(4),sub15_2v,sub25_1v
      external gg_Hgaga_Hpi,gg_Hgaga_Hpi_gvec
      ndmax=2

c---- calculate both initial-initial dipoles
c---- note that we do not require the gg dipoles, so the v-type
c---- entries are left as dummies
      call dips(1,p,1,5,2,sub15_2,sub15_2v,msq15_2,msq15_2v,
     & gg_Hgaga_Hpi,gg_Hgaga_Hpi_gvec)
      call dips(2,p,2,5,1,sub25_1,sub25_1v,msq25_1,msq25_1v,
     & gg_Hgaga_Hpi,gg_Hgaga_Hpi_gvec)

      do j=-nf,nf
      do k=-nf,nf
      msq(1:2,j,k)=0._dp

! Quark-channel subtractions removed
      if     ((j .ne. 0) .and. (k == 0)) then
         msq(1,j,k)=2._dp*cf
     &   *(msq15_2(0,0)*sub15_2(gq)+msq15_2v(0,0)*sub15_2v)
      elseif ((j == 0) .and. (k .ne. 0)) then
         msq(2,j,k)=2._dp*cf
     &   *(msq25_1(0,0)*sub25_1(gq)+msq25_1v(0,0)*sub25_1v)
      elseif ((j == 0) .and. (k == 0)) then
         msq(1,j,k)=2._dp*xn
     &   *(msq15_2(j,k)*sub15_2(gg)+msq15_2v(j,k)*sub15_2v)
         msq(2,j,k)=2._dp*xn
     &   *(msq25_1(j,k)*sub25_1(gg)+msq25_1v(j,k)*sub25_1v)
!         write(6,*) 'dipole 1',2._dp*xn
!     &   *(msq15_2(j,k)*sub15_2(gg)),2._dp*xn*msq15_2v(j,k)*sub15_2v
!         write(6,*) 'dipole 2',2._dp*xn
!     &   *(msq25_1(j,k)*sub25_1(gg)),2._dp*xn*msq25_1v(j,k)*sub25_1v
      endif

      enddo
      enddo

      return
      end


