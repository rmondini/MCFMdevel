      implicit none
      include 'types.f'
      real(dp):: y12,y23,y31,soft1(-1:1),soft2(-1:3),al,be,
     & I0JSTW,I1JSTW
      integer::j

      y12=1._dp
      open(unit=60,file='ggg_abelian.dat')
      open(unit=61,file='ggg_nonabelian.dat')
      open(unit=62,file='qgq_abelian.dat')
      open(unit=63,file='qgq_nonabelian.dat')
      open(unit=64,file='qag_abelian.dat')
      open(unit=65,file='qag_nonabelian.dat')
      do j=1,199
      y31=dfloat(j)/200._dp
      y23=1._dp-y31
c      write(6,*) 2*y31,soft2(3)/4._dp,soft2(2)/3._dp,
c     & soft2(1)/2._dp,soft2(0),soft2(-1)
      call soft_ab_ggg(2,y12,y23,y31,soft1,soft2)
      write(60,*) 2*y31,soft2(3)/4._dp,soft2(2)/3._dp,
     & soft2(1)/2._dp,soft2(0),soft2(-1)
      call soft_nab_ggg(2,y12,y23,y31,soft1,soft2)
      write(61,*) 2*y31,soft2(3)/4._dp,soft2(2)/3._dp,
     & soft2(1)/2._dp,soft2(0),soft2(-1)
      call soft_ab_qgq(2,y12,y23,y31,soft1,soft2)
      write(62,*) 2*y31,soft2(3)/4._dp,soft2(2)/3._dp,
     & soft2(1)/2._dp,soft2(0),soft2(-1)
      call soft_nab_qgq(2,y12,y23,y31,soft1,soft2)
      write(63,*) 2*y31,soft2(3)/4._dp,soft2(2)/3._dp,
     & soft2(1)/2._dp,soft2(0),soft2(-1)
      call soft_ab_qag(2,y12,y23,y31,soft1,soft2)
      write(64,*) 2*y31,soft2(3)/4._dp,soft2(2)/3._dp,
     & soft2(1)/2._dp,soft2(0),soft2(-1)
      call soft_nab_qag(2,y12,y23,y31,soft1,soft2)
      write(65,*) 2*y31,soft2(3)/4._dp,soft2(2)/3._dp,
     & soft2(1)/2._dp,soft2(0),soft2(-1)

      enddo
      stop
      end
