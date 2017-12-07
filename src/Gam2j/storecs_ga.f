      subroutine storecs_ga(mcs)
      implicit none
      include 'types.f'
c-- this routine transfers the information on the colour structure
c-- for the Z2jet matrix elements into separate arrays for each
c-- incoming parton case
      
      include 'mmsq_cs_ga.f'
      integer:: i,j,k
      real(dp):: mcs(0:2)
      
      do i=0,2
!        do j=1,2
!          do k=1,2
         mcs(i)=mmsq_cs_ga(i)
!     mcs(i,j,k)=mmsq_cs(i,j,k)
!          enddo
!        enddo
      enddo
      return
      end
