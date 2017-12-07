      function ggHexpand(expandWilson)
      implicit none
      include 'types.f'
      include 'nflav.f'
      include 'masses.f'
      include 'scale.f'
      include 'qcdcouple.f'
      integer:: expandWilson
      real(dp):: ggHexpand
      
!      write(6,*) 'C1',ason2pi*11._dp
!      write(6,*) 'C2',ason2pi**2*(
!     &    1933._dp/18._dp - 67._dp/12._dp*real(nflav,dp)
!     &    -(19._dp/2._dp + 8._dp/3._dp*real(nflav,dp))*log(mt**2/musq))

      if (expandWilson < 1) then
        ggHexpand=1._dp
        return
      endif
      
      if (expandWilson == 1) then
        ggHexpand=1._dp+ason2pi*11._dp
        return
      endif
      
      if (expandWilson == 2) then
        ggHexpand=1._dp+ason2pi*11._dp+ason2pi**2*(
     &    1933._dp/18._dp - 67._dp/12._dp*real(nflav,dp)
     &    -(19._dp/2._dp + 8._dp/3._dp*real(nflav,dp))*log(mt**2/musq))
        return
      endif
      
      write(6,*) 'Unexpected parameter: expandWilson = ',expandWilson
      stop
      
      return
      end
