      subroutine maketaucut_inc(pparton,pincj,scetjets,isub,passed)
c---- CW. July 16 Modified version to allow for inclusive style analysis 
c---- J. Campbell, April 2015
c----
c---- Adapted from code written by R. Boughezal et al. for
c---- 1-jettiness calculation in W+jet events
c----
c---- Implements cut on 1-jettiness for processes with nqcdjets=1
c---- Implements cut on 0-jettiness for processes with nqcdjets=0
c----
c---  Only keeps events with tau > taucut (passed via common block)
c---  tau < taucut will be treated separately using SCET
c----
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'npart.f'
      include 'nqcdjets.f'
      include 'taucut.f'
      include 'plabel.f'
      include 'first.f'
!      include 'scetincj.f'
      logical passed
      integer i,j,isub,ihard,scetjets
      real(dp) pparton(mxpart,4),pincj(mxpart,4),getdynamictau,
     &     ptjet(3),nn1(4),nn2(4),nn3(4),tau,taua,taub,tauj
      real(dp):: ppin(4),ppout(4),Q(4),Qrest(4),Qsq
      real(dp) tauhist
      common/tauhist/tauhist
!$omp threadprivate(/tauhist/)
      integer, save :: ipp
!$omp threadprivate(ipp)
      logical failedndp 
      common/failedndp/failedndp 
!$omp threadprivate(/failedndp/) 

c--- determine beginning of parton entries in plabel
      if (first) then
        first=.false.
        ipp=3
        do while ((ipp < mxpart ) .and. (plabel(ipp) .ne. 'pp'))
          ipp=ipp+1
        enddo
        if (ipp == mxpart) then
          write(6,*) 'Could not identify partons in maketaucut.f'
          stop
        endif
c        write(6,*) 'found ipp=',ipp  
      endif
       
      passed=.false.
      tauhist=zip
!----- in this routine we will use pincj as pincj (obtained from clustering) 
!      pincj(:,:)=pincj(:,:)
!      jets=scetjets
      
!      if(scetjets.ne.1) then 
!      call writeout(pincj) 
!      write(6,*) scetjets
!      pause
!      endif
c--- if there are not enough partons (i.e. we're doing virtual) then return
      if ((npart .eq. ipp-3+ntau) .and. (isub .eq. 0)) then
        passed=.true.
        return
      endif

!--- routine is not designed for more than 3 jets
      if (scetjets .gt. 3) then
        write(6,*) 'Error: >3 jets found in maketaucut.f'
        write(6,*) ' scetjets=',scetjets
        stop
      endif

!--- routine is not designed for nqcdjets>2
      if (nqcdjets .gt. 2) then
        write(6,*) 'Error: unimplemented nqcdjets in maketaucut.f'
        write(6,*) ' nqcdjets=',nqcdjets
        stop
      endif

c--- reference directions for beam axes
      nn1(:)=(/0._dp,0._dp, 1._dp,1._dp/)
      nn2(:)=(/0._dp,0._dp,-1._dp,1._dp/)

c---  compute hard jet reference direction for 1-jettiness
      if (ntau > 0) then
        ptjet(:)=0._dp
        do i=0,scetjets-1
          ptjet(i+1)=sqrt(pincj(ipp+i,1)**2+pincj(ipp+i,2)**2)
        enddo
        ihard=ipp
        if (scetjets > 1) then
          if (ptjet(2) .gt. max(ptjet(1),ptjet(3))) ihard=ipp+1
        endif
        if (scetjets > 2) then
          if (ptjet(3) .gt. max(ptjet(1),ptjet(2))) ihard=ipp+2
        endif
        nn3(1:3)=pincj(ihard,1:3)
     &    /sqrt(pincj(ihard,1)**2+pincj(ihard,2)**2+pincj(ihard,3)**2)
        nn3(4)=1.0_dp
      endif

      tau=0._dp
c--- compute N-jettiness tau, by dotting the final state QCD parton
c--- 4-momentum into ni, for i=1,2 and 3 (ntau>0 only)
      do j = ipp,npart+2-isub     
      taua =pparton(j,4)*nn1(4)-pparton(j,1)*nn1(1)
     &     -pparton(j,2)*nn1(2)-pparton(j,3)*nn1(3)
      taub=pparton(j,4)*nn2(4)-pparton(j,1)*nn2(1)
     &    -pparton(j,2)*nn2(2)-pparton(j,3)*nn2(3)
      if (ntau .eq. 0) then
        tau=tau+min(abs(taua),abs(taub))
      else
        tauj=pparton(j,4)*nn3(4)-pparton(j,1)*nn3(1)
     &      -pparton(j,2)*nn3(2)-pparton(j,3)*nn3(3)
        tau=tau+min(abs(taua),abs(taub),abs(tauj))
      endif
      enddo

      if ((1 == 2) .and. (tauboost)) then ! begin clumsy commenting-out

c--- alternative method: compute all quantities in Born frame where
c--- entire Born system is at rest
      Q(:)=pparton(3,:)
      if (ipp > 4) Q(:)=Q(:)+pparton(4,:)
      if (ipp > 5) Q(:)=Q(:)+pparton(5,:)
      if (ipp > 6) Q(:)=Q(:)+pparton(6,:)
      if (ntau > 0) Q(:)=Q(:)+pincj(ihard,:)
      Qsq=Q(4)**2-Q(1)**2-Q(2)**2-Q(3)**2
      Qrest(:)=zip
      Qrest(4)=sqrt(Qsq)
      ppin(:)=-pparton(1,:)
      call boostx(ppin,Q,Qrest,ppout)
      nn1(1:3)=ppout(1:3)/sqrt(ppout(1)**2+ppout(2)**2+ppout(3)**2)
      nn1(4)=1.0_dp
      ppin(:)=-pparton(2,:)
      call boostx(ppin,Q,Qrest,ppout)
      nn2(1:3)=ppout(1:3)/sqrt(ppout(1)**2+ppout(2)**2+ppout(3)**2)
      nn2(4)=1.0_dp
      if (ntau > 0) then
        ppin(:)=pincj(ihard,:)
        call boostx(ppin,Q,Qrest,ppout)
        nn3(1:3)=ppout(1:3)/sqrt(ppout(1)**2+ppout(2)**2+ppout(3)**2)
        nn3(4)=1.0_dp
      endif
      
      tau=0._dp
c--- compute N-jettiness tau, by dotting the final state QCD parton
c--- 4-momentum into ni, for i=1,2 and 3 (ntau>0 only)
      do j = ipp,npart+2-isub     
        ppin(:)=pparton(j,:)
        call boostx(ppin,Q,Qrest,ppout)
        taua=ppout(4)*nn1(4)-ppout(1)*nn1(1)
     &      -ppout(2)*nn1(2)-ppout(3)*nn1(3)
        taub=ppout(4)*nn2(4)-ppout(1)*nn2(1)
     &      -ppout(2)*nn2(2)-ppout(3)*nn2(3)
        if (ntau .eq. 0) then
          tau=tau+min(abs(taua),abs(taub))
        else
          tauj=ppout(4)*nn3(4)-ppout(1)*nn3(1)
     &        -ppout(2)*nn3(2)-ppout(3)*nn3(3)
          tau=tau+min(abs(taua),abs(taub),abs(tauj))
      endif
!      write(6,*) 'new taua,taub',taua,taub
      enddo
!      write(6,*)
      
      endif ! end commenting-out
      tauhist=tau

c--- check to make sure no NaN
      if (tau .ne. tau) then
        call writeout(pparton)
        write(6,*) 'maketaucut.f:  tau=',tau
        stop
      endif
c      write(6,*) 'tau',tau
      
c--- if tau is too small, implement the jettiness cut to remove the event
      if (tau .lt. getdynamictau(pparton)) return

c--- at this point cut has been passed
      passed=.true.
c      if (failedndp) then 
c         write(6,*) 'WARNING, phase space point has passed, '
c         write(6,*) 'both tau cut AND has failed ndp check'
c      endif

      return
      end
      
