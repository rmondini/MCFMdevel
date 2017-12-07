      subroutine maketaucut_dec(pparton,pincj,scetjets,isub,passed)
!====== C Williams March 2017
!======Special version of the routine maketaucuts designed for
!=====Hadronic decays of a vector boson, this routine ignores the beams
!===== in the construction of tau 
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
      logical passed,ttest(2)
      integer jj
      real(dp) :: ycut
      integer i,j,isub,ihard,scetjets
      real(dp) pparton(mxpart,4),pincj(mxpart,4),
     &     ptjet(3),nn1(4),nn2(4),nn3(4),tau,taua,taub,tauj
      real(dp) tauhist,s56
      common/tauhist/tauhist
!$omp threadprivate(/tauhist/)
      integer, save :: ipp
!     $omp threadprivate(ipp)
      integer jets
      logical failedndp 
      common/failedndp/failedndp
      real(dp) :: ptmax,ptmax2
      integer ihard2
      logical is_hadronic
!$omp threadprivate(/failedndp/) 
      if(ntau.ne.2) then
         write(6,*) 'maketaucut_dec.f is designed only for'
         write(6,*) 'use when ntau=2 found ntau=',ntau
         stop
      endif
      ihard=0
      ihard2=0
   
      
c--- determine beginning of parton entries in plabel
      if (first) then
         first=.false.
         ipp=3
         do while ((ipp < mxpart ) .and.(.not.(is_hadronic(ipp))))
            ipp=ipp+1
         enddo
         if (ipp == mxpart) then
            write(6,*) 'Could not identify partons in maketaucut.f'
            stop
         endif
c        write(6,*) 'found ipp=',ipp  
      endif

      passed=.false.
!      return
!      if(scetjets < 2) then
!         call writeout(pparton)
!         write(6,*) '*********** '
!         call writeout(pincj)
!         write(6,*) scetjets
!         pause
!      endif
      
!---  Demand at least (ntau) hard jets after implementing the jet
!---  algorithm, otherwise the event is removed
!      if (scetjets .lt. ntau) return

c---  if there are not enough partons (i.e. we're doing virtual) then return
      if ((npart .eq. ipp-3+ntau) .and. (isub .eq. 0)) then
!        passed=.true.
!        return
      endif

      goto 12

      
      
!----- find the hardest jet first
      ptjet(:)=0._dp
      do i=1,scetjets
         ptjet(i)=sqrt(pincj(i+ipp-1,1)**2+pincj(i+ipp-1,2)**2)
      enddo
      ihard=1
      ihard2=0
      ptmax=ptjet(1)
      ptmax2=0._dp
      do i=2,scetjets
         if(ptjet(i).gt.ptmax) then
            ihard=i
            ptmax=ptjet(i)
         endif
      enddo

      do i=1,scetjets
         if((ptjet(i).gt.ptmax2).and.(i.ne.ihard)) then
            ptmax2=ptjet(i)
            ihard2=i
         endif
      enddo
 
      ihard=ipp+ihard-1
      ihard2=ipp+ihard2-1

!      goto 13
      
!-----Durham option      
 12   continue
!      do jj=1,2
!         if(jj==1) ycut=0.1_dp
!     if(jj==2) ycut=0.01_dp
      ycut=0.01_dp
      s56=(-pparton(1,4)-pparton(2,4)-pparton(3,4)-pparton(4,4))**2
      do i=1,3
         s56=s56-(-pparton(1,i)-pparton(2,i)-pparton(3,i)-pparton(4,i))**2
      enddo
      call jade_hbb_fullp(pparton,ycut,jets,pincj,isub,s56)
      if(jets < 2) then
         write(6,*) 'PANIC!!! jade jets < 2 ',jets,s56
         call writeout(pparton)
         stop
!         pause
      endif

       if(jets .ne. 3) then
         passed=.false.
         return
      endif
      passed=.true.
      return
      ihard=1
      ihard2=0
      ptmax=pincj(1,4)
      ptmax2=0._dp
      do i=2,jets
         if(pincj(i,4).gt.ptmax) then
            ihard=i
            ptmax=pincj(i,4)
         endif
      enddo

      do i=1,jets
         if((pincj(i,4).gt.ptmax2).and.(i.ne.ihard)) then
            ptmax2=pincj(i,4)
            ihard2=i
         endif
      enddo
! 13   continue
      
      nn1(1:3)=pincj(ihard,1:3)
     &     /sqrt(pincj(ihard,1)**2+pincj(ihard,2)**2+pincj(ihard,3)**2)
      nn1(4)=1.0_dp

      nn2(1:3)=pincj(ihard2,1:3)
     &   /sqrt(pincj(ihard2,1)**2+pincj(ihard2,2)**2+pincj(ihard2,3)**2)
      nn2(4)=1.0_dp

!      ihard=5
!      ihard2=6
!      nn1(1:3)=pparton(ihard,1:3)
!     &     /sqrt(pparton(ihard,1)**2+pparton(ihard,2)**2+pparton(ihard,3)**2)
!      nn1(4)=1.0_dp

!      nn2(1:3)=pparton(ihard2,1:3)
!     &   /sqrt(pparton(ihard2,1)**2+pparton(ihard2,2)**2+pparton(ihard2,3)**2)
!      nn2(4)=1.0_dp


      tau=0._dp

    
      do j = ipp,npart+2-isub
      
      taua =pparton(j,4)*nn1(4)-pparton(j,1)*nn1(1)
     &     -pparton(j,2)*nn1(2)-pparton(j,3)*nn1(3)
      taub=pparton(j,4)*nn2(4)-pparton(j,1)*nn2(1)
     &    -pparton(j,2)*nn2(2)-pparton(j,3)*nn2(3)
      tau=tau+min(abs(taua),abs(taub))
 !     write(6,*) j,taua,taub,tau
      enddo

!      if((ihard2 < 5).and.(tau.lt.0.1_dp)) then 
!      call writeout(pparton)
!      write(6,*) scetjets,ipp
!      call writeout(pincj)
!      
!      write(6,*) ptjet
!      write(6,*) ptmax,ptmax2,ihard,ihard2

 !     pause
 

!      write(6,*) nn1(4),nn1(1),nn1(2),nn1(3)
!      write(6,*) nn2(4),nn2(1),nn2(2),nn2(3)
!      write(6,*) tau
      
!      pause
!      endif
c--- check to make sure no NaN
      if (tau .ne. tau) then
        call writeout(pparton)
        write(6,*) 'maketaucut.f:  tau=',tau
        stop
      endif

      if(tau.lt.taucut) then
         ttest(jj)=.false.
      else
         ttest(jj)=.true.
      endif

!      enddo
!      if(ttest(1).neqv.ttest(2)) then
!         write(6,*) 'found idff'
!         write(6,*) ttest(1),ttest(2)
!         pause
!      endif
         
c--- if tau is too small, implement the jettiness cut to remove the event
      if (tau .lt. taucut) return

c--- at this point cut has been passed
      passed=.true.

      return
      end
