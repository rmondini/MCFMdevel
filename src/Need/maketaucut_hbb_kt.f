!---  HEALTH WARNING: THIS ROUTINE IS DESIGNED FOR ALL FSR FROM HIGGS!!
!---- CW March 17
      subroutine maketaucut_hbb(pparton,isub,itag,passed)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'npart.f'
      include 'taucut.f'
      include 'plabel.f'
      include 'first.f'
      integer itest
      logical passed
      integer isub,j
      real(dp) :: pparton(mxpart,4),pjetboost(mxpart,4)
      real(dp) :: ppartonboost(mxpart,4),ptest(4)
      real(dp) :: phiggs(4),lambda,beta(3),betasq
      real(dp) :: ptjet(mxpart),ptmax,ptmax2
      integer ihard,ihard2,npp
      real(dp) :: nn1(4),nn2(4)
      real(dp) :: pfinal(mxpart,4),phad(mxpart,4),phadboost(mxpart,4)
      logical is_hadronic
      integer i,ipp,itag
      real(dp) :: s56,gamma,tau,taua,taub
      integer jets,maxjet,iter,jetindex(mxpart)
      integer nmin1,nmin2
      real(dp) :: dijmin,dkmin,Rmin
      integer ipow,nk
      integer ajet,nu
      logical jetmerge_hbb
      common/jetmerge_hbb/jetmerge_hbb
!$omp threadprivate(/jetmerge_hbb/)


      
      passed=.false.
      jets=0
      maxjet=0
      ipow=-1
      jetindex(:)=0
      Rmin=0.5_dp
      jetmerge_hbb=.false.
      
!-----stage 0 - get Higgs momentum
      if(itag==34) then
         phiggs(:) = -pparton(1,:)-pparton(2,:)-pparton(3,:)-pparton(4,:)
      else
         write(6,*) 'Unrecognized itag in Higgs=>bb taucut routine'
         write(6,*) 'itag = ',itag
         stop
      endif

      s56=phiggs(4)**2-phiggs(1)**2-phiggs(2)**2-phiggs(3)**2

!-----define beta for boost
      beta(1)=phiggs(1)/phiggs(4)
      beta(2)=phiggs(2)/phiggs(4)
      beta(3)=phiggs(3)/phiggs(4)

      phad(:,:)=zip
      phadboost(:,:)=zip
      pjetboost(:,:)=zip
      i=0      
!-----stage 1 - identify hadrons in pparton and boost to higgs rest frame
      do j=3,npart+2-isub
         if(is_hadronic(j)) then
            i=i+1
            maxjet=maxjet+1
            jetindex(maxjet)=i
            phad(i,:)=pparton(j,:)
         endif
         call myboost(phad(i,:),phadboost(i,:),beta)
         pjetboost(i,:)=phadboost(i,:)
      enddo
      npp=i
    
 
!      if(npart==5.and.maxjet==2) then
!         write(6,*) '**********************'
!         call writeout(pparton)
!         call writeout(phadboost)
!         call writeout(pjetboost)
!         write(6,*) npart,isub
!         write(6,*) jetindex
!         stop
!      endif

      
!      call writeout(pparton)
!      call writeout(phad)
!      call writeout(phadboost)
!      write(6,*) beta
!      write(6,*) s56
!      stop
      
!-----stage 2 re-cluster in Higgs rest-frame and idenfity two leading-jets
!---- this is a local anti-kt clustering with paramters defined as inputs
!---- therefore has no bearing on any other jet-couting which may (or may not)
!---- happen
      
c--- skip clustering if we only have one parton  
!      if (maxjet.le.1) goto 2

!----- if maxjet =2 then we have two back-to-back partons, no need to check tau
      if(npp==2) then
         passed=.true.
         return
      endif

      goto 12
      iter=0
c--- loops through all the iterations of the algorithm      
    1 iter=iter+1

c      write(*,*) 'iter ',iter
c      write(*,*) 'jets ',jets
c      write(*,*) 'maxjet ',maxjet

c--- step1: find (i,j) pair with lowest measure of all non-jets so far
      call findmind(phadboost,pjetboost,iter
     &     ,maxjet,dijmin,nmin1,nmin2,ipow)

c--- step2: find jet K with lowest Et
      call findminet(phadboost,pjetboost,iter,maxjet,dkmin,nk,ipow)
      dkmin=dkmin*Rmin

c--- step3: compare the two ...      
      if (dijmin < dkmin) then
c---  ... if we should combine, go ahead
c        write(*,*) 'Clustered ',nmin1,nmin2
         jetmerge_hbb=.true.
    
        call combine(pjetboost,nmin1,nmin2)
c--- combined object goes into nmin1, now shuffle nmin2 off the end 
        call swapjet(pjetboost,jetindex,nmin2,maxjet)        
        maxjet=maxjet-1
        iter=iter-1
c        do i=1,maxjet
c          do j=1,4
c            write(*,*) 'qjet(',i,',',nu,') = ',qjet(i,nu)
c          enddo
c        enddo
      else
c---  ... we've finished a jet
        jets=jets+1
c        write(*,*) 'Now swapping ',jets,' and ',nk
        call swapjet(pjetboost,jetindex,jets,nk)
      endif

c--- in the next iteration we search for jets in pjet from iter+1...maxjet
c--- so if this condition isn't true then there's one jet left at maxjet

      if (iter < maxjet-1) goto 1
      jets=jets+1


      ajet=0
      pfinal(:,:)=zip
      do i=1,jets
         ajet=ajet+1
         do nu=1,4
            pfinal(jetindex(ajet),nu)=pjetboost(i,nu)
         enddo
      enddo

      ptjet(:)=0._dp
      
      do i=1,jets
         ptjet(i)=pfinal(i,4)
      enddo
      ptmax=0._dp
      ptmax2=0._dp
      do i=1,jets
         if(ptjet(i).gt.ptmax) then
            ihard=i
            ptmax=ptjet(i)
         endif
      enddo

      do i=1,jets
         if((ptjet(i).gt.ptmax2).and.(i.ne.ihard)) then
            ihard2=i
            ptmax2=ptjet(i)
         endif
      enddo
            
      
 !     if(jets==2) then
 !        ihard=1
 !        ihard2=2
 !        ptmax2=ptmax
!     endif

 12   continue
      call jade_hbb(phadboost,0.01_dp,itest,pfinal,npp)
!      write(6,*) itest
! 12   call jade_hbb(phadboost,npart,0.01_dp,jets,isub,pfinal)
! 12   continue
      
!      jets=itest
!      write(6,*) itest,jets
!      call writeout(pfinal)
!      pause
      ptmax=0._dp
      ptmax2=0._dp
      do i=1,itest
         if(pfinal(i,4).gt.ptmax) then
            ihard=i
            ptmax=pfinal(i,4)
         endif
      enddo

      if(itest > 1) then 
      do i=1,itest
         if((pfinal(i,4).gt.ptmax2).and.(i.ne.ihard)) then
            ihard2=i
            ptmax2=pfinal(i,4)
         endif
      enddo
      else
         ihard2=ihard
      endif

      
      nn1(1:3)=pfinal(ihard,1:3)
     &   /sqrt(pfinal(ihard,1)**2+pfinal(ihard,2)**2+pfinal(ihard,3)**2)
      nn1(4)=1.0_dp

      nn2(1:3)=pfinal(ihard2,1:3)
     &/sqrt(pfinal(ihard2,1)**2+pfinal(ihard2,2)**2+pfinal(ihard2,3)**2)
      nn2(4)=1.0_dp


    
!-----stage 3 compute tau and comapre to taucut
      tau=0._dp
      do j = 1,npp     
      taua =phadboost(j,4)*nn1(4)-phadboost(j,1)*nn1(1)
     &     -phadboost(j,2)*nn1(2)-phadboost(j,3)*nn1(3)
      taub=phadboost(j,4)*nn2(4)-phadboost(j,1)*nn2(1)
     &    -phadboost(j,2)*nn2(2)-phadboost(j,3)*nn2(3)
      tau=tau+min(abs(taua),abs(taub))

      if((npp >=3).and.(jets==2)) then 
!         call writeout(pparton)
!         call writeout(phadboost)
!         call writeout(pfinal)
!         write(6,*) s56,jets,npp,ihard,ihard2
 !        write(6,*) tau,taucut
!         write(6,*) j,taua,taub,tau
      endif
      enddo


!       if((npp >=3).and.(jets==2)) pause
!      write(6,*) nn1(4),nn1(1),nn1(2),nn1(3)
!      write(6,*) nn2(4),nn2(1),nn2(2),nn2(3)
!      write(6,*) taua,taub

!      pause
c--- check to make sure no NaN
      if (tau .ne. tau) then
         passed=.false.
         return 
        call writeout(pparton)
        write(6,*) 'maketaucut.f:  tau=',tau
        stop
      endif

c--- if tau is too small, implement the jettiness cut to remove the event
      if (tau .lt. taucut) then
         passed=.false.
         return
!     call writeout(pparton)
!         call writeout(phadboost)
!         write(6,*) nn1
!         write(6,*) nn2
!         write(6,*) jets
!         write(6,*) npp,tau,taub,tau
!         pause
      endif
      
c--- at this point cut has been passed
      passed=.true.

    

!---- stage 4 return
!      write(6,*) 'not here'
      return
 2    write(6,*) 'not enough partons for H=>bb decay routine'
      write(6,*) 'found jets = ',maxjet
      end


      subroutine myboost(pin,pout,beta)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      real(dp) :: pin(4),pout(4),gamma,beta(3),betasq

      pout(:)=zip
      betasq=beta(1)**2+beta(2)**2+beta(3)**2
      gamma=one/sqrt(one-betasq)
      pout(4)=gamma*(pin(4)-beta(3)*pin(3)-beta(2)*pin(2)-beta(1)*pin(1))

      pout(1)= (one + ((-one + gamma)*beta(1)**2)/betasq)*pin(1) + 
     -  ((-one + gamma)*beta(1)*beta(2)*pin(2))/betasq + 
     -  ((-one + gamma)*beta(1)*beta(3)*pin(3))/betasq - 
     -     gamma*beta(1)*pin(4)

      pout(2)=((-one + gamma)*beta(1)*beta(2)*pin(1))/betasq + 
     -  (one + ((-one + gamma)*beta(2)**2)/betasq)*pin(2) + 
     -  ((-one + gamma)*beta(2)*beta(3)*pin(3))/betasq - 
     -     gamma*beta(2)*pin(4)

      pout(3)= ((-one + gamma)*beta(1)*beta(3)*pin(1))/betasq + 
     -  ((-one + gamma)*beta(2)*beta(3)*pin(2))/betasq + 
     -  (one + ((-one + gamma)*beta(3)**2)/betasq)*pin(3) - 
     -  gamma*beta(3)*pin(4)

      
      return
      end
