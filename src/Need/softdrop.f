!---- improved version of softdrop routine
!---- based on Algorithm outlined in 1402.2657

!---- first kt cluster and build array npj(jetid,partonsinj)
!---- then go over jets, check pt/eta and then soft drop

!-----q= input partons
!-----qkt = output kt jets
!-----qsoft = output softdropped jets
      
      subroutine softdrop(q,qkt,qsoft,Rmin,isub,ipow)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'softdrop_var.f'
      real(dp) :: q(mxpart,4),qjet(mxpart,4),qsoft(mxpart,4)
      real(dp) :: qkt(mxpart,4),qjs(mxpart,4)
      real(dp) :: Rmin
      integer isub,ipow
      integer ns1,ns2,j
      integer :: npj(mxpart,mxpart)
      include 'nf.f'
      include 'npart.f'
      include 'jetcuts.f'
      include 'jetlabel.f'
      include 'plabel.f'
      include 'is_functions_com.f'
      real(dp):: pt,dijmin,dkmin,ayrap
      logical:: jetmerge,failed,is_hadronic
      common/jetmerge/jetmerge
!$omp threadprivate(/jetmerge/)
      integer ix,iter,nmin1,nmin2,nu,nk
      integer i,maxjet,ajet
      integer jetindex(mxpart),jarr(mxpart)
      integer nj1(mxpart),nj2(mxpart)
      real(dp):: pj1(4),pj2(4)
      real(dp) :: pjtemp(4),pt1,pt2
      real(dp) :: beta,R0,R12,zcut,r_sd
      logical pass_softdrop,singleton
      integer, save :: nevent,neventmax
      character*6 mystring
      
 !     include 'first.f'
      
 !     if(first) then
 !        first=.false.
 !        nevent=0
 !        neventmax=10001
 !     endif
!----- hard code soft drop params for now

      beta=2.0_dp
      R0=one
      zcut=0.1_dp
      
      npj(:,:)=0
      npjr(:,:)=0
      pinSDj(:,:)=0
      qsoft(:,:)=zip
      qjet(:,:)=zip
      qkt(:,:)=zip
      qjs(:,:)=zip
      firstjet=-1

      jets=0
      maxjet=0
      jetmerge=.false.
      
      do i=3,npart+2-isub
         if (is_hadronic(i)) then
            if(firstjet .lt. 0) firstjet=i
            maxjet=maxjet+1
            jetindex(maxjet)=i
            npj(maxjet,1)=i 
            jetlabel(maxjet)=plabel(i)
            do nu=1,4
               qjet(maxjet,nu)=q(i,nu)
            enddo
         endif
      enddo

c---  for no partons, just switch q into qfinal
      if (maxjet == 0) then
         do i=1,mxpart
            do nu=1,4
               qsoft(i,nu)=q(i,nu)
            enddo
         enddo
         jets=0
         return
      endif
      
c---  skip clustering if we only have one parton, identify it in npj
      if (maxjet == 1)         goto 2
            
      ix=1
      iter=0
c      write(6,*) 'init ncombine ',ncombine(1:3)
c--- loops through all the iterations of the algorithm      
    1 iter=iter+1
      
c--- step1: find (i,j) pair with lowest measure of all non-jets so far
      call findmind(q,qjet,iter,maxjet,dijmin,nmin1,nmin2,ipow)
 
c--- step2: find jet K with lowest Et
      call findminet(q,qjet,iter,maxjet,dkmin,nk,ipow)
      dkmin=dkmin*Rmin
!
c      write(*,*) 'dijmin = ',dijmin,' with ',nk,' value of dk = ',dkmin

c---  step3: compare the two ...      
      if (dijmin < dkmin) then
c---  ... if we should combine, go ahead
c        write(*,*) 'Clustered ',nmin1,nmin2,jetindex(nmin1),jetindex(nmin2)
         jetmerge=.true.
         
         call combine(qjet,nmin1,nmin2)
!---  combine npj
         call combine_npj(npj,nmin1,nmin2,jetindex)
         
c---  combined object goes into nmin1, now shuffle nmin2 off the end 
        call swapjet_sd(qjet,jetindex,npj,nmin2,maxjet)

        
        
        maxjet=maxjet-1
        iter=iter-1
        
      else
         
       
c---  ... we've finished a jet
         jets=jets+1
c     write(*,*) 'Now swapping ',jets,' and ',nk
!----- add nk to the npj if its not alreay in there 
         call swapjet_sd(qjet,jetindex,npj,jets,nk)
      endif

c--- in the next iteration we search for jets in pjet from iter+1...maxjet
c--- so if this condition isn't true then there's one jet left at maxjet

      if (iter < maxjet-1) then
         goto 1
      endif
      
 2    continue      
      jets=jets+1
!
!-----copy over qjet for soft drop
      qjs(:,:)=qjet(:,:)
      
!         call writeout(q)
!         write(6,*) 'before pt/yrap check',jets

 

      

!---- At this point qsoft is just kt clustered jets, now copy over and soft drop jets
 !     qkt(:,:)=qsoft(:,:)
      
 !     if(jets==0) return
!     write(6,*) jets
!      call writeout(q)
!      write(6,*) '^^^^ inp events'
!      call writeout(qjs)
!      write(6,*) 'before soft drop',jets

      do i=1,jets
         ix=1
         pjtemp(:)=zip
 !----- stage 1, j        
         jarr(:)=npj(i,:)
         
!         write(6,*) 'initial jarr ',jarr(1:5)
         call clustCA_sd(jarr,q,nj1,pj1,nj2,pj2,singleton)
 113     continue
         if(singleton) then
            qjs(i,:)=pj1(:)
            pinSDj(i,:)=nj1(:)
            cycle               !--- nothing to do for this jet
         endif
!-------have at least two partons in jet so soft drop it
         pt1=sqrt(pj1(1)**2+pj1(2)**2)
         pt2=sqrt(pj2(1)**2+pj2(2)**2)
         R12=r_sd(pj1,pj2)
         if(pass_softdrop(pt1,pt2,zcut,R12,R0,beta)) then
!------jet has passed, keep combined object
            pjtemp(:)=pj1(:)+pj2(:)
            qjs(i,:)=pjtemp(:)
            ns1 = count(nj1/=0)
            ns2 = count(nj2/=0)
            pinSDj(i,:)=nj1(:)
            do j=ns1+1,ns1+ns2
               pinSDj(i,j)=nj2(j-ns1)
            enddo
            cycle
         else
!-------soft drop harder subjet
            if(pt1 > pt2) then
               jarr(:)=nj1(:)
               call clustCA_sd(jarr,q,nj1,pj1,nj2,pj2,singleton)
            else
               jarr(:)=nj2(:)
               call clustCA_sd(jarr,q,nj1,pj1,nj2,pj2,singleton)
            endif
            ix=ix+1
            goto 113
         endif
      enddo

     
     
c--- restore incoming partons
      do i=1,2
         do nu=1,4
            qsoft(i,nu)=q(i,nu)
            qkt(i,nu)=q(i,nu)
         enddo
      enddo
c---  set all other momenta to zero and restore leptons
      do i=3,npart+2
         do nu=1,4
            qsoft(i,nu)=0._dp
            qkt(i,nu)=0._dp
            if (.not.(is_hadronic(i))) then
               qsoft(i,nu)=q(i,nu)
               qkt(i,nu)=q(i,nu)
            endif
         enddo
      enddo
!---- now check pt of jets to see if they are kept
      ajet=0
      do i=1,jets
        if ((pt(i,qjet) >= ptjetmin) .and.
     &       (ayrap(i,qjet) >= etajetmin) .and.
     &        (ayrap(i,qjet) <= etajetmax)) then
           ajet=ajet+1
           npjr(ajet,:)=npj(i,:)
           do nu=1,4
              qsoft(jetindex(ajet),nu)=qjs(i,nu)
              qkt(jetindex(ajet),nu)=qjet(i,nu)
           enddo
        endif
      enddo
       
!      pause
c--- if no jets are removed by eta and pt cuts, then jets=ajet
      if (ajet < jets) then
        do i=ajet+1,jets
          do nu=1,4
            qsoft(jetindex(i),nu)=0._dp
            qkt(jetindex(i),nu)=0._dp
          enddo
        enddo
        jets=ajet
      endif




!====== debug

!      write(*,*) 'phase-space point soft drop'

!      write(*,*)
!      do i=1,8
!         write(*,*) q(i,1),q(i,2),q(i,3),q(i,4)
!      enddo

!      write(*,*)
!      do i=1,8
!         write(*,*) qjet(i,1),qjet(i,2),qjet(i,3),qjet(i,4)
!      enddo

!      write(*,*)
!      do i=1,8
!         write(*,*) qkt(i,1),qkt(i,2),qkt(i,3),qkt(i,4)
!      enddo

!      write(*,*)
!      do i=1,8
!         write(*,*) qsoft(i,1),qsoft(i,2),qsoft(i,3),qsoft(i,4)
!      enddo

!      write(*,*)
!      write(*,*) 'npjr'
!      do i=1,4
!         write(*,*) npjr(i,1),npjr(i,2),npjr(i,3),npjr(i,4)
!      enddo

!      write(*,*)
!      write(*,*) 'n of jets in qkt = ',jets

!      write(*,*)
!      write(*,*) '*******************'

!      stop

!====== debug










      return
      
!------test SD
      if(jets > 0) then 
      nevent=nevent+1
      if(nevent < neventmax) then
         if(nevent < 10) write(mystring,11) nevent
         if((nevent>= 10).and.(nevent < 100)) write(mystring,12) nevent
         if((nevent>= 100).and.(nevent < 1000)) write(mystring,13) nevent
         if((nevent>= 1000).and.(nevent < 10000)) write(mystring,14) nevent
         
         open(unit=10,status='unknown',file='inp4_testevent_'//trim(mystring))
         open(unit=100,status='unknown',file='out4_testevent_'//trim(mystring))
         do i=3,npart+2-isub
            if(is_hadronic(i)) then 
               write(10,*) q(i,1),q(i,2),q(i,3),q(i,4)
            endif
         enddo
! debug, this is just true for VH example
         do i=5,jets+4              
            write(100,*) qsoft(i,1),qsoft(i,2),qsoft(i,3),qsoft(i,4)
         enddo
         close(10)
         close(100)
 11      format (i1)
 12      format (i2)
 13      format (i3)
 14      format (i4)
      endif
      

      endif
      

















!      if(nevent==2833) then
!         write(6,*) jets
!         call writeout(qsoft)
!         stop
!      endif
!      write(6,*) nevent
!      if(nevent==neventmax) stop

      
!------ test writeout 
 !     write(6,*) 'number of jets',jets
 !     write(6,*) 'partons in jets :'
 !     do i=1,jets
 !        ptemp(:)=zip
 !        write(6,*) 'jet ',i,' :',npj(i,1:4)
 !        do nu=1,mxpart
 !           if(npj(i,nu).ne.0) then
 !              ptemp(:)=ptemp(:)+q(npj(i,nu),:)
 !           endif
 !        enddo
!!         write(6,*) 'jet      ',i," : ",qsoft(jetindex(i),:) 
!!         write(6,*) 'rebuilt  ',i," : ",ptemp(:)
 !        if(abs(ptemp(4)/qsoft(jetindex(i),4)-1d0).gt.0.00001_dp) goto 1111
 !     enddo
!      pause
      
      return
      
 1111 call writeout(qsoft)
      write(6,*) '^^^^ CLUSTERED MOM^^^^'
      call writeout(q) 
      write(6,*) '^^^^ INPUT  MOM^^^^'
      end


      subroutine combine_npj(npj,nmin1,nmin2,jetindex)
!-----this routine adds jetindex(nmin2) to npj(nmin1)
!---- and deletes npj(nmin2) while shuffling everyone else up
      implicit none
      include 'mxpart.f'
      integer  npj(mxpart,mxpart),jetindex(mxpart)
      integer  nmin1,nmin2
      integer i,j,nonZ,atemp(mxpart),nonz2


      
      
      if(nmin2 < nmin1) then
         write(6,*) 'got my combos wrong'
         stop
      endif

      nonZ=1
      nonZ2=0
!---- find where to add
      do i=1,mxpart
         if(npj(nmin1,i).ne.0) nonZ=nonZ+1
         if(npj(nmin2,i).ne.0) nonZ2=nonZ2+1
      enddo

      do j=1,mxpart
         if(j.eq.nmin1) then
            do i=1,nonZ2
               npj(j,nonZ-1+i) = npj(nmin2,i)
            enddo
         elseif(j.eq.nmin2) then
            npj(j,:)=0
         endif
      enddo
          
      return 
      end    

      subroutine swapjet_sd(pjet,jetindex,npj,i,j)
      implicit none
      include 'types.f'
c--- swaps jets i..j in pjet
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'jetlabel.f'
      integer:: i,j,k,jetindex(mxpart),itemp
      real(dp):: pjet(mxpart,4),tmp
      integer  npj(mxpart,mxpart),ntmp(mxpart)
      character*2 chartmp
 
c--- escape if we're trying to swap the same jets
      if (i == j) return

      do k=1,4
        tmp=pjet(i,k)
        pjet(i,k)=pjet(j,k)
        pjet(j,k)=tmp
      enddo
 
      chartmp=jetlabel(i)
      jetlabel(i)=jetlabel(j)
      jetlabel(j)=chartmp

!      itemp=jetindex(i)
!      jetindex(i)=jetindex(j)
!      jetindex(j)=itemp


      ntmp(:)=npj(i,:)
      npj(i,:)=npj(j,:)
      npj(j,:)=ntmp(:)
      
      
      return
      end

      function pass_softdrop(pt1,pt2,zcut,R12,R0,beta)
      implicit none
      include 'types.f'
      logical pass_softdrop

      real(dp) :: pt1,pt2,zcut,R12,R0,beta
      real(dp) :: lhs,rhs

      lhs=min(pt1,pt2)/(pt1+pt2)

      rhs=zcut*(R12/R0)**beta
      
      pass_softdrop=.false.
      if(lhs.gt.rhs) pass_softdrop=.true.

      return
      end

        function r_sd(p1,p2)
      implicit none
      include 'types.f'
      real(dp):: r_sd
c---- calculate the jets separation between p(i) and p(j)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      real(dp):: p1(4),p2(4),r1,r2,dely,delphi,ei,ej,pti2,ptj2
      integer:: i,j

      
      
      pti2=p1(1)**2+p1(2)**2
      ptj2=p2(1)**2+p2(2)**2

      ei=sqrt(pti2+p1(3)**2)
      ej=sqrt(ptj2+p2(3)**2)

      ei=p1(4)
      ej=p2(4)
      
      r1= (ei+p1(3))*(ej-p2(3))/
     &   ((ej+p2(3))*(ei-p1(3)))
!      r1= (p1(4)+p1(3))*(p2(4)-p2(3))/
!     &   ((p2(4)+p2(3))*(p1(4)-p1(3)))

      dely=0.5_dp*log(r1)

      r2= (p1(1)*p2(1)+p1(2)*p2(2))/sqrt(pti2*ptj2)
      if (r2 > +0.9999999_dp) r2=+1._dp
      if (r2 < -0.9999999_dp) r2=-1._dp
      delphi=acos(r2)
      r_sd=sqrt(dely**2+delphi**2)
      
      return
      end
