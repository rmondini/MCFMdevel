!---- CW an implementation of the soft drop routine of
!---- 1402.2657

      subroutine softdrop_old(q,qsoft,Rmin,isub,ipow)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      real(dp) :: psoft(4)
      integer :: partoninj(mxpart,mxpart)
      integer npartsinjet,jetid,ncombine(mxpart)
      integer ix,ij
      include 'nf.f'
      include 'npart.f'
      include 'jetcuts.f'
      include 'jetlabel.f'
      include 'plabel.f'
      include 'is_functions_com.f'

      real(dp):: q(mxpart,4),qjet(mxpart,4),qfinal(mxpart,4),qsoft(mxpart,4)
      real(dp):: pt,Rmin,dijmin,dkmin,ayrap
      integer:: i,nu,iter,nmin1,nmin2,maxjet,nk,
     & ajet,jetindex(mxpart),isub,ipow
      logical:: jetmerge,failed,is_hadronic
      common/jetmerge/jetmerge
!$omp threadprivate(/jetmerge/)
      integer checklen_softd,SDjets
      real(dp) :: pt2,pt1,R0,z_cut,pj1(4),pj2(4)
      integer jetindexO(mxpart),firstJ
  
      integer, save :: nevent,neventmax
      character*4 mystring
      logical storenmin1,storenmin2,printout
      include 'first.f'

      if(first) then
         first=.false.
         nevent=1
         neventmax=99
      endif

      
      
      call writeout(q)
      write(6,*) ' ^^^^^ INPUT MOMENTUM ^^^^^'

 !     call writeout(q)
 !     write(6,*) '  ^^^^^^ INPUT MOMENTUM &^^^^^^^^'
      ncombine(:)=1
      npartsinjet=0
      jetid=0
!----- logic is partonsinjet( jetid,partonid 1-3) 3 max partons in jet at NNLO.
      partoninj(:,:)=0
      psoft(:)=zip
      qjet(:,:)=zip
!------stage 1, cluster with kt algorithm, but keep tabs on parton id which went into jet
      firstJ=-1
      jets=0
      maxjet=0
      jetmerge=.false.
      
      do i=1,mxpart
         do nu=1,4
            qsoft(i,nu)=0._dp
         enddo
      enddo

c---  pick out jets: note that we search to npart+2-isub, to get the
c---  number of particles right. Note that isub=0 for all calls except
c---  the dipole contributions, where isub=1.   
      do i=3,npart+2-isub
         if (is_hadronic(i)) then
            if(firstJ <0) firstJ=i
            maxjet=maxjet+1
            jetindex(maxjet)=i
            jetlabel(maxjet)=plabel(i)
            do nu=1,4
               qjet(maxjet,nu)=q(i,nu)
            enddo
         endif
      enddo

      
      jetindexO(:)=jetindex(:)
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

c---  skip clustering if we only have one parton  
      if (maxjet == 1) goto 2
      
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
c      write(*,*) 'Comparing pair (',nmin1,',',nmin2,') value of'
c      write(*,*) 'dijmin = ',dijmin,' with ',nk,' value of dk = ',dkmin

c---  step3: compare the two ...      
      if (dijmin < dkmin) then
c---  ... if we should combine, go ahead
c        write(*,*) 'Clustered ',nmin1,nmin2,jetindex(nmin1),jetindex(nmin2)
         jetmerge=.true.
         
        call combine(qjet,nmin1,nmin2)
        

!-----add nmin1 and nmin2 to partonsinj, provided they are not there already!
        ncombine(jets+1)=ncombine(jets+1)+1
        storenmin1=.true.
        storenmin2=.true.
c        write(6,*)  'JJ',jetindex(nmin1),jetindex(nmin2)
        if(ncombine(jets+1).ge.2) then
           ij=ncombine(jets+1) 
!---- check whether nmin1 and nmin2 are already stored. (note nmin1 is actually the only
!----- one which will cause problems
           do i=1,ij
              if(partoninj(jets+1,i)==jetindex(nmin1)) storenmin1=.false.
              if(partoninj(jets+1,i)==jetindex(nmin2)) storenmin2=.false.
           enddo
        endif
        if(storenmin1) then
           partoninj(jets+1,ix) = jetindex(nmin1)
           ix=ix+1
        endif
        if(storenmin2) then
           partoninj(jets+1,ix) = jetindex(nmin2)
           ix=ix+1
        endif

        
c---  combined object goes into nmin1, now shuffle nmin2 off the end 
        call swapjet_sd(qjet,jetindex,nmin2,maxjet)

        
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
!         write(6,*) ncombine(1:3)
!         write(6,*) jets,jetindex(jets),jetindex(nk)
         if(ncombine(jets)==1)  partoninj(jets,1)=jetindex(nk)
            
            
!---- reset ix
         ix=1
c         write(6,*) 'jet index before',jetindex(1:3)
c        write(*,*) 'Now swapping ',jets,' and ',nk
         call swapjet_sd(qjet,jetindex,jets,nk)
c         write(6,*) 'jet index after',jetindex(1:3)
c         pause
      endif

c--- in the next iteration we search for jets in pjet from iter+1...maxjet
c--- so if this condition isn't true then there's one jet left at maxjet

      if (iter < maxjet-1) goto 1
      
 2    continue      
      jets=jets+1
      if(ncombine(jets)==1)  partoninj(jets,1)=jetindex(nk)
      
c--- restore incoming partons
      do i=1,2
        do nu=1,4
          qsoft(i,nu)=q(i,nu)
        enddo
      enddo
c--- set all other momenta to zero and restore leptons
      do i=3,npart+2
        do nu=1,4
          qsoft(i,nu)=0._dp
          if (.not.(is_hadronic(i))) then
            qsoft(i,nu)=q(i,nu)
          endif
        enddo
      enddo


            
c----remove jets that are below the pT threhold or which lie outside
c----the observable rapidity region
     
c      write(*,*) 'AFTER CLUSTERING: Obtained ',jets,' jets'
     
c---  restore jets
      ajet=0
c      write(6,*) jets
!      write(6,*) ptjetmin,etajetmin,etajetmax
!      write(6,*) ncombine(1),ncombine(2),ncombine(3)
!      write(6,*) jetindex(1),partoninj(1,:)
!      write(6,*) jetindex(2),partoninj(2,:)
!      write(6,*) jetindex(3),partoninj(3,:)
!      call writeout(qjet)
!     pause
      printout=.false.
      do i=1,jets
        ajet=ajet+1

c        write(6,*) 'jet i = ',i,qjet(i,1:4)
 !        write(6,*) 'jet i =2 ',qjet(2,1:4)
 !        write(6,*) 'jet i =3 ',qjet(3,1:4)
         
         write(*,*) 'Jet ',i,'(',jetlabel(i),')',jetindex(i),jetindexO(i)
         write(*,*) 'pt: ',pt(i,qjet),' vs min. ',ptjetmin
        write(*,*) 'aeta: ',ayrap(i,qjet),' vs min. ',etajetmin
        write(*,*) 'aeta: ',ayrap(i,qjet),' vs max. ',etajetmax
        write(6,*) 'ncombine = ',ncombine(i)
        write(6,*) 'partons in j =',partoninj(i,1:3)
!        if ((pt(i,qjet) >= ptjetmin) .and.
!     &      (ayrap(i,qjet) >= etajetmin) .and.
!     &      (ayrap(i,qjet) <= etajetmax)) then 
!        do nu=1,4
!          qfinal(jetindex(ajet),nu)=qjet(i,nu)
!     enddo
        if(ncombine(ajet)==1) then
           
!---- single parton jet, nothing to do for soft drop
           qsoft(jetindexO(ajet),:)=q(partoninj(ajet,1),:)
       else
!----- call softdrop check
!---- jet has passed kt-clustering, so do soft drop on it and its constituents
          call softdrop_check(ajet,ncombine(ajet),partoninj,q,psoft,SDjets)
          qsoft(jetindexO(ajet),:)=psoft(:)
       endif

       if(pt(i,qjet).ne.pt(jetindexO(ajet),qsoft)) then
          printout=.true.
       endif
       
!       jetlabel(ajet)=jetlabel(i)
!      endif
      enddo


      
c--- if no jets are removed by eta and pt cuts, then jets=ajet
      if (ajet < jets) then
        do i=ajet+1,jets
          do nu=1,4
            qsoft(jetindex(i),nu)=0._dp
          enddo
        enddo
        jets=ajet
      endif

      nevent=nevent+1
      if(nevent < neventmax) then
         if(nevent < 10) write(mystring,11) nevent
         if((nevent>= 10).and.(nevent < 100)) write(mystring,12) nevent
         open(unit=10,status='unknown',file='inp_testevent_'//trim(mystring))
         open(unit=100,status='unknown',file='out_testevent_'//trim(mystring))
         write(6,*) jets
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
 11      format (i1)
 12      format (i2)
      endif

      call writeout(qsoft)
      write(6,*) ' ^^^^^ SOFT DROPPED MOMENTUM ^^^^^'
      write(6,*) jets
      if(nevent ==97) pause
 !     if(printout) pause
      return
      end


      subroutine softdrop_check(jetid,nij,pij,q,qsoft,SDjets) 
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      integer jetid,nij,pij(mxpart,mxpart),SDjets
      real(dp) :: q(mxpart,4),qsoft(4)
      integer checklen_softd
      integer i,pijtemp(mxpart),nijtemp,hangcheck
      real(dp) :: beta,z_cut,R_sd,R0,R12
      real(dp) :: pj1(4),pj2(4),pt1,pt2
      integer :: pj1id(mxpart),pj2id(mxpart)
      logical pass_softdrop
      
!------loop over partons in jet, each time we either pass soft drop
!------condition, or decompose the subjet further using C/A clustering.
      SDjets=0
      qsoft(:)=zip
      beta=2._dp
      z_cut=0.1_dp
      R0=1.0_dp
      
!-----CA cluster (if nij > 2), else trivial clustering

!---- just one parton in the jet, add it to qsoft and return
      if(nij ==1) then
         qsoft(:)=q(pij(jetid,1),:) 
         return
      endif

      nijtemp=nij
      pijtemp(:)=pij(jetid,:)
      hangcheck=1
 12   continue

 !     write(6,*) hangcheck,nijtemp,pijtemp(1:nijtemp)
      
      if(hangcheck > 10) goto 13
!---- CHECK nijtemp !=1 (if so can just return), note this is impossible on first pass nij =2 to get here, so 
      if(nijtemp==1) then
         qsoft(:)=q(pijtemp(1),:)
         return
      else                      !--- do CA splitting
         write(6,*) 'nijtemp = ',nijtemp
         write(6,*) 'pijtemp = ',pijtemp(1:3)
         call CA_clust_SoftDrop(q,pijtemp,nijtemp,pj1,pj2,pj1id,pj2id)

         write(6,*) 'inp pj1 = ',pj1
         write(6,*) 'partons in pj1',pj1id(1:3)
         write(6,*)  'inp pj2 = ',pj2
         write(6,*) 'partons in pj2',pj2id(1:3)
!---- note in the above routine pj1id always has more partons then pj2id, but pj1 may not be harder than pj2
!------calculate R12 and pts
         R12=r_sd(pj2,pj1)
         pt1=sqrt(pj1(1)**2+pj1(2)**2)
         pt2=sqrt(pj2(1)**2+pj2(2)**2)
         write(6,*) 'R12 between subjets =',R12
         write(6,*) 'passed softdtop ',pass_softdrop(pt1,pt2,z_cut,R12,R0,beta)
         write(6,*) z_cut,R0,R12,beta,pt1,pt2
!---- check soft drop
         if(pass_softdrop(pt1,pt2,z_cut,R12,R0,beta)) then
!---- keep combined jet object
            qsoft(:)=pj1(:)+pj2(:)
            return
         else
!-----failed, cluster again on harder subjet          
            if(pt1 > pt2) then
               nijtemp=checklen_softd(pj1id)
               pijtemp(:)=pj1id(:)
            else
               nijtemp=checklen_softd(pj2id)
               pijtemp(:)=pj2id(:)
            endif
         endif
      endif
     
      hangcheck=hangcheck+1
      goto 12
 13   continue 
!-----should never get here
      write(6,*) 'Error in Soft drop , did many clusterings',hangcheck 
      stop
      end
      


      subroutine CA_clust_SoftDrop(p,partid,npartid,pj1,pj2,pj1id,pj2id)
!-----routine returns pj1 and pj2, the last two CA clustered objects.
!----- routine clusters partons outlined in partid
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'

      real(dp) :: p(mxpart,4),pca(mxpart,4),pj1(4),pj2(4)
      integer partid(mxpart),i1,i2,n,i,npartid,j
      real(dp) :: dij,dijmin,d
      integer :: nj1(mxpart),nj2(mxpart),clustered
      integer :: pj1id(mxpart),pj2id(mxpart)
      integer maxjet,k,tempx(mxpart)
      integer checklen_softd,iter,itemp,jetid
      real(dp) ::  tmp,pt
     
      nj1(:)=0
      nj2(:)=0
      pj1id(:)=0
      pj2id(:)=0
      pj1(:)=zip
      pj2(:)=zip
      pca(:,:)=zip
      maxjet=npartid

      if(maxjet < 2) then
         write(6,*) 'Error in C/A routine for Soft drop',maxjet,' < 2'
         stop
      endif
      
      if(npartid==2) then
!-----nothing to do, since we only have two partons
         pj1id(1)=partid(1)
         pj1(:)=p(partid(1),:)
         pj2id(1)=partid(2)
         pj2(:)=p(partid(2),:)
         return
      endif
         
      do i=1,npartid
         pca(i,:)=p(partid(i),:)
      enddo
     
      
      clustered=0
      iter=0
 1    continue
      iter=iter+1
   
!---- stage one find dij min and compare to 1
      
      do i=iter,maxjet
         do j=i+1,maxjet
           
            d=dij(p,pca,i,j,0)
            if((i==1).and.(j==2)) then
               dijmin=d
               i1=i
               i2=j
            elseif( d < dijmin) then 
               dijmin=d
               i1=i
               i2=j
            endif
         enddo
      enddo

      write(6,*) 'i1,i2 = ',i1,i2

!---- check that indeed i1 and i2 do pass clustering condition
!---- 
      if(dijmin < 1) then
         call combine(pca,i1,i2)
!---- check i1 and i2 are currently stored in nj1, else start new jet nj2 (or add to nj1)
         call checkandadd_nj(nj1,nj2,partid(i1),partid(i2))
         do k=1,4
            tmp=pca(i2,k)
            pca(i2,k)=pca(maxjet,k)
            pca(maxjet,k)=tmp
            itemp=partid(i2)
            partid(i2)=partid(maxjet)
            partid(maxjet)=itemp
         enddo
         
 !        write(6,*) nj1
 !        write(6,*) nj2
 !        pause
         clustered=clustered+1
         maxjet=maxjet-1
         iter=iter-1
      else
         write(6,*) 'not sure what to do with that dijmin=',dijmin
         stop 
      endif

      
      if(maxjet==2) then

         
!----- We're done, have the last two C/A clusterings 
!---- combinatorics easier in main routine.

!-----options are : nj1 and nj2 are both filled by non-trivial clusterings
!-----or only one (nj1) has been filled and the second jet is a parton
         if(checklen_softd(nj2)==0) then
            do i=1,npartid
               if((pt(1,pca)==pt(partid(i),p)).or.(pt(2,pca)==pt(partid(i),p))) then
!-----we found the parton to go in nj2
               nj2(1)=partid(i)
            endif
         enddo
      endif
         
         if(checklen_softd(nj1) < checklen_softd(nj2)) then
            tempx(:)=nj1(:)
            itemp = i1
            nj1(:)=nj2(:)
            i1=i2
            i2=itemp
            nj2(:)=tempx(:)
         endif
!-define jets, should always be stored in j1 and j2 
         pj1(:)=pca(1,:)
         pj1id(:)=nj1(:)
         
         pj2(:)=pca(2,:)
         pj2id(:)=nj2(:)
         return
      else
!----- go round again
         goto 1
      endif
         
      
      return
      end
      

      subroutine checkandadd_nj(n1,n2,i1,i2)
      implicit none
      include 'types.f'
      include 'mxpart.f'

      integer n1(mxpart),n2(mxpart)
      integer i1,i2
      integer k,nonzeron1,nonzeron2
      logical i1inn1,i2inn1
      
!---- checks whether i1 and i2 are in n1 -
      
!---- if one is, add the other to the end of n1
!---- if neither is then, if both n1==n2==0 then add to n1
!-----else add to n2

      i1inn1=.false.
      i2inn1=.false.

      nonzeron1=0
      nonzeron2=0
      
      do k=1,mxpart
         if(n1(k).ne.0) nonzeron1=nonzeron1+1
         if(n2(k).ne.0) nonzeron2=nonzeron2+1
         if(n1(k)==i1) i1inn1=.true.
         if(n2(k)==i2) i2inn1=.true.
      enddo

      if(nonzeron1==0.and.nonzeron2==0) then
         n1(nonzeron1+1)=i1
         n1(nonzeron1+2)=i2
         return
      endif

      if(.not.(i1inn1).and..not.(i2inn1)) then
         n2(nonzeron2+1)=i1
         n2(nonzeron2+2)=i2
         return
      endif

      if(i1inn1.and. .not.(i2inn1)) then
         n1(nonzeron1+1)=i2
      endif
         
      if(i2inn1.and. .not.(i1inn1)) then
         n1(nonzeron1+1)=i1
      endif

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

      r1= (ei+p1(3))*(ej-p2(3))/
     &   ((ej+p2(3))*(ei-p1(3)))
!      r1= (p1(4)+p1(3))*(p2(4)-p2(3))/
!     &   ((p2(4)+p2(3))*(p1(4)-p1(3)))

      dely=0.5_dp*log(r1)

      r2= (p1(1)*p2(1)+p1(2)*p2(2))/sqrt(pti2*ptj2)
      if (r2 > +0.9999999_dp) r2=+1._dp
      if (r2 < -0.9999999_dp) r2=-1._dp
      delphi=acos(r2)
      write(6,*) 'delphi, dely',delphi,dely
      r_sd=sqrt(dely**2+delphi**2)
      
      return
      end
      

      function checklen_softd(nij)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      integer checklen_softd
      integer nij(mxpart),i

      checklen_softd=0
      do i=1,mxpart
         if(nij(i).ne.0) checklen_softd=checklen_softd+1
      enddo
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
!      write(6,*) 'soft D lhs = ',lhs
!      write(6,*) 'soft D rhs = ',rhs
!      write(6,*) pt1,pt2
!      write(6,*) R12
      
      pass_softdrop=.false.
      if(lhs.gt.rhs) pass_softdrop=.true.

      return
      end

  
      subroutine swapjet_sd(pjet,jetindex,i,j)
      implicit none
      include 'types.f'
c--- swaps jets i..j in pjet
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'jetlabel.f'
      integer:: i,j,k,jetindex(mxpart),itemp
      real(dp):: pjet(mxpart,4),tmp
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

      itemp=jetindex(i)
      jetindex(i)=jetindex(j)
      jetindex(j)=itemp
      return
      end
