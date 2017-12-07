!--- CW Nov 17
      subroutine clustCA_sd(npi,p,nj1,pj1,nj2,pj2,singleton)
      implicit none
!---- routine for Soft drop, takes in npi (parton info for jet construction)
!-----p (parton level moemnta)
!---  returns n1 p1 and n2 p2 which are the last two stages of CA clustering
!---- n1 = partons in jet 1, n2 = partons in p2

      
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      integer npi(mxpart)
      real(dp) :: p(mxpart,4),pca(mxpart,4),pj1(4),pj2(4),pji(4)
      real(dp) :: dij,dijmin,d
      integer :: nj1(mxpart),nj2(mxpart),ntemp(mxpart)
      integer maxjet,len_n
      integer checklen_softd,iter,itemp,jetid
      real(dp) ::  tmp,pt
      logical singleton
      integer  :: clust1(mxpart,mxpart)
      integer i,i1,i2,j
      integer, save :: nevent,neventmax
      character*6 mystring
      logical is_hadronic
      include 'npart.f'
      include 'first.f'
      
       if(first) then
          first=.false.
          nevent=0
          neventmax=10001
       endif

      nj1(:)=0
      nj2(:)=0
      pj1(:)=zip
      pj2(:)=zip
      clust1(:,:)=0
      pca(:,:)=zip
      
      singleton=.false. 
      len_n=checklen_softd(npi)

      
      if(len_n ==1) then
!---- no clustering to do, just fill n1 and return
         pj1(:)=p(npi(1),:)
         nj1(:)=npi(:)
         singleton=.true. 
         return
      elseif(len_n==2) then
!-----trivial clustering, return both
         pj1(:)=p(npi(1),:)
         pj2(:)=p(npi(2),:)
         nj1(1)=npi(1)
         nj2(1)=npi(2)
         return
      endif

      pji(:)=zip
!---- have some work to do, at least three partons,
!---- copy over relevant partons into pca first and define maxjet
      do i=1,len_n
         pji(:)=pji(:)+p(npi(i),:)
         clust1(i,1)=npi(i)
         pca(i,:)=p(npi(i),:)
      enddo
      
      iter=0
      maxjet=len_n
      
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
     
!---- check that indeed i1 and i2 do pass clustering condition
!---- 
      if(dijmin < 1) then
         call combine(pca,i1,i2)
!===== merge in array too 
         call combine_npj(clust1,i1,i2,npi)
         call swapjet_sd(pca,npi,clust1,i2,maxjet)
!====== move combined object into i1
         maxjet=maxjet-1
         iter=iter-1
      else
!-----return event as singleton for now, keep input and output to compare to gav's code
         nevent=nevent+1
         if(nevent < neventmax) then
            if(nevent < 10) write(mystring,11) nevent
            if((nevent>= 10).and.(nevent < 100)) write(mystring,12) nevent
            if((nevent>= 100).and.(nevent < 1000)) write(mystring,13) nevent
            if((nevent>= 1000).and.(nevent < 10000)) write(mystring,14) nevent
            
            open(unit=10,status='unknown',file='Iwe2_testevent_'//trim(mystring))
            open(unit=100,status='unknown',file='Owe2_testevent_'//trim(mystring))
            do i=3,npart+2
               if(is_hadronic(i)) then 
                write(100,*) p(i,1),p(i,2),p(i,3),p(i,4)
              endif
            enddo
            write(10,*) pji(1),pji(2),pji(3),pji(4)
            close(10)
            close(100)
 11         format (i1)
 12         format (i2)
 13         format (i3)
 14         format (i4)
         endif
         pj1(:)=p(npi(1),:)
         nj1(:)=npi(:)
         singleton=.true. 
         return
      endif

      
      if(maxjet==2) then
!-----are done, since we want the last two clustered bits returned
         nj1(:)=clust1(1,:)
         nj2(:)=clust1(2,:)
         pj1(:)=pca(1,:)
         pj2(:)=pca(2,:)         
      else
!------go around again
         goto 1
      endif

      

!-------check output
!      write(6,*) 'partons in jet 1',nj1(1:5)
!      write(6,*) 'partons in jet 2',nj2(1:5)
!      write(6,*) 'momentum jet 1 :',pj1(:)
!      write(6,*) 'momentum jet 2 :',pj2(:)

      
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
