!==== C.W Feb 16, ATLAS cuts for Higgs studies, see for instance
!==== 1407.4222
      function gencuts_ATLAS_gaga(p)
      implicit none
      include 'types.f'
      logical gencuts_ATLAS_gaga
      include 'constants.f'
      include 'mxpart.f'
      include 'first.f' 
      real(dp) :: p(mxpart,4),mgg,ptgaga,pth,pts,ygh,ygs,ygg
      common/gaga_observables/mgg,ptgaga,pth,pts,ygh,ygs,ygg
!$omp threadprivate(/gaga_observables/)
      integer i
      real(dp):: ptfrac1,ptfrac2
      real(dp):: pt,pttwo,yrap,yraptwo
!-----observables for BSM interest
      real(dp):: gaga_highm_ptgaga(0:3),gaga_highm_ptgh(0:3)
      real(dp):: gaga_highm_ptgs(0:3),gaga_highm_ygg(0:3)
      real(dp):: gaga_highm_yh(0:3),gaga_highm_ys(0:3)
      integer itag
      common/gaga_highmgg_obs/gaga_highm_ptgaga,gaga_highm_ptgh,
     &     gaga_highm_ptgs,gaga_highm_ygg,gaga_highm_yh,gaga_highm_ys
!$omp threadprivate(/gaga_highmgg_obs/)

      
!=====initalize out of bounds
      
      gaga_highm_ptgaga(:)=-1._dp
      gaga_highm_ptgh(:)=-1._dp
      gaga_highm_ptgs(:)=-1._dp
      gaga_highm_ygg(:)=-10000._dp
      gaga_highm_yh(:)=-10000._dp
      gaga_highm_ys(:)=-10000._dp
      
      gencuts_ATLAS_gaga=.false.

      ptfrac1=0.35_dp
      ptfrac2=0.25_dp
!=====basic photon cuts from input file, this routine calculates
!==== mgaga and checks that

!=====pt(ga,hard)/mgaga > ptcut1, and pt(ga,soft)/mgaga > ptcut2

!=====initalize histo observables
      pts=-1._dp
      pth=-1._dp
      mgg=-1._dp
      ptgaga=-1._dp
      ygh=-1000._dp
      ygs=-1000._dp
      ygg=-1000._dp
      
      
      if(first) then
         first=.false.
         write(6,*) '******** Additional photon cuts applied *****'
         write(6,*) '* pt(hard)/m34 > ',ptfrac1,'              *'
         write(6,*) '* pt(soft)/m34 > ',ptfrac2,'              *'
         write(6,*) '*********************************************'
         write(6,*) '***** ATLAS Rapidity conditions applied *****'
         write(6,*) '* |eta(gamma) | <  ',2.37_dp,'              *'
         write(6,*) '* ATLAS Crack (1.37, 1.56) excluded         *'
         write(6,*) '*********************************************'
      endif
      
      
      mgg=zip
      do i=1,3
         mgg=mgg-(p(3,i)+p(4,i))**2
      enddo
      i=4
      mgg=mgg+(p(3,i)+p(4,i))**2
      mgg=sqrt(mgg)
      
!===== routine also calculates quantities for histogram 
      if(pt(3,p).gt.pt(4,p)) then
         pth=pt(3,p)
         pts=pt(4,p)
         ygh=yrap(3,p)
         ygs=yrap(4,p)
         
      else
         pts=pt(3,p)
         pth=pt(4,p)
         ygh=yrap(4,p)
         ygs=yrap(3,p)
      endif

      if((pth/mgg < ptfrac1).or.(pts/mgg < ptfrac2)) then
         goto 101
      endif

!-----rapidity and crack measurements
      do i=3,4
         if(abs(yrap(i,p)) > 2.37_dp) goto 101
         if((abs(yrap(i,p)) > 1.37_dp).and.(abs(yrap(i,p))  < 1.56_dp))
     &        goto 101
      enddo
         
!======calculate stuff for histos
      ptgaga=pttwo(3,4,p)
      ygg=yraptwo(3,4,p) 
      

      itag=2
!========now organize for mgg histos
      if((mgg >= 120._dp).and.(mgg  <= 130._dp)) then
         itag=1
      endif
      

      gaga_highm_ptgaga(0)=ptgaga
      gaga_highm_ptgaga(itag)=ptgaga
      
      gaga_highm_ptgh(0)=pth
      gaga_highm_ptgh(itag)=pth
      gaga_highm_ptgs(0)=pts
      gaga_highm_ptgs(itag)=pts

      gaga_highm_ygg(0)=ygg
      gaga_highm_ygg(itag)=ygg
      
      gaga_highm_yh(0)=ygh
      gaga_highm_yh(itag)=ygh

      gaga_highm_ys(0)=ygs
      gaga_highm_ys(itag)=ygs

      
      return
!=====zero histos if fail
 101  continue
      gencuts_ATLAS_gaga=.true.
      pts=-1._dp
      pth=-1._dp
      mgg=-1._dp
      ptgaga=-1._dp
      ygh=-1000._dp
      ygs=-1000._dp
      ygg=-1000._dp
      
      
      return 
      end
      
