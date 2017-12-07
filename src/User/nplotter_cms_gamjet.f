      subroutine nplotter_cms_gamjet(p,wt,wt2,switch,nd)
      implicit none
      include 'types.f'
c--- Variable passed in to this routine:
c
c---      p:  4-momenta of particles in the format p(i,4)
c---          with the particles numbered according to the input file
c---          and components labelled by (px,py,pz,E)
c
c---     wt:  weight of this event
c
c---    wt2:  weight^2 of this event
c
c--- switch:  an integer:: equal to 0 or 1, depending on the type of event
c---                0  --> lowest order, virtual or real radiation
c---                1  --> counterterm for real radiation
!-- nd determines whether we are analysing a photon dipole and hence have
!---to rescale accordingly
      
      include 'vegas_common.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'histo.f'
      include 'jetlabel.f'
      include 'frag.f'
      include 'phot_dip.f'
      include 'outputflags.f'
      include 'kprocess.f' 
      real(dp):: p(mxpart,4),wt,wt2,ptV2j
      real(dp):: yrap,pt,r,mjj,logmin,logmax,logbin,pttwo,yraptwo
      real(dp):: yV,y4,y5,ptV,pt4,pt5,r45,r35,r34,s34,m34
      real(dp):: ptj,yj,phi_gj,mgj,cos_star,r2,yg,etarap,ptmax
      integer:: switch,n,nplotmax,jtag
      integer tag
      integer:: nd,i
      logical, save::first=.true.
      common/nplotmax/nplotmax
ccccc !$omp threadprivate(first,/nplotmax/)
      real(dp) :: ptreg1,ptreg2,ptreg3,ptreg4,ayj
  
************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************

      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag=tagbook
        yV=1d3
        y4=1d3
        ptV=1d3
        pt4=1d3
        ptreg1=-1d8
        ptreg2=-1d8
        ptreg3=-1d8
        ptreg4=-1d8
        ptV2j = -1d8

c---Intiailise jet 
        y5=1d3
        pt5=1d3
c--- If re5 is not changed by the NLO value, it will be out of
c--- the plotting range
        r34=1d3
        r35=1d3
        r45=1d3
        mjj=1d7 
        jets=2
        goto 99
      else
c--- Add event in histograms
        tag=tagplot
      endif

************************************************************************
*                                                                      *
*     DEFINITIONS OF QUANTITIES TO PLOT                                *
*                                                                      *
************************************************************************
    
!      ptmax=zip
!      jtag=-1
!---- find hardest jet
!      do i = 4,3+jets
!        if(pt(i,p) > ptmax) then
!          ptmax=pt(i,p)
!          jtag=i
!        endif
!      enddo
!      ptj=pt(jtag,p)
!      yj=yrap(jtag,p)
!      ayj=abs(yj)
!      yg=yrap(3,p)

!---- invariant mass of photon + jet 
!      mgj=(p(3,4)+p(jtag,4))**2-(p(3,1)+p(jtag,1))**2
!     &   -(p(3,2)+p(jtag,2))**2-(p(3,3)+p(jtag,3))**2
!      enddo
!      mgj=sqrt(mgj)
!      cos_star=tanh((yg-yj)/2._dp)
!
!      r2= (p(3,1)*p(jtag,1)+p(3,2)*p(jtag,2))
!     &     /(ptj*ptV)
!      if (r2 > +0.9999999_dp) r2=+1._dp
!      if (r2 < -0.9999999_dp) r2=-1._dp
!      
!      phi_gj=acos(r2)

!--- vector boson pt and rapidity (accounting also for Z+jet case)
      if((kcase==kZ_1jet).or.(kcase==kZ_2jet)) then 
         ptV=pttwo(3,4,p)
         yV=yraptwo(3,4,p)
      else
        ptV=pt(3,p)
        yV=etarap(3,p)
      endif
         
************************************************************************
*                                                                      *
*     FILL HISTOGRAMS                                                  *
*                                                                      *
************************************************************************

c--- Call histogram routines
   99 continue

c--- Book and fill ntuple if that option is set, remembering to divide
c--- by # of iterations now that is handled at end for regular histograms
      if (creatent .eqv. .true.) then
        call bookfill(tag,p,wt/real(itmx,dp))  
      endif

c--- "n" will count the number of histograms
      n=nextnplot              

c--- Syntax of "bookplot" routine is:
c
c---   call bookplot(n,tag,titlex,var,wt,wt2,xmin,xmax,dx,llplot)
c
c---        n:  internal number of histogram
c---      tag:  "book" to initialize histogram, "plot" to fill
c---   titlex:  title of histogram
c---      var:  value of quantity being plotted
c---       wt:  weight of this event (passed in)
c---      wt2:  weight of this event (passed in)
c---     xmin:  lowest value to bin
c---     xmax:  highest value to bin
c---       dx:  bin width
c---   llplot:  equal to "lin"/"log" for linear/log scale   

!====== pt(V) in 20 equal bins of log10(ptV) from 100 to 800 GeV
      logmin=2._dp ! log10(100._dp)
      logmax=log10(800._dp) ! log10(800._dp)
      logbin=(logmax-logmin)/20._dp 
      
      call bookplot(n,tag,'l10 pt(V) for rat',log10(ptV),wt,wt2
     &     ,logmin,logmax,logbin,'log')
      call ebookplot(n,tag,log10(ptV),wt)
      n=n+1
!      call bookplot(n,tag,'pt(V) for rat',ptV,wt,wt2,100._dp,800._dp
!     &,10._dp,'log')
!      n=n+1

      
!-----first 19 bins can be obtianed from above, 
      logmax=1093._dp
      logmin=721._dp
      logbin=logmax-logmin
      call bookplot(n,tag,'l10 pt(V) for fig2 a',ptV,wt,wt2
     &     ,logmin,logmax,logbin,'log')
      call ebookplot(n,tag,ptV,wt)
      n=n+1
!      call bookplot(n,tag,'pt(V) for fig 2 a',ptV,wt,wt2,100._dp,721._dp
!     &,10._dp,'log')
!      n=n+1
      

!-----nj >= 2 pt spectrum

      if(jets>=2) then
         ptV2j=ptV
      else
         ptV2j=zip
      endif

      logmin=2._dp              ! log10(100._dp)
      logmax=log10(800._dp)     ! log10(800._dp)
      logbin=(logmax-logmin)/20._dp 
      call bookplot(n,tag,'l10 pt(V) nj >= 2',log10(ptV2j),wt,wt2
     &     ,logmin,logmax,logbin,'log')
      call ebookplot(n,tag,log10(ptV2j),wt)
      n=n+1
!      call bookplot(n,tag,'pt(V) for fig 2 a',ptV2j,wt,wt2,100._dp,721._dp
!     &,10._dp,'log')
!      n=n+1
      logmax=1093._dp
      logmin=721._dp
      logbin=logmax-logmin
!-----overflow bin to 1 TeV
      call bookplot(n,tag,'l10 pt(V) nj >=2 b',ptV2j,wt,wt2
     &     ,logmin,logmax,logbin,'log')
      call ebookplot(n,tag,ptV2j,wt)
      n=n+1
!      call bookplot(n,tag,'pt(V) for fig 2 a',ptV2j,wt,wt2,721._dp,1000._dp
!     &,10._dp,'log')
!      n=n+1
      

************************************************************************
*                                                                      *
*     FINAL BOOKKEEPING                                                *
*                                                                      *
************************************************************************

c--- We have over-counted the number of histograms by 1 at this point
      n=n-1

c--- Ensure the built-in maximum number of histograms is not exceeded    
      call checkmaxhisto(n)

c--- Set the maximum number of plots, on the first call
      if (first) then
        first=.false.
        nplotmax=n
      endif
 
     
      return 
      end
