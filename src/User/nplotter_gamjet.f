      subroutine nplotter_gamjet(p,wt,wt2,switch,nd)
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
      real(dp):: p(mxpart,4),wt,wt2
      real(dp):: yrap,pt,r,mjj
      real(dp):: y3,y4,y5,pt3,pt4,pt5,r45,r35,r34,s34,m34
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
        y3=1d3
        y4=1d3
        pt3=1d3
        pt4=1d3
        ptreg1=-1d8
         ptreg2=-1d8
         ptreg3=-1d8
         ptreg4=-1d8

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
c--- Photons order based on pt
    
      pt3=pt(3,p)
      y3=etarap(3,p)

      ptmax=zip
      jtag=-1
!---- find hardest jet
      do i = 4,3+jets
         if(pt(i,p).gt.ptmax) then
            ptmax=pt(i,p)
            jtag=i
         endif
      enddo
      ptj=pt(jtag,p)
      yj=yrap(jtag,p)
      yg=yrap(3,p)
      ayj=abs(yj)

!---- invariant mass of photon + jet 
      mgj=(p(3,4)+p(jtag,4))**2
      do i=1,3
         mgj=mgj-(p(3,i)+p(jtag,i))**2
      enddo
      mgj=sqrt(mgj)
      cos_star=tanh((yg-yj)/2._dp)

   
      r2= (p(3,1)*p(jtag,1)+p(3,2)*p(jtag,2))
     &     /(ptj*pt3)
      if (r2 > +0.9999999_dp) r2=+1._dp
      if (r2 < -0.9999999_dp) r2=-1._dp
      
      phi_gj=acos(r2)

  
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


!====== pt photon 
      
      call bookplot(n,tag,'pt_phot 1',pt3,wt,wt2,45._dp,55._dp
     &     ,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt_phot 2',pt3,wt,wt2,55._dp,100._dp
     &     ,15._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt_phot 3',pt3,wt,wt2,100._dp,150._dp
     &     ,25._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt_phot 4',pt3,wt,wt2,150._dp,200._dp
     &     ,50._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt_phot 5',pt3,wt,wt2,200._dp,400._dp
     &     ,200._dp,'log')
      n=n+1

!----- pt jet
      call bookplot(n,tag,'pt_jet 1',ptj,wt,wt2,40._dp,55._dp
     &     ,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt_jet 2',ptj,wt,wt2,55._dp,100._dp
     &     ,15._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt_jet 3',ptj,wt,wt2,100._dp,150._dp
     &     ,25._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt_jet 4',ptj,wt,wt2,150._dp,200._dp
     &     ,50._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt_jet 5',ptj,wt,wt2,200._dp,400._dp
     &     ,200._dp,'log')
      n=n+1

!----- y jet
      call bookplot(n,tag,'|y_jet|',ayj,wt,wt2,0._dp,2.37_dp
     &     ,.237_dp,'lin')
      n=n+1

!----- y jet
      call bookplot(n,tag,'|y_phot|',y3,wt,wt2,-3._dp,3._dp
     &     ,.2_dp,'lin')
      n=n+1
      
!----- phi gamj jet
      call bookplot(n,tag,'delta phi_gam j',phi_gj,wt,wt2,0._dp,3.2_dp
     &     ,.32_dp,'lin')
      n=n+1

!----- mgj
      call bookplot(n,tag,'m gam j 1',mgj,wt,wt2,161._dp,200._dp
     &     ,39._dp,'lin')
      n=n+1
      call bookplot(n,tag,'m gam j 2',mgj,wt,wt2,200._dp,400._dp
     &     ,100._dp,'lin')
      n=n+1
      call bookplot(n,tag,'m gam j 3',mgj,wt,wt2,400._dp,600._dp
     &     ,200._dp,'lin')
      n=n+1
      call bookplot(n,tag,'m gam j 4',mgj,wt,wt2,600._dp,1000._dp
     &     ,400._dp,'lin')
      n=n+1
      
!-----cos*
      call bookplot(n,tag,'cos gam j 1',abs(cos_star),wt,wt2,0._dp,0.7_dp
     &     ,0.1_dp,'log')
      n=n+1
      call bookplot(n,tag,'cos gam j 2',abs(cos_star),wt,wt2,0.7_dp,.83_dp
     &     ,0.13_dp,'log')
      n=n+1
      

      
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
