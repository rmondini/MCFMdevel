      subroutine nplotter_dirgam(p,wt,wt2,switch,nd)
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
      include 'cplx.h'
      include 'histo.f'
      include 'jetlabel.f'
      include 'frag.f'
      include 'taucut.f' 
      include 'phot_dip.f'
      include 'outputflags.f'
      
      real(dp):: p(mxpart,4),wt,wt2
      real(dp):: yrap,pt,r,mjj
      real(dp):: y3,y4,y5,pt3,pt4,pt5,r45,r35,r34,s34,m34
      integer:: switch,n,nplotmax
      integer tag
      integer:: nd
      logical, save::first=.true.
      common/nplotmax/nplotmax
ccccc !$omp threadprivate(first,/nplotmax/)
      real(dp) :: ptreg1,ptreg2,ptreg3,ptreg4
      real(dp) tauhist
      integer mjets
      common/tauhist/tauhist
!$omp threadprivate(/tauhist/)

      
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
        tauhist=-1d8
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
        jets=-1
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
    
      if(jets<=1) then 
         pt3=pt(3,p)
         pt4=pt(4,p)
         y3=yrap(3,p)
         y4=yrap(4,p)
         r34=R(p,3,4)
         s34=2._dp*(p(4,4)*p(3,4)-p(4,1)*p(3,1)-p(4,2)*p(3,2)
     &           -p(4,3)*p(3,3))
         pt5=1d7
         y5=1d7
         r45=1d7
         r35=1d7
      else
         pt3=pt(3,p)
         y3=yrap(3,p)
         if(pt(4,p)>pt(5,p)) then 
            pt4=pt(4,p)
            y4=yrap(4,p)
            pt5=pt(5,p) 
            y5=yrap(5,p) 
            r34=R(p,3,4)
            r35=R(p,3,5)
            s34=2._dp*(p(4,4)*p(3,4)-p(4,1)*p(3,1)-p(4,2)*p(3,2)
     &              -p(4,3)*p(3,3))
         else
            pt4=pt(5,p) 
            y4=yrap(5,p) 
            pt5=pt(4,p) 
            y5=yrap(4,p) 
            r34=R(p,3,5)
            r35=R(p,3,4)
            s34=2._dp*(p(5,4)*p(3,4)-p(5,1)*p(3,1)-p(5,2)*p(3,2)
     &              -p(5,3)*p(3,3))
         endif
         r45=R(p,4,5)
         mjj=2._dp*(p(4,4)*p(5,4)-p(4,1)*p(5,1)-p(4,2)*p(5,2)
     &           -p(4,3)*p(5,3))
  
         mjj=sqrt(mjj)
      endif
                  
      m34=sqrt(s34)

      ptreg1=-1d8
      ptreg2=-1d8
      ptreg3=-1d8
      ptreg4=-1d8
!---- ATLAS 8 TeV
      if(abs(y3) < 0.6_dp) then
         ptreg1=pt3
      elseif((abs(y3) >= 0.6_dp).and.(abs(y3) < 1.37_dp)) then 
         ptreg2=pt3
      elseif((abs(y3) >= 1.56_dp).and.(abs(y3) < 1.81_dp)) then 
         ptreg3=pt3
      elseif((abs(y3) >= 1.81_dp).and.(abs(y3) < 2.37_dp)) then 
         ptreg4=pt3
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

      if(usescet) then 
!----- tau 1 for studies, low tau region
      call bookplot(n,tag,'tauhist low',tauhist,wt,wt2,taucut,
     &   taucut*100._dp, taucut*5._dp ,'log') 
      n=n+1
!----- tau 1 overflow
      call bookplot(n,tag,'tauhist over',tauhist,wt,wt2,
     &   taucut*100._dp,2000._dp, 50._dp ,'log') 
      n=n+1
      endif

!      write(6,*) real(jets)
!      pause
!---  njets
      call bookplot(n,tag,'jets',real(jets,dp),wt,wt2,-0.5_dp,4.5_dp,1._dp
     &     ,'lin') 
      n=n+1
      
!-----each region requires three histograms 1 with bins of 5 GeV covers
!---- 25-525 GeV (with the last complete bin being the 400-470 Bin
!---- the second does 400-1400 in 10 GeV Bins
!---- for region 1 there is an overflow bin of 1100-1500 GeV 

      
      call bookplot(n,tag,'pt_phot_inc',pt3,wt,wt2,25._dp,2500._dp
     &     ,25._dp,'log')
      n=n+1      

      call bookplot(n,tag,'pt_phot_y0_inc',ptreg1,wt,wt2,0._dp,2500._dp
     &     ,25._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'pt_phot_y1_inc',ptreg2,wt,wt2,0._dp,2500._dp
     &     ,25._dp,'log')
      n=n+1

      call bookplot(n,tag,'pt_phot_y2_inc',ptreg3,wt,wt2,0._dp,2500._dp
     &     ,25._dp,'log')
      n=n+1

      call bookplot(n,tag,'pt_phot_y3_inc',ptreg4,wt,wt2,0._dp,2500._dp
     &     ,25._dp,'log')
      n=n+1

!---- region 1       
      call bookplot(n,tag,'pt_phot_y0_1',ptreg1,wt,wt2,25._dp,525._dp
     &     ,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt_phot_y0_2',ptreg1,wt,wt2,400._dp,1400._dp
     &     ,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt_phot_y0_3',ptreg1,wt,wt2,1100._dp,1500._dp
     &     ,400._dp,'log')
      n=n+1

      

!---- region 2
      call bookplot(n,tag,'pt_phot_y1_1',ptreg2,wt,wt2,25._dp,525._dp
     &     ,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt_phot_y1_2',ptreg2,wt,wt2,400._dp,1400._dp
     &     ,10._dp,'log')
      n=n+1


!---- region 3
      call bookplot(n,tag,'pt_phot_y2_1',ptreg3,wt,wt2,25._dp,525._dp
     &     ,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt_phot_y2_2',ptreg3,wt,wt2,400._dp,1400._dp
     &     ,10._dp,'log')
      n=n+1


      call bookplot(n,tag,'pt_phot_y3_1',ptreg4,wt,wt2,25._dp,525._dp
     &     ,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'pt_phot_y3_2',ptreg4,wt,wt2,400._dp,1400._dp
     &     ,10._dp,'log')
      n=n+1

      
c$$$      call bookplot(n,tag,'pt_phot_inc',pt3,wt,wt2,25._dp,2500._dp
c$$$     &     ,25._dp,'log')
c$$$      n=n+1
c$$$      
c$$$
c$$$      call bookplot(n,tag,'pt_phot_y0_inc',ptreg1,wt,wt2,0._dp,2500._dp
c$$$     &     ,25._dp,'log')
c$$$      n=n+1
c$$$      
c$$$!====== 10 GeV Bins for ATLAS first rapidity region
c$$$
c$$$      
c$$$      call bookplot(n,tag,'pt_phot_y0_1',ptreg1,wt,wt2,25._dp,85._dp
c$$$     &     ,10._dp,'log')
c$$$      n=n+1
c$$$!===== 20 GeV Bins, ditto 
c$$$      call bookplot(n,tag,'pt_phot_y0_2',ptreg1,wt,wt2,85._dp,125._dp
c$$$     &     ,20._dp,'log')
c$$$      n=n+1
c$$$!===== 25 GeV Bins, ditto 
c$$$      call bookplot(n,tag,'pt_phot_y0_3',ptreg1,wt,wt2,125._dp,200._dp
c$$$     &     ,25._dp,'log')
c$$$      n=n+1
c$$$!===== 50 GeV Bins, ditto 
c$$$      call bookplot(n,tag,'pt_phot_y0_4',ptreg1,wt,wt2,200._dp,400._dp
c$$$     &     ,50._dp,'log')
c$$$      n=n+1
c$$$!===== 70 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y0_5',ptreg1,wt,wt2,400._dp,470._dp
c$$$     &     ,70._dp,'log')
c$$$      n=n+1
c$$$!===== 80 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y0_6',ptreg1,wt,wt2,470._dp,550._dp
c$$$     &     ,80._dp,'log')
c$$$      n=n+1
c$$$!=====100 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y0_7',ptreg1,wt,wt2,550._dp,750._dp
c$$$     &     ,100._dp,'log')
c$$$      n=n+1
c$$$!=====150 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y0_8',ptreg1,wt,wt2,750._dp,900._dp
c$$$     &     ,150._dp,'log')
c$$$      n=n+1
c$$$!=====200 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y0_9',ptreg1,wt,wt2,900._dp,1100._dp
c$$$     &     ,200._dp,'log')
c$$$      n=n+1
c$$$!=====400 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y0_10',ptreg1,wt,wt2,1100._dp,1500._dp
c$$$     &     ,400._dp,'log')
c$$$      n=n+1
c$$$
c$$$
c$$$      call bookplot(n,tag,'pt_phot_y1_inc',ptreg2,wt,wt2,0._dp,2000._dp
c$$$     &     ,50._dp,'log')
c$$$      n=n+1
c$$$!====== 10 GeV Bins for ATLAS second rapidity region
c$$$      
c$$$      call bookplot(n,tag,'pt_phot_y1_1',ptreg2,wt,wt2,25._dp,85._dp
c$$$     &     ,10._dp,'log')
c$$$      n=n+1
c$$$!===== 20 GeV Bins, ditto 
c$$$      call bookplot(n,tag,'pt_phot_y1_2',ptreg2,wt,wt2,85._dp,125._dp
c$$$     &     ,20._dp,'log')
c$$$      n=n+1
c$$$!===== 25 GeV Bins, ditto 
c$$$      call bookplot(n,tag,'pt_phot_y1_3',ptreg2,wt,wt2,125._dp,200._dp
c$$$     &     ,25._dp,'log')
c$$$      n=n+1
c$$$!===== 50 GeV Bins, ditto 
c$$$      call bookplot(n,tag,'pt_phot_y1_4',ptreg2,wt,wt2,200._dp,400._dp
c$$$     &     ,50._dp,'log')
c$$$      n=n+1
c$$$!===== 70 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y1_5',ptreg2,wt,wt2,400._dp,470._dp
c$$$     &     ,70._dp,'log')
c$$$      n=n+1
c$$$!===== 80 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y1_6',ptreg2,wt,wt2,470._dp,550._dp
c$$$     &     ,80._dp,'log')
c$$$      n=n+1
c$$$!=====100 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y1_7',ptreg2,wt,wt2,550._dp,750._dp
c$$$     &     ,100._dp,'log')
c$$$      n=n+1
c$$$!=====150 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y1_8',ptreg2,wt,wt2,750._dp,900._dp
c$$$     &     ,150._dp,'log')
c$$$      n=n+1
c$$$!=====200 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y1_9',ptreg2,wt,wt2,900._dp,1100._dp
c$$$     &     ,200._dp,'log')
c$$$      n=n+1
c$$$!=====400 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y1_10',ptreg2,wt,wt2,1100._dp,1500._dp
c$$$     &     ,400._dp,'log')
c$$$      n=n+1
c$$$
c$$$
c$$$            
c$$$      call bookplot(n,tag,'pt_phot_y2_inc',ptreg3,wt,wt2,0._dp,2000._dp
c$$$     &     ,50._dp,'log')
c$$$      n=n+1
c$$$!====== 10 GeV Bins for ATLAS third rapidity region
c$$$      
c$$$      call bookplot(n,tag,'pt_phot_y2_1',ptreg3,wt,wt2,25._dp,85._dp
c$$$     &     ,10._dp,'log')
c$$$      n=n+1
c$$$!===== 20 GeV Bins, ditto 
c$$$      call bookplot(n,tag,'pt_phot_y2_2',ptreg3,wt,wt2,85._dp,125._dp
c$$$     &     ,20._dp,'log')
c$$$      n=n+1
c$$$!===== 25 GeV Bins, ditto 
c$$$      call bookplot(n,tag,'pt_phot_y2_3',ptreg3,wt,wt2,125._dp,200._dp
c$$$     &     ,25._dp,'log')
c$$$      n=n+1
c$$$!===== 50 GeV Bins, ditto 
c$$$      call bookplot(n,tag,'pt_phot_y2_4',ptreg3,wt,wt2,200._dp,400._dp
c$$$     &     ,50._dp,'log')
c$$$      n=n+1
c$$$!===== 70 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y2_5',ptreg3,wt,wt2,400._dp,470._dp
c$$$     &     ,70._dp,'log')
c$$$      n=n+1
c$$$!===== 80 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y2_6',ptreg3,wt,wt2,470._dp,550._dp
c$$$     &     ,80._dp,'log')
c$$$      n=n+1
c$$$!=====100 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y2_7',ptreg3,wt,wt2,550._dp,750._dp
c$$$     &     ,100._dp,'log')
c$$$      n=n+1
c$$$!=====150 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y2_8',ptreg3,wt,wt2,750._dp,900._dp
c$$$     &     ,150._dp,'log')
c$$$      n=n+1
c$$$!=====200 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y2_9',ptreg3,wt,wt2,900._dp,1100._dp
c$$$     &     ,200._dp,'log')
c$$$      n=n+1
c$$$!=====400 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y2_10',ptreg3,wt,wt2,1100._dp,1500._dp
c$$$     &     ,400._dp,'log')
c$$$      n=n+1
c$$$
c$$$         call bookplot(n,tag,'pt_phot_y3_inc',ptreg4,wt,wt2,0._dp,2000._dp
c$$$     &     ,50._dp,'log')
c$$$      n=n+1
c$$$
c$$$!====== 10 GeV Bins for ATLAS final rapidity region
c$$$      
c$$$      call bookplot(n,tag,'pt_phot_y3_1',ptreg4,wt,wt2,25._dp,85._dp
c$$$     &     ,10._dp,'log')
c$$$      n=n+1
c$$$!===== 20 GeV Bins, ditto 
c$$$      call bookplot(n,tag,'pt_phot_y3_2',ptreg4,wt,wt2,85._dp,125._dp
c$$$     &     ,20._dp,'log')
c$$$      n=n+1
c$$$!===== 25 GeV Bins, ditto 
c$$$      call bookplot(n,tag,'pt_phot_y3_3',ptreg4,wt,wt2,125._dp,200._dp
c$$$     &     ,25._dp,'log')
c$$$      n=n+1
c$$$!===== 50 GeV Bins, ditto 
c$$$      call bookplot(n,tag,'pt_phot_y3_4',ptreg4,wt,wt2,200._dp,400._dp
c$$$     &     ,50._dp,'log')
c$$$      n=n+1
c$$$!===== 70 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y3_5',ptreg4,wt,wt2,400._dp,470._dp
c$$$     &     ,70._dp,'log')
c$$$      n=n+1
c$$$!===== 80 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y3_6',ptreg4,wt,wt2,470._dp,550._dp
c$$$     &     ,80._dp,'log')
c$$$      n=n+1
c$$$!=====100 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y3_7',ptreg4,wt,wt2,550._dp,750._dp
c$$$     &     ,100._dp,'log')
c$$$      n=n+1
c$$$!=====150 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y3_8',ptreg4,wt,wt2,750._dp,900._dp
c$$$     &     ,150._dp,'log')
c$$$      n=n+1
c$$$!=====200 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y3_9',ptreg3,wt,wt2,900._dp,1100._dp
c$$$     &     ,200._dp,'log')
c$$$      n=n+1
c$$$!=====400 GeV Bins, ditto
c$$$      call bookplot(n,tag,'pt_phot_y3_10',ptreg3,wt,wt2,1100._dp,1500._dp
c$$$     &     ,400._dp,'log')
c$$$      n=n+1
c$$$      
c$$$      call bookplot(n,tag,'pt_phot',pt3,wt,wt2,0._dp,200._dp,5._dp,'log')
c$$$      n=n+1
c$$$      call bookplot(n,tag,'y_phot',y3,wt,wt2,-3._dp,3._dp,0.1_dp,'lin')
c$$$      n=n+1
c$$$      call bookplot(n,tag,'pt_j1',pt4,wt,wt2,0._dp,200._dp,5._dp,'log')
c$$$      n=n+1 
c$$$      call bookplot(n,tag,'y_j1',y4,wt,wt2,-5._dp,5._dp,0.2_dp,'lin')
c$$$      n=n+1  
c$$$      call bookplot(n,tag,'DeltaR34',r34,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
c$$$      n=n+1
c$$$     
c$$$!    m34      
c$$$      call bookplot(n,tag,'m_gamj',m34,wt,wt2,0._dp,300._dp,20._dp,'log')
c$$$      n=n+1
c$$$      if(jets>1) then 
c$$$      call bookplot(n,tag,'y_j2',y5,wt,wt2,-3._dp,3._dp,0.2_dp,'lin')
c$$$      n=n+1
c$$$      call bookplot(n,tag,'pt_j2',pt5,wt,wt2,0._dp,200._dp,2._dp,'log')
c$$$      n=n+1
c$$$      call bookplot(n,tag,'mjj',mjj,wt,wt2,0._dp,150._dp,5._dp,'log')
c$$$      n=n+1
c$$$      call bookplot(n,tag,'DeltaR35',r35,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
c$$$      n=n+1
c$$$      call bookplot(n,tag,'DeltaR45',r45,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
c$$$      n=n+1
c$$$      endif
  
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
      
      

