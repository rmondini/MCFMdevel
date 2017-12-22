      subroutine nplotter_VHbbar_boost(p,wt,wt2,switch,nd)
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
      
      include 'vegas_common.f'
      include 'constants.f'
      include 'nf.f'
      include 'kprocess.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'histo.f'
      include 'outputflags.f'
      include 'interference.f'
      include 'jetlabel.f'
      include 'softdrop_var.f'
      integer nd
      real(dp):: p(mxpart,4),wt,wt2
      integer:: switch,n
      integer i,tag
      logical, save::first=.true.
      integer nplotmax
      real(dp) :: lptjet,mjet,mjetbsd,rhojet,n21jet,ptV,ptjet(mxpart)
      common/nplotmax/nplotmax
      common/observables_VHbb_boost/lptjet,mjet,mjetbsd,rhojet,n21jet,ptV,ptjet
!$omp threadprivate(/observables_VHbb_boost/)

************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************

      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag=tagbook
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

!=== they are common observables so they are already defined

!==== debug
!      if(nplj.eq.4) then
!      if(nj.eq.2) then

!      write(*,*)
!      write(*,*) '*******************'
!      write(*,*) 'phase-space point in nplotter'

!      write(*,*)
!      do i=1,8
!         write(*,*) p(i,1),p(i,2),p(i,3),p(i,4)
!      enddo
!      write(*,*) 'n of SD jets = ',nj
!      write(*,*)
!      write(*,*) 'id leading-pt jet = ',idlptjet
!      write(*,*) 'n of particles in leading jet = ',nplj
!      write(*,*) 'leading pt sd jet (GeV) = ',lptjet
!      write(*,*) 'mjet sd (GeV) = ',mjet
!      write(*,*) 'mjet before sd (GeV) = ',mjetbsd
!      write(*,*) 'rhojet = ',rhojet
!      write(*,*) 'n21jet = ',n21jet
!      write(*,*) '*******************'
!      write(*,*)

!      endif
!==== debug

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

!=====plots

      call bookplot(n,tag,'mSD leading jet (40-201 GeV)',
     &     mjet,wt,wt2,40._dp,201._dp,7._dp,'log')
      n=n+1

      call bookplot(n,tag,'mSD leading jet (0-500 GeV)',
     &     mjet,wt,wt2,0._dp,500._dp,5._dp,'log')
      n=n+1

      call bookplot(n,tag,'mass leading jet before soft drop (40-201 GeV)',
     &     mjetbsd,wt,wt2,40._dp,201._dp,7._dp,'log')
      n=n+1

      call bookplot(n,tag,'mass leading jet before soft drop (0-500 GeV)',
     &     mjetbsd,wt,wt2,0._dp,500._dp,5._dp,'log')
      n=n+1

      call bookplot(n,tag,'pt leading jet (0-1000 GeV)',
     &     lptjet,wt,wt2,0._dp,1000._dp,10._dp,'log')
      n=n+1

      call bookplot(n,tag,'ptV (400-1000 GeV)',
     &     ptV,wt,wt2,400._dp,1000._dp,10._dp,'log')
      n=n+1

      call bookplot(n,tag,'rho jet (-6,-2)',
     &     rhojet,wt,wt2,-6._dp,-2._dp,0.1_dp,'log')
      n=n+1


      call bookplot(n,tag,'n21 (0,0.06)',
     &     n21jet,wt,wt2,0._dp,0.06_dp,0.06_dp,'log')
      n=n+1
      call bookplot(n,tag,'n21 (0,0.06)',
     &     n21jet,wt,wt2,0.06_dp,0.4_dp,0.02_dp,'log')
      n=n+1
      
      call bookplot(n,tag,'n21 (0,0.5)',
     &     n21jet,wt,wt2,0._dp,0.5_dp,0.02_dp,'log')
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

