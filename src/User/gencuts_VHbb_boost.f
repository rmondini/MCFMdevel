!---- Roberto Mondini. Dec 2017

      function gencuts_VHbb_boost(p,pkt,psoft)

!===== cuts for V(l1,l2) H(b,bbar)
!===== label:    (3 , 4)
!===== boosted Higgs CMS analysis

!===== p is array of momenta (including incoming and leptons) before soft drop
!===== psoft is array of momenta (including incoming and leptons) after soft drop

      implicit none
      include 'types.f'
      include 'mpicommon.f'
      logical gencuts_VHbb_boost
      include 'constants.f'
      include 'mxpart.f'
      include 'kprocess.f'
      include 'nwz.f'
      include 'jetcuts.f'
      include 'clustering.f'
      include 'jetlabel.f' 
      include 'first.f'
      include 'energy.f'
      include 'softdrop_var.f'
      real(kind=dp):: p(mxpart,4),pkt(mxpart,4),psoft(mxpart,4)
      logical chooseptVcut,chooseptktjet
      real(dp) :: pt,yrap,onemassmod,n21
      real(dp) :: ptlmin,ylmax,metmin
      real(dp) :: ptbbjetmin,ybbjetmax,rhojetmin,rhojetmax
      real(dp) :: e2corr, e23corr
      integer i
      integer lepid,metid,j,njets
      integer :: jets_orig
      real(dp) :: lptjet,mjet,mjetbsd,rhojet,n21jet,ptV,ptjet(mxpart),e21val
      common/observables_VHbb_boost/lptjet,mjet,mjetbsd,rhojet,n21jet,ptV,ptjet,e21val
!$omp threadprivate(/observables_VHbb_boost/) 

      jets_orig=jets
      noj=jets_orig
      nj=jets

      ptjet(:)=0

      gencuts_VHbb_boost=.false.

!==== .true. means using ptV for selection pt>450GeV
!==== .false. means using pt of b-tagged leading jet for selection pt>450GeV
      chooseptVcut=.true.

!==== .true. means using pt of anti-kt jet for rhojet
!==== .false. means using pt of sd jet for rhojet
      chooseptktjet=.true.

!====== debug

!      write(*,*) 'phase-space point gencuts'

!      write(*,*)
!      do i=1,8
!         write(*,*) p(i,1),p(i,2),p(i,3),p(i,4)
!      enddo

!      write(*,*)
!      do i=1,8
!         write(*,*) psoft(i,1),psoft(i,2),psoft(i,3),psoft(i,4)
!      enddo

!      write(*,*)
!      write(*,*) 'npjr'
!      do i=1,4
!         write(*,*) npjr(i,1),npjr(i,2),npjr(i,3),npjr(i,4)
!      enddo

!      write(*,*)
!      write(*,*) 'n of SD jets = ',jets

!      write(*,*)
!      write(*,*) '*******************'

!      stop

!====== debug


!======= lepton cuts
      ptlmin=15._dp
      ylmax=2.5_dp
      metmin=15._dp

!======= CMS analysis cuts
      ptbbjetmin=450._dp
      ybbjetmax=2.5_dp
      rhojetmin=-6.0_dp
      rhojetmax=-2.1_dp

      if(first) then
         first=.false.
         if (rank == 0) then
            write(6,*)  '****************** V(l1,l2)H(b,bbar) cuts ********************'
            write(6,99) '*        pt(lepton)      >   ',ptlmin,
     &                ' GeV            *'
            write(6,99) '*        |eta(lepton)|      <   ',ylmax,
     &     '            *'
            write(6,99) '*        pt(V/jet)      >   ',ptbbjetmin,
     &                ' GeV            *'
            write(6,99) '*        |eta(jet)|      <   ',ybbjetmax,
     &     '            *'
            write(6,99) '*        rho(jet)      >   ',rhojetmin,
     &     '            *'
            write(6,99) '*        rho(jet)      <   ',rhojetmax,
     &     '            *'
            write(6,*)  '**************************************************'
         endif
      endif








!=====lepton cuts W 
!      if((kcase.eq.kWHbbar).or.(kcase.eq.kWH1jet)
!     &     .or.(kcase.eq.kWHbbdk)) then
!         if(nwz==1) then
!=========lepton is p4
!            lepid=4
!            metid=3
!         elseif(nwz==-1) then 
!=========lepton is p3
!            lepid=4
!            metid=3
!         else
!            write(6,*) 'unrecognized nwz =',nwz
!            stop
!         endif
                 
!=========lepton is p4
!         if((pt(lepid,p).lt.ptlmin).or.(abs(yrap(lepid,p)).gt.ylmax))then
!            gencuts_VHbb_boost=.true.
!            goto 999
!         endif
!=========met is p3
!         if(pt(metid,p).lt.metmin) then
!            gencuts_VHbb_boost=.true.
!            goto 999
!         endif





!=====lepton cuts Z
      if((kcase.eq.kZHbbnn).or.(kcase.eq.kZHbbjt).or.
     &  (kcase.eq.kZHn2ac)) then
!======check both leptons
         do lepid=3,4
            if((pt(lepid,p).lt.ptlmin).or.(abs(yrap(lepid,p)).gt.ylmax))then
               gencuts_VHbb_boost=.true.
               goto 999
            endif
         enddo
      endif

!=====find id of b-tagged largest-pt jet ("leading")

      idlptjet=0
      lptjet=0._dp

!==== safety check, make sure there is a jet to play with
      if(jets==0) then 
         gencuts_VHbb_boost=.true.
         goto 999
      endif
         
      do j=1,jets

         if( any(npjr(j,:)==5) .or. any(npjr(j,:)==6) ) then 
            ptjet(j)=pt(j+4,pkt)
            if(ptjet(j).gt.lptjet) then
               idlptjet=j
               if(chooseptktjet) then
                  lptjet=ptjet(j)
               else
                  lptjet=pt(idlptjet+4,psoft)
               endif
            endif
         endif

      enddo

!-----check that we have a b-tagged jet
      if(idlptjet==0) then
         gencuts_VHbb_boost=.true.
         goto 999
      endif
         
      
      nplj = count(pinSDj(idlptjet,:)/=0)

 
!=====mass of jet before soft drop

      mjetbsd=onemassmod(idlptjet+4,pkt)

!=====pt cut: pt>450 GeV

      ptV=zip
      do j=1,2
         ptV=ptV+(p(3,j)+p(4,j))**2
      enddo
      ptV=sqrt(ptV)

      if(chooseptVcut) then
         if(ptV.lt.ptbbjetmin) then
            gencuts_VHbb_boost=.true.
            goto 999
         endif
      else
         if(lptjet.lt.ptbbjetmin) then
            gencuts_VHbb_boost=.true.
            goto 999
         endif
      endif

!=====|yrap leading jet|<2.5

      if(abs(yrap(idlptjet+4,pkt)).gt.ybbjetmax) then
         gencuts_VHbb_boost=.true.
         goto 999
      endif

!=====rhojet cut: -6.0<rhojet<-2.1
!=====rhojet=log(mjet**2/ptjet**2)

      mjet=onemassmod(idlptjet+4,psoft)

  
      
      rhojet=log(mjet**2/lptjet**2)

      if((rhojet.lt.rhojetmin).or.(rhojet.gt.rhojetmax)) then
         gencuts_VHbb_boost=.true.
         goto 999
      endif

 !     if(count(pinSDj(idlptjet,:)/=0)/=count(npjr(idlptjet,:)/=0)) then
 !        write(6,*) pinSDj(idlptjet,:)
 !        write(6,*) npjr(idlptjet,:)
 !        
 !        pause
 !     endif
      
!=====n21ddt cut: n21ddt=n21-n21(26%) <= 0

      n21jet=n21(p,psoft)
      e21val=e2corr(p)
!      write(6,*) n21(p,psoft)/n21_old(p,psoft)
      
!      if(n21jet .gt. 0._dp) then
!         gencuts_VHbb_boost=.true.
!         goto 999
!      endif










!======= debug

!      if(nplj.eq.3) then

!      if(jets.eq.3) then

!      write(*,*)
!      write(*,*) '*******************'
!      write(*,*) 'phase-space point after all cuts'

!      write(*,*)
!      do i=1,8
!         write(*,*) p(i,1),p(i,2),p(i,3),p(i,4)
!      enddo
!      write(*,*)
!      do i=1,8
!         write(*,*) psoft(i,1),psoft(i,2),psoft(i,3),psoft(i,4)
!      enddo
!      write(*,*)
!      write(*,*) 'npjr'
!      do i=1,4
!         write(*,*) npjr(i,1),npjr(i,2),npjr(i,3),npjr(i,4)
!      enddo
!      write(*,*)
!      write(*,*) 'n of SD jets = ',jets
!      write(*,*)
!      write(*,*) 'id leading-pt jet = ',idlptjet
!      write(*,*) 'n of particles in leading jet = ',nplj
!      write(*,*) 'leading pt kt jet (GeV) = ',pt(idlptjet+4,pkt)
!      write(*,*) 'leading pt sd jet (GeV) = ',pt(idlptjet+4,psoft)
!      write(*,*) 'mjet sd (GeV) = ',mjet
!      write(*,*) 'mjet before sd (GeV) = ',mjetbsd
!      write(*,*) 'yrap leading jet before sd = ',yrap(idlptjet+4,pkt)
!      write(*,*) 'rhojet = ',rhojet
!      write(*,*) 'n21jet = ',n21jet
!      write(*,*) '*******************'
!      write(*,*)

!      endif

!      stop

!======= debug











!=====reset jets to not disrupt any future routines
      jets=jets_orig       
      return
 99   format(1x,a29,f6.2,a17)
 999  continue
      jets=jets_orig 
      return 
      end
      
