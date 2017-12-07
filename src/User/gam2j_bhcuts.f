!===== CW July 16
!===== routine to 1106.1423

!----- jets are initially clustered elsewhere with a definition
!----- pt_j > 30 and |eta| < 5

!===== this routine calculates the MET and H_T as defined in 1106.1423
!===== and imposes cuts on these variables

!===== note that photon cuts are defined as usual (pt > 100, |eta| < 2.5)
      subroutine gam2j_bhcuts(p,passed)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'jetcuts.f'
      include 'jetlabel.f'
      include 'kprocess.f'
      real(dp):: p(mxpart,4)
      logical passed
      real(dp):: ht,met,ptjdef
      integer set,i,high_pt_j,ntagjet
      real(dp) :: htmin,metmin,metvec(4),phi_jmet,phimin
      real(dp) :: ptj,r2,r1,ejet,emet,dely,rmetj,rmin,ptj1,ptj2,ptj3
      real(dp) :: etarap,yrap,pt
      real(dp) :: tagjet(3,4)
      real(dp) :: pttwo
      integer jtag,j1,j2,j3
      
!==== set chooses which region we are in (c.f. 1106.1423, Section IV)
      set=3

      passed=.true.
      
      if(jets<2) then
! should not get here if jet cuts are working properly
         passed =.false.
         return
      endif
      
      if(set==1) then
         htmin=300._dp
         metmin=250._dp
      elseif(set==2) then
         htmin=500._dp
         metmin=150._dp
      elseif(set==3) then
         htmin=300._dp
         metmin=150._dp
      else
         write(6,*) 'unrecognized set in gam2j_bhcuts.f',set
         stop
      endif
      phimin=0.5_dp
      if ((kcase==kgam_2j) .or. (kcase==kgam_3j)) then 
         rmin=0.4_dp
         jtag=4
         j1=4
         j2=5
         j3=6
      elseif ((kcase==kZ_2jet) .or. (kcase==kZ_3jet)) then 
         rmin=0._dp
         jtag=5
         j1=5
         j2=6
         j3=7
!-----do pt check on 34 here
         if(pttwo(3,4,p).lt.100._dp) then
            passed=.false.
            return
         endif
      endif
      
!--------define H_T as the scalar sum of transverse mom for jets with pt > 50
!-------and eta < 2.5
!----- at the same time require the presence of at least two tagging jets

      ht=zip
      high_pt_j=0
      metvec(:)=zip
      tagjet(:,:)=zip

!-----note that we assume jets are properly defined elsewhere (>30 GeV eta < 5)
!---- such that met = - vector sum of all jets 
      do i = jtag,jtag-1+jets
         
         if((pt(i,p) > 50._dp) .and.
     &        (abs(etarap(i,p)) < 2.5_dp)) then
            ht=ht+pt(i,p)
            high_pt_j=high_pt_j+1
         endif
         metvec(:)=metvec(:)-p(i,:)
      enddo

      if ((kcase==kgam_2j) .or. (kcase==kZ_2jet)) then
        ntagjet=2
      else
        ntagjet=3
      endif
      
!-----check we have at least "ntagjet" tagging jets
      if(high_pt_j < ntagjet) then
         passed=.false.
         return
      endif

!-----check ht > htmin
      if(ht<htmin) then
         passed=.false.
         return
      endif

!----- now calculate met
      met=sqrt(metvec(1)**2+metvec(2)**2)
      if(met < metmin) then
         passed=.false.
         return
      endif

!==== calculate tagging jets
      if(jets==2) then
         tagjet(1,:)=p(j1,:)
         tagjet(2,:)=p(j2,:)
      elseif(jets==3) then
         ptj1=pt(j1,p)
         ptj2=pt(j2,p)
         ptj3=pt(j3,p)
         if((ptj1.gt.ptj2).and.(ptj1.gt.ptj3)) then
            if(ptj2.gt.ptj3) then
               tagjet(1,:)=p(j1,:)
               tagjet(2,:)=p(j2,:)
               tagjet(3,:)=p(j3,:)
            else
               tagjet(1,:)=p(j1,:)
               tagjet(2,:)=p(j3,:)
               tagjet(3,:)=p(j2,:)
            endif
         elseif((ptj2.gt.ptj1).and.(ptj2.gt.ptj3)) then
            if(ptj1.gt.ptj3) then
               tagjet(1,:)=p(j2,:)
               tagjet(2,:)=p(j1,:)
               tagjet(3,:)=p(j3,:)
            else
               tagjet(1,:)=p(j2,:)
               tagjet(2,:)=p(j3,:)
               tagjet(3,:)=p(j1,:)
            endif
         elseif((ptj3.gt.ptj1).and.(ptj3.gt.ptj2)) then
           if(ptj1.gt.ptj2) then
               tagjet(1,:)=p(j3,:)
               tagjet(2,:)=p(j1,:)
               tagjet(3,:)=p(j2,:)
            else
               tagjet(1,:)=p(j3,:)
               tagjet(2,:)=p(j2,:)
               tagjet(3,:)=p(j1,:)
            endif
         endif
      else
         write(6,*) 'Illegal # of jets in gam2j_bhcuts: ',jets
         stop       
      endif

!==== TAG JET CUTS
!======phi(MET,j) and R(MET,j) check
      do i=1,ntagjet
         ptj=sqrt(tagjet(i,1)**2+tagjet(i,2)**2)
         r2= (tagjet(i,1)*metvec(1)+tagjet(i,2)*metvec(2))
     &        /(ptj*met)
         if (r2 > +0.9999999_dp) r2=+1._dp
         if (r2 < -0.9999999_dp) r2=-1._dp
         
         phi_jmet=acos(r2)
         
         if (i == 3) then ! relaxed phi cut for third jet in three-jet case only
           phimin=phimin-0.2_dp
         endif
         if(phi_jmet<phimin) then
           passed=.false.
           return
         endif

c--- note that dR > dphi, so that if rmin < phimin then
c--- above cut on dphi automatically guarantees dR cut
         if (rmin > phimin) then
         
           ejet=sqrt(ptj**2+tagjet(i,3)**2)
           emet=sqrt(met**2+metvec(3)**2)
           
           r1= (ejet+tagjet(i,3))*(emet-metvec(3))/
     &          ((emet+metvec(3))*(ejet-tagjet(i,3)))
           dely=0.5_dp*log(r1)

           rmetj=sqrt(dely**2+phi_jmet**2)
           
           if(rmetj < rmin) then
!              write(6,*) rmetj,rmin
              passed=.false.
              return
           endif
         
         endif
         
      enddo

      return
      end
