!---- CW Jan 17 Cutting routine for Z/ga ratio
      subroutine gencuts_ZgaRat(p,inchtcut,passed)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'jetcuts.f'
      include 'jetlabel.f'
      include 'kprocess.f'
      real(dp) :: p(mxpart,4)
      logical passed,inchtcut
      real(dp) :: ptVmin,mllmin,mllmax
      real(dp) :: yVmax,ptlmin,ylmax
      real(dp) :: ptV,yV,ht,htcut,mll
      real(dp) :: pt,pttwo,yrap,yraptwo,R
      integer j,jtag
      
      passed=.true.

      if (inchtcut) then
         htcut=300._dp
      else
         htcut=-1._dp 
      endif
!------cut variables
      ptVmin=100._dp
      yVmax=1.4_dp

!--- specific Z+j variables
      mllmin=71._dp
      mllmax=111._dp 
      ptlmin=20._dp
      ylmax=2.4_dp
      mll=zip

!------ first calculate desired quantitiues
      if((kcase==kZ_1jet).or.(kcase==kZ_2jet)) then
         ptV=pttwo(3,4,p)
         yV=yraptwo(3,4,p)
         mll=(p(3,4)+p(4,4))**2
     &      -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2
         mll=sqrt(mll)
         jtag=5
      elseif((kcase==kgamjet).or.(kcase==kgam_2j)) then
         ptV=pt(3,p)
         yV=yrap(3,p)
         jtag=4
      else
         write(6,*) 'Unrecongized process in Z/ga ratio routine',kcase
         stop
      endif

      if (ptV < ptVmin) then
         passed=.false.
         return
      endif

      if(abs(yV) > yVmax) then
         passed=.false.
         return
      endif
         
      
!--------define H_T as the scalar sum of transverse mom for jets with pt > 30
!-------and eta < 2.4 (jets are clustered elsewhere)
      ht=zip
      do j=jtag,jtag-1+jets
         ht=ht+pt(j,p)
      enddo

      if(ht < htcut) then
         passed=.false.
         return
      endif

!-----processes specific cuts
      if((kcase==kgamjet).or.(kcase==kgam_2j)) then
!---  photon cuts on R(3,jets)
         do j=4,3+jets
            if(R(3,j,p) < 0.5_dp) then
               passed=.false.
               return
            endif
         enddo
      elseif((kcase==kZ_1jet).or.(kcase==kZ_2jet)) then
! ---- cuts on
!     mll
!     ptl > 30 |eta(l)| < 2.4
!     R(ell,jets)
         if((mll > mllmax) .or. (mll < mllmin)) then
            passed=.false.
            return 
         endif
         if((pt(3,p) < ptlmin) .or. (pt(4,p) < ptlmin)) then 
            passed=.false.
            return 
         endif
         if((abs(yrap(3,p)) > ylmax) .or. (abs(yrap(4,p)) > ylmax)) then 
            passed=.false.
            return 
         endif         
         do j=5,4+jets
            if((R(3,j,p) < 0.5_dp) .or. (R(4,j,p)< 0.5_dp)) then
               passed=.false.
               return 
            endif
         enddo
      endif

      return
      end
