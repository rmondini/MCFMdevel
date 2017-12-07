      function gencuts_ATLAS_gamjet(p)
!---- cutting routine for ATLAS gamjet analysis 1307.6795 
!
!---- it is assumed that basic jet and photon cuts are done in the input
!---- file, and that this routine will just do the cracks and advanced cuts
      implicit none
      include 'types.f'
      logical gencuts_ATLAS_gamjet
      include 'constants.f'
      include 'mxpart.f'
      include 'jetlabel.f'
      real(dp) :: p(mxpart,4)
      real(dp) :: cos_star_max,yrap,etarap,pti
      real(dp) :: mgj_min,cos_star, ptmax,yj,mgj,etag,pt
      integer jtag,i

      gencuts_ATLAS_gamjet=.false.

      mgj_min=161._dp
      cos_star_max = 0.83_dp

! rapidity of photon
      etag=etarap(3,p)
!---- remove ATLAS crack
      if( (abs(etag) > 1.37_dp) .and. (abs(etag) < 1.52_dp)) then
        gencuts_ATLAS_gamjet=.true.
        return
      endif

!---- find hardest jet for additional cuts
      jtag=-1
      ptmax=zip
      do i=4,3+jets
        pti=pt(i,p)
        if(pti > ptmax) then
          jtag=i
          ptmax=pti
        endif
      enddo

      yj=yrap(jtag,p)
!-----ensure hardest jet is in fiducial vol
      if(abs(yj) > 2.37_dp) then
        gencuts_ATLAS_gamjet=.true.
        return
      endif

!      return

!-----cut on |eta_gam + y _jet|
      if(abs(etag+yj) > 2.37_dp) then
         gencuts_ATLAS_gamjet=.true.
         return
      endif     
     
!---- invariant mass of photon + jet 
      mgj=(p(3,4)+p(jtag,4))**2
      do i=1,3
         mgj=mgj-(p(3,i)+p(jtag,i))**2
      enddo
      mgj=sqrt(max(mgj,zip))
      if(mgj < mgj_min) then
        gencuts_ATLAS_gamjet=.true.
        return
      endif

!---- cos theta*
      cos_star=tanh((etag-yj)/2._dp)
      if(abs(cos_star) > cos_star_max) then
        gencuts_ATLAS_gamjet=.true.
        return
      endif
      
      return
      end
