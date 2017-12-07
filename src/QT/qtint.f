      function qtint(r,wgt)
      implicit none
      include 'types.f'
      real(dp):: qtint
      include 'constants.f'
      include 'mxpart.f'
      include 'vegas_common.f'
      include 'sprods_com.f'
      include 'kprocess.f'
      include 'x1x2.f'
      include 'energy.f'
      include 'npart.f'
      include 'kpart.f'
      include 'dynamicscale.f'
      include 'initialscales.f'
      include 'first.f'
      include 'nf.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'facscale.f'
      include 'taucut.f'
      real(dp):: p(mxpart,4),pjet(mxpart,4),r(mxdim),W,xmsq,
     & val,val2,pswt,xjac,wgt,QB(2),z1,z2,flux,BrnRat
      real(dp):: msq(-nf:nf,-nf:nf),Kqq,Kgg,Kfac
      integer j,k,imode
      logical:: bin,includedipole
      integer, save :: iorder
      common/bin/bin
      common/BrnRat/BrnRat
!$omp threadprivate(iorder)

      qtint=0._dp

      W=sqrts**2

      call gen_lops(r,p,pswt,*999)
      call dotem(npart+2,p,s)

c----reject event if any s(i,j) is too small
c      call smalls(s,npart,*999)
c----reject event if any tau is too small
      call smalltau(p,npart,*999)

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out
      if (includedipole(0,p) .eqv. .false.) then
        goto 999
      endif
      
      if (dynamicscale) call scaleset(initscale,initfacscale,p)
     
      xx(1)=-2._dp*p(1,4)/sqrts
      xx(2)=-2._dp*p(2,4)/sqrts

      z1=r(ndim-1)**2
      z2=r(ndim)**2
      xjac=four*sqrt(z1*z2)

      flux=fbGeV2/(two*xx(1)*xx(2)*W)

      QB(1)=-two*p(1,4)
      QB(2)=-two*p(2,4)

c--- determine order of calculation on first call
      if (first) then
        first=.false.
        if (kpart==knnlo) then
          iorder=2
        elseif (kpart==ksnlo) then
          iorder=1
        else
          write(6,*) 'Error in qtint: kpart=',kpart
          stop
        endif
        if (coeffonly) iorder=-iorder
        call CalcParam(nf)
      endif

      call LoadAllPDF(xx(1),xx(2),z1,z2,facscale)

c--- Calculate the required matrix elements
      if     ((kcase==kW_only) .or. (kcase==kZ_only)
     &   .or. (kcase==kgamgam)) then
         imode=0
         if (kcase==kW_only) call qqb_w(p,msq)
         if (kcase==kZ_only) call qqb_z(p,msq)
         if (kcase==kgamgam) then
           call qqb_gamgam(p,msq)
           imode=1
         endif
         xmsq=0d0
         do j=-nf,nf
         do k=-nf,nf
         if (msq(j,k) .ne. 0d0) then
           Kfac=Kqq(iorder,j,k,xx(1),xx(2),z1,z2,
     &              as,scale,facscale,p(1:4,:),taucut**2,imode)
           xmsq=xmsq+msq(j,k)*Kfac
         endif
         enddo
         enddo
!         write(6,*) 'xmsq',xmsq
!         pause
!      elseif ((kcase==kggfus0) .or. (kcase==kHigaga)) then
!        call lumxmsq_h(p,xx,z1,z2,QB,iorder,xmsq)
!      elseif ((kcase==kWHbbar) .or. (kcase==kWHgaga)
!     &   .or. (kcase==kWH__WW)) then
!        call lumxmsq_wh(p,xx,z1,z2,QB,iorder,xmsq)
!      elseif ((kcase==kZHbbar) .or. (kcase==kZHgaga)
!     &   .or. (kcase==kZH__WW)) then
!        call lumxmsq_zh(p,xx,z1,z2,QB,iorder,xmsq)
!      elseif (kcase==kgamgam) then
!         call lumxmsq_gaga(p,xx,z1,z2,QB,iorder,xmsq)
      endif

      qtint=flux*xjac*pswt*xmsq/BrnRat

      call getptildejet(0,pjet)
      val=qtint*wgt 
      val2=val**2 
      
      if (val == val) then
        continue
      else
        write(6,*) 'Discarded NaN, val=',val
        write(6,*) 'flux,xjac,pswt,xmsq,BrnRat',flux,xjac,pswt,xmsq,BrnRat
!        write(6,*) 'wgt',wgt
!        write(6,*) '      r(1)=',r(1),'_dp'
!        write(6,*) '      r(2)=',r(2),'_dp'
!        write(6,*) '      r(3)=',r(3),'_dp'
!        write(6,*) '      r(4)=',r(4),'_dp'
!        write(6,*) '      r(5)=',r(5),'_dp'
!        write(6,*) '      r(6)=',r(6),'_dp'
        call flush(6)
        goto 999
      endif

      if (bin) then
         call nplotter(pjet,val,val2,0)
      endif

      return

 999  continue
      qtint=zip
      
      return
      end


