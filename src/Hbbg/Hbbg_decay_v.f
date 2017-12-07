!--- CW March 17, Routine to compute the H=> bb g virtual matrix element squared
      subroutine Hbbg_decay_v(p,i1,i2,i3,msqv)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'hbbparams.f'
      include 'epinv.f'
      include 'scale.f'
      include 'scheme.f'
      include 'nflav.f'
      real(dp) :: msqv,p(mxpart,4)
      integer i1,i2,i3
      complex(dp) :: vamp_lc(2,2),vamp_slc(2,2),loamps(2,2)
      integer h1,h2,imax
      real(dp) :: fac,faclo
      real(dp):: msqlo,msquv,hdecaylo
      logical nfonly
      common/nfonly/nfonly
!$omp threadprivate(/nfonly/)

      nfonly=.false.
!==== DEBUG
!==== KC testing code
!      epinv=0._dp
!      musq=1.1_dp**2
!      include 'phasespace_KC.f'
!      call spinorz(3,p,za,zb)
!      call Hbbg_vamp_lc(1,2,3,za,zb,vamp_lc)
!      call Hbbg_vamp_slc(1,2,3,za,zb,vamp_slc)
!      write(6,*) '*********** LC AMPLITUDES **************'
!      do h1=1,2
!         do h2=1,2
!            write(6,*) h1,h2,vamp_lc(h1,h2)*im
!         enddo
!      enddo
!      write(6,*) '*********** SLC AMPLITUDES **************'
!      do h1=1,2
!         do h2=1,2
!            write(6,*) h1,h2,vamp_slc(h1,h2)*im
!         enddo
!      enddo
!      write(6,*) '*****************************************'
!      stop
!==== END DEBUG
      


      scheme='dred'
      msqv=zip
      msqlo=zip
     
!---- LO factor 
      faclo=V*mb_eff**2*gsq*gwsq/wmass**2/four
!==== NLO correction
      fac = faclo*ason2pi
      
      imax=max(i1,i2,i3)
      
      call spinoru(imax,p,za,zb)
      call Hbbg_decay(p,i1,i2,i3,hdecaylo)
     
      call Hbbg_loamp(i1,i2,i3,za,zb,loamps)
      call Hbbg_vamp_lc(i1,i2,i3,za,zb,vamp_lc)
      call Hbbg_vamp_slc(i1,i2,i3,za,zb,vamp_slc)

     
      do h1=1,2
         do h2=1,2
            msqv=msqv-1._dp*real(vamp_lc(h1,h2)*conjg(loamps(h1,h2)),dp)*xn
     &           -1._dp/xn*real(vamp_slc(h1,h2)*conjg(loamps(h1,h2)),dp)
            msqlo=msqlo+real(loamps(h1,h2)*conjg(loamps(h1,h2)),dp)
         enddo
      enddo

      msqv=msqv*fac
      
      msquv=zip
      
!---- UV renormalization of alpha_s (note finite renorm to go to FDH) 
!----- b0 
      msquv=-epinv*(11._dp*xn/6._dp-(real(nflav,dp))/3._dp)+xn/6._dp
!-----Z_b
      msquv=msquv-3._dp*cf*epinv-cf
      
      msquv=msquv*hdecaylo*ason2pi
      msqv=msquv+msqv

!---- this bit for nf only
!      if(nfonly) then
!          msquv=epinv*((real(nflav,dp))/3._dp)
!          msqv=zip
!          msqv = msquv*hdecaylo*ason2pi
!      endif

!-----this bit for no nf
!      msqv=msqv-epinv*((real(nflav,dp))/3._dp)*hdecaylo*ason2pi
         
         
!---- DEBUG - Checking code
!----
!      write(6,*) '*********** LO AMPLITUDES **************'
!      do h1=1,2
!         do h2=1,2
!            write(6,*) h1,h2,loamps(h1,h2)
!         enddo
!      enddo
!      write(6,*) '*********** LC AMPLITUDES **************'
!      do h1=1,2
!         do h2=1,2
!            write(6,*) h1,h2,vamp_lc(h1,h2)*im
!         enddo
!      enddo
!      write(6,*) '*********** SLC AMPLITUDES **************'
!      do h1=1,2
!         do h2=1,2
!            write(6,*) h1,h2,vamp_slc(h1,h2)*im
!         enddo
!      enddo
!      write(6,*) '*****************************************'
!      stop
!----- END DEBUG 
      
      
      return
      end
