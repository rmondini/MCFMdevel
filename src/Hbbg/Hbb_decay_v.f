      subroutine Hbb_decay_v(p,i1,i2,msqv)
!----- integrated dipoles and virtual ME for H=>bb 
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
 !     complex(dp) :: vamp_lc(2,2),vamp_slc(2,2),loamps(2,2)
      integer h1,h2,imax
      real(dp) :: fac,s56
      real(dp) :: vpole,msq_z,msqlo,msquv
      real(dp)  :: intdip,xl12,dot,dipHD_qq
      complex(dp) :: lnrat

      scheme='dred'

      
      call spinoru(max(i1,i2),p,za,zb)
      s56=s(i1,i2)
!----- lo
      msqlo=xn*gwsq*mb_eff**2/(4._dp*wmass**2)*2._dp*s56
!---  virt FDH
      if(scheme.eq.'dred') then 
         vpole=-real((CF*two*(2*epinv**2
     &        + 2*epinv*lnrat(musq,-s56) + lnrat(musq,-s56)**2))/2._dp,dp)
!---- UV counter term FDH 
         msquv = -3._dp*cf*epinv-cf
      elseif(scheme.eq.'tH-V') then 
!-----virt CDR
      vpole= two*(-CF - (3*CF*epinv)/2. - CF*epinv**2 + (7*CF*Pi**2)/12. - 
     -  CF*epinv*Log(musq/s56) - (CF*Log(musq/s56)**2)/2._dp)
!     UV counter term CDR (inc in vpole
      msquv=zip
      endif
      
!---- int dips
      xl12=log(musq/(+two*dot(p,i1,i2)))
      msq_z=cf*dipHD_qq(xl12)*two  


      msqv=(msquv+vpole+msq_z)*msqlo*ason2pi

      return
      end


      subroutine Hbb_decay_z(p,i1,i2,msq_z)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'scale.f'
      include 'qcdcouple.f' 
      include 'nf.f'
      include 'epinv.f'
      include 'epinv2.f'
      real(dp) :: p(mxpart,4),msq_z
      integer i1,i2,i3,i
      real(dp)::ff_qq,ff_gg,ff_gq,fi_gq
      real(dp)::hbbgidips(5),hdecay
      real(dp):: xl12,xl13,xl23
      real(dp):: dipHD_qq,dipHD_qg,dipHD_gg
      real(dp):: dot
      msq_z=zip
      
!      call Hbbg_decay(p,i1,i2,i3,hdecay)
      
      hbbgidips(:)=zip

    
      
      
      return 

      end
