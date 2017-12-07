! This routine returns the approximate two-loop amplitudes for g+g -> H -> ga+ga
!
! It multiplies the one-loop amplitude that is supplied by the correction
! factor that is obtained in the effective theory (mt -> infinity)
      subroutine get_ggHgaga_amps_twoloop(amp_1loop,amp_2loop)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'nf.f'
      include 'masses.f'
      include 'nflav.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'scheme.f'
      real(dp):: ss
      complex(dp):: amp_1loop(2,2,2,2),amp_2loop(2,2,2,2),I1ggtogamgam,
     & lnrat,xlog

      ss=s(1,2)
      xlog=lnrat(musq,-ss)
c--- This is taken from hep-ph/0109078 Eq.(2.11); note however that the log
c--- proportional to the beta-function coefficient is added in Eq. (2.11)
c--- and subtracted again in Eq. (4.5), therefore we omit it here.
      I1ggtogamgam=
     & -xn*((epinv**2+epinv*xlog+0.5_dp*xlog**2)
     & +(11._dp/6._dp-real(nflav,dp)/(3._dp*xn))*epinv)

! This factor has been derived by taking the difference between I1ggtogamgam
! given above (i.e. from Eq. (2.11) and the log in Eq. (4.5)) and the constant
! in gg_h_v.f, after factors of im*pi in that expression have been restored
      amp_2loop(:,:,:,:)=amp_1loop(:,:,:,:)
     & *ason2pi*(I1ggtogamgam+xn*two)
     
      return
      end


