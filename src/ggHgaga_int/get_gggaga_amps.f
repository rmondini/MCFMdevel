! This routine returns the one-loop amplitudes for g+g -> ga+ga,
! optionally including the effect of the top-quark mass
!
! It is based on the routine gggaga_mt, which computes the amplitude squared
      subroutine get_gggaga_amps(p,massless,amp_tot)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'nf.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'sprods_com.f'
      real(dp):: p(mxpart,4),fac_top,fac_light,fac,Qsum,
     & fac_t,fac_b,fac_c
      complex(dp):: amp_massless(2,2,2,2),amp_mass(2,2,2,2)
      complex(dp):: amp_t(2,2,2,2),amp_b(2,2,2,2),amp_c(2,2,2,2)
      complex(dp):: amp_tot(2,2,2,2)
      logical massless

!!! NEW IMPLEMENTATION OF t,b,c MASSIVE LOOPS

      amp_massless(:,:,:,:)=czip
      amp_t(:,:,:,:)=czip
      amp_b(:,:,:,:)=czip
      amp_c(:,:,:,:)=czip
     
      Qsum=Q(1)**2+Q(2)**2+Q(3)**2
    
      fac_light=Qsum*2._dp*ason2pi*esq
      fac_t=Q(2)**2*2._dp*ason2pi*esq
      fac_b=Q(1)**2*2._dp*ason2pi*esq
      fac_c=fac_t
      
      call spinoru(4,p,za,zb)

!        call fill_amp_gggaga_mass(1,2,3,4,za,zb,amp_t,mt**2)
!        call fill_amp_gggaga_mass(1,2,3,4,za,zb,amp_b,mbsq)
!        call fill_amp_gggaga_mass(1,2,3,4,za,zb,amp_c,mcsq)
      if (massless) then
        continue
      else
        call fill_amp_gggaga_mass(1,2,3,4,za,zb,amp_t,168.2_dp**2)
        call fill_amp_gggaga_mass(1,2,3,4,za,zb,amp_b,2.78_dp**2)
        call fill_amp_gggaga_mass(1,2,3,4,za,zb,amp_c,0.72_dp**2)
      endif
      
      call fill_amp_gggaga_mass(1,2,3,4,za,zb,amp_massless,zip)
      
      amp_tot(:,:,:,:)=
     &      fac_light*amp_massless(:,:,:,:)
     &     +fac_t*amp_t(:,:,:,:)
     &     +fac_b*amp_b(:,:,:,:)
     &     +fac_c*amp_c(:,:,:,:)

      if (massless) then
! Only return massless loops
        amp_tot(:,:,:,:)=
     &        (Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2+Q(5)**2)
     &       *2._dp*ason2pi*esq*amp_massless(:,:,:,:)
      endif
      
      amp_tot=amp_tot*im

      return
      
!!! END NEW IMPLEMENTATION

      amp_massless(:,:,:,:)=czip
      amp_mass(:,:,:,:)=czip
     
      Qsum=Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2+Q(5)**2
    
      fac_light=Qsum*2._dp*ason2pi*esq
      fac_top=Q(2)**2*2._dp*ason2pi*esq
      
      call spinoru(4,p,za,zb)

!      if (includetop) then
        call fill_amp_gggaga_mass(1,2,3,4,za,zb,amp_mass,mt**2)
!      endif
      
      call fill_amp_gggaga_mass(1,2,3,4,za,zb,amp_massless,zip)
      
      amp_tot(:,:,:,:)=
     &      fac_light*amp_massless(:,:,:,:)
     &     +fac_top*amp_mass(:,:,:,:)

      amp_tot=amp_tot*im

      return
      end

