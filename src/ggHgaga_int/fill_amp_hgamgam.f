      subroutine fill_amp_hgamgam(sH,za,zb,k3,k4,amp)
c--- amplitudes for H(kH) -> ga(k3)+ga(k4) with kH^2 = sH
      implicit none
      include 'types.f'
      complex(dp):: amp(2,2)
      include 'constants.f'
      include 'nf.f'
      include 'cplx.h'
      include 'mxpart.f'
      include 'masses.f'
      include 'msbarmasses.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'couple.f'
      include 'kpart.f'
      include 'first.f'
      include 'zprods_decl.f'
      integer k3,k4
      complex(dp):: Iw,Iq,Ftriangle,bit
      real(dp):: prefac,sH
      real(dp):: x_t,x_b,x_w,x,mt_eff,mb_eff,massfrun
      save mt_eff,mb_eff
!$omp threadprivate(mt_eff,mb_eff)

C---statement functions
      Iq(x)=four*x*(ctwo+(four*x-one)*Ftriangle(x))
      Iw(x)=-ctwo*(cplx1(six*x+one)
     & +(six*x*(two*x-one))*Ftriangle(x))
C---end statement functions

      if (first) then
c--- run mt to appropriate scale
        if (kpart==klord) then
          mb_eff=massfrun(mb_msbar,hmass,amz,1)
          mt_eff=massfrun(mt_msbar,hmass,amz,1)
        else
          mb_eff=massfrun(mb_msbar,hmass,amz,2)
          mt_eff=massfrun(mt_msbar,hmass,amz,2)
        endif
        first=.false.
      endif

c---maybe it would be better to add esq at a higher scale.
      prefac=esq/(four*pi)*sH*sqrt(Gf/(eight*rt2*pi**2))
      prefac=prefac*rt2 ! statistical factor should be applied in wrapper routine
      x_b=mb_eff**2/sH
      x_t=mt_eff**2/sH
      x_w=wmass**2/sH
      bit=im*prefac*(xn*(Q(1)**2*Iq(x_b)+Q(2)**2*Iq(x_t))+Iw(x_w))
     
      amp(:,:)=czip
      amp(1,1)=bit*za(k3,k4)/zb(k3,k4)/rt2
      amp(2,2)=bit*zb(k3,k4)/za(k3,k4)/rt2

      return
      end

