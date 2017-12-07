! This routine returns the two-loop amplitudes for g+g -> ga+ga
!
! It is based on the routine virtgamgam, which computes the interference
! of one- and two-loop amplitudes
      subroutine get_gggaga_amps_twoloop(p,amp_1loop,amp_2loop,just1loop)
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
      include 'zprods_decl.f'
      include 'epinv.f'
      include 'scheme.f'
      integer:: h1,h2,h3,h4
      logical just1loop
      real(dp):: p(mxpart,4),fac,Qsum,ss,tt,uu
      complex(dp):: amp_1loop(2,2,2,2),amp_2loop(2,2,2,2),
     & m1(2,2,2,2),m2fin(2,2,2,2),
     & I1ggtogamgam,xlog,FL(2,2,2,2),FSL(2,2,2,2),lnrat,
     & phase(2,2,2,2)

      scheme='dred'
c--- initialize gg 2-loop matrix elements
      call spinoru(4,p,za,zb)

! These are the definitions in hep-ph/0109078 and hep-ph/0206194
! Note that they differ by the interchange of t and u from our normal choice
      ss=s(1,2)
      tt=s(2,3)
      uu=s(1,3)

      xlog=lnrat(musq,-ss)
c--- This is taken from hep-ph/0109078 Eq.(2.11); note however that the log
c--- proportional to the beta-function coefficient is added in Eq. (2.11)
c--- and subtracted again in Eq. (4.5), therefore we omit it here.
      I1ggtogamgam=
     & -xn*((epinv**2+epinv*xlog+0.5_dp*xlog**2)
     & +(11._dp/6._dp-real(nflav,dp)/(3._dp*xn))*epinv)

      call M1fill(ss,tt,uu,m1)

c--- testing the 1-loop matrix element contribution up to Order(ep^2)
c      call oneloopep(s,t,u,IxM1ep)

      if (just1loop) then
        FL=czip
        FSL=czip
      else
c--- pure 2-loop contribution
        call twoloop(ss,tt,uu,FL,FSL)
      endif
      
c--- Compare m1 with m1ep
c      do h1=1,2
c      do h2=1,2
c      do h3=1,2
c      do h4=1,2
c      write(6,*) h1,h2,h3,h4,m1(h1,h2,h3,h4)*I1ggtogamgam,
c     &           IxM1ep(h1,h2,h3,h4)
c      enddo
c      enddo
c      enddo
c      enddo
c      pause

! construct phases according to Eq. (3.5)
      phase(2,2,2,2)=im*zb(1,2)*zb(3,4)/(za(1,2)*za(3,4))
      phase(1,1,1,1)=im*za(1,2)*za(3,4)/(zb(1,2)*zb(3,4))
      
      phase(1,2,2,2)=im*za(1,2)*za(1,4)*zb(2,4)/(za(3,4)*za(2,3)*za(2,4))
      phase(2,1,1,1)=im*zb(1,2)*zb(1,4)*za(2,4)/(zb(3,4)*zb(2,3)*zb(2,4))
      
      phase(2,1,2,2)=im*za(2,3)*za(2,4)*zb(3,4)/(za(1,4)*za(3,1)*za(3,4))
      phase(1,2,1,1)=im*zb(2,3)*zb(2,4)*za(3,4)/(zb(1,4)*zb(3,1)*zb(3,4))
      
      phase(2,2,1,2)=im*za(3,2)*za(3,4)*zb(2,4)/(za(1,4)*za(2,1)*za(2,4))
      phase(1,1,2,1)=im*zb(3,2)*zb(3,4)*za(2,4)/(zb(1,4)*zb(2,1)*zb(2,4))

      phase(1,1,2,2)=im*za(1,2)*zb(3,4)/(zb(1,2)*za(3,4))
      phase(2,2,1,1)=im*zb(1,2)*za(3,4)/(za(1,2)*zb(3,4))
      
      phase(1,2,1,2)=im*za(1,3)*zb(2,4)/(zb(1,3)*za(2,4))
      phase(2,1,2,1)=im*zb(1,3)*za(2,4)/(za(1,3)*zb(2,4))
      
      phase(2,1,1,2)=im*za(2,3)*zb(1,4)/(zb(2,3)*za(1,4))
      phase(1,2,2,1)=im*zb(2,3)*za(1,4)/(za(2,3)*zb(1,4))
      
c--- the following amplitude has been added by me (permutation of above)
      phase(2,2,2,1)=im*za(4,2)*za(4,3)*zb(2,3)/(za(1,3)*za(2,1)*za(2,3))
      phase(1,1,1,2)=im*zb(4,2)*zb(4,3)*za(2,3)/(zb(1,3)*zb(2,1)*zb(2,3))

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      m2fin(h1,h2,h3,h4)=xn*FL(h1,h2,h3,h4)-FSL(h1,h2,h3,h4)/xn
      amp_1loop(h1,h2,h3,h4)=phase(h1,h2,h3,h4)*m1(h1,h2,h3,h4)
      amp_2loop(h1,h2,h3,h4)=phase(h1,h2,h3,h4)
     & *(m1(h1,h2,h3,h4)*I1ggtogamgam+m2fin(h1,h2,h3,h4))

      enddo
      enddo
      enddo
      enddo
      
c----remaining factor [avegg*V*statfac] is applied in wrapping routine
      Qsum=Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2+Q(5)**2 
      fac=two*Qsum*esq*ason2pi

      amp_1loop=amp_1loop*fac
      amp_2loop=amp_2loop*fac*ason2pi

      return
      end


