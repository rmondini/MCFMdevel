      function ampsq_1gam3g(i1,i2,i3,j4,j5,j6,za,zb)
      implicit none
      include 'types.f'
      real(dp):: ampsq_1gam3g
c--- Matrix element squared for the process
c---       0  -->  qb(j5) + q(j6) + g(i1) + g(i2) + g(i3) + gam(j4)
c---
c--- Taken from "Multi-Photon Amplitudes for Next-to-Leading Order QCD"
c---  V. Del Duca, W. Kilgore and F. Maltoni, hep-ph/9910253
c---
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      integer:: i1,i2,i3,j1,j2,j3,j4,j5,j6,h1,h2,h3,h4,h5,icol,ii(3),h(3)
      complex(dp):: amp(6,2,2,2,2,2),m(6),tempm0,
     & amp_1gam3g_mppppm,amp_1gam3g_pmpppm,
     & amp_1gam3g_ppmppm,amp_1gam3g_pppmpm,amp_1gam3g_mmpppm,
     & amp_1gam3g_pmmppm,amp_1gam3g_mpmppm
      integer, parameter:: ii1(6)=(/1,1,2,2,3,3/)
      integer, parameter:: ii2(6)=(/2,3,1,3,1,2/)
      integer, parameter:: ii3(6)=(/3,2,3,1,2,1/)
      
      ii(1)=i1
      ii(2)=i2
      ii(3)=i3
      
c--- loop over color ordering of gluons
      do icol=1,6
      
      j1=ii(ii1(icol))
      j2=ii(ii2(icol))
      j3=ii(ii3(icol))
        
c--- ordering of labels in amp is as follows:
c---  (gluon j1, gluon j2, gluon j3, photon j4, antiquark j5)
c--- for consistency with function names (and quark j6=3-j5)

c--- trivial amplitude
      amp(icol,2,2,2,2,2)=czip
c--- basic MHV amplitudes
      amp(icol,1,2,2,2,2)=amp_1gam3g_mppppm(j1,j2,j3,j4,j5,j6,za,zb)
      amp(icol,2,1,2,2,2)=amp_1gam3g_pmpppm(j1,j2,j3,j4,j5,j6,za,zb)
      amp(icol,2,2,1,2,2)=amp_1gam3g_ppmppm(j1,j2,j3,j4,j5,j6,za,zb)
      amp(icol,2,2,2,1,2)=amp_1gam3g_pppmpm(j1,j2,j3,j4,j5,j6,za,zb)

c--- non-MHV amplitudes
      amp(icol,1,1,2,2,2)=amp_1gam3g_mmpppm(j1,j2,j3,j4,j5,j6,za,zb)
      amp(icol,2,1,1,2,2)=amp_1gam3g_pmmppm(j1,j2,j3,j4,j5,j6,za,zb)
      amp(icol,1,2,1,2,2)=amp_1gam3g_mpmppm(j1,j2,j3,j4,j5,j6,za,zb)

c--- parity and charge conjugation
c--- (complex conjugate, interchange q, qb)
      amp(icol,1,1,1,1,2)=czip
      
      amp(icol,2,1,1,1,2)=-amp_1gam3g_mppppm(j1,j2,j3,j4,j6,j5,zb,za)
      amp(icol,1,2,1,1,2)=-amp_1gam3g_pmpppm(j1,j2,j3,j4,j6,j5,zb,za)
      amp(icol,1,1,2,1,2)=-amp_1gam3g_ppmppm(j1,j2,j3,j4,j6,j5,zb,za)
      amp(icol,1,1,1,2,2)=-amp_1gam3g_pppmpm(j1,j2,j3,j4,j6,j5,zb,za)
     
      amp(icol,2,2,1,1,2)=-amp_1gam3g_mmpppm(j1,j2,j3,j4,j6,j5,zb,za)
      amp(icol,1,2,2,1,2)=-amp_1gam3g_pmmppm(j1,j2,j3,j4,j6,j5,zb,za)
      amp(icol,2,1,2,1,2)=-amp_1gam3g_mpmppm(j1,j2,j3,j4,j6,j5,zb,za)

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      amp(icol,h1,h2,h3,h4,1)=amp(icol,3-h1,3-h2,3-h3,3-h4,2)
!      write(6,*) h1,h2,h3,h4,1,amp(icol,h1,h2,h3,h4,1)
!      write(6,*) h1,h2,h3,h4,2,amp(icol,h1,h2,h3,h4,2)
      enddo
      enddo
      enddo
      enddo
      
      enddo
      
!      do h1=1,2
!      do h2=1,2
!      do h3=1,2
!      do h4=1,2
!      amp(1,h1,h2,h3,h4,1)=-amp(6,h3,h2,h1,h4,2)
!      amp(2,h1,h2,h3,h4,1)=-amp(4,h3,h2,h1,h4,2)
!      amp(3,h1,h2,h3,h4,1)=-amp(5,h3,h2,h1,h4,2)
!      amp(4,h1,h2,h3,h4,1)=-amp(2,h3,h2,h1,h4,2)
!      amp(5,h1,h2,h3,h4,1)=-amp(3,h3,h2,h1,h4,2)
!      amp(6,h1,h2,h3,h4,1)=-amp(1,h3,h2,h1,h4,2)
!      enddo
!      enddo
!      enddo
!      enddo

c--- note: obvious redundancy in this routine, but might be
c--- worth checking relations for use in virtual
      
      ampsq_1gam3g=0._dp
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do h5=1,2
      
      h(1)=h1
      h(2)=h2
      h(3)=h3

c--- ensure all color permutations have helicities permuted appropriately too   
      do icol=1,6
      m(icol)=amp(icol,h(ii1(icol)),h(ii2(icol)),h(ii3(icol)),h4,h5)
      enddo
      
      tempm0=zip
      do icol=1,6
c--- leading color contribution
      ampsq_1gam3g=ampsq_1gam3g+real(m(icol)*conjg(m(icol)),dp)
      tempm0=tempm0+m(icol)
      enddo    

c--- subleading color
      ampsq_1gam3g=ampsq_1gam3g-real(
     &  (m(1)+m(3)+m(4))*conjg(m(1)+m(3)+m(4))
     & +(m(2)+m(5)+m(6))*conjg(m(2)+m(5)+m(6))
     & +(m(1)+m(2)+m(3))*conjg(m(1)+m(2)+m(3))
     & +(m(4)+m(5)+m(6))*conjg(m(4)+m(5)+m(6))
     & +(m(1)+m(2)+m(5))*conjg(m(1)+m(2)+m(5))
     & +(m(3)+m(4)+m(6))*conjg(m(3)+m(4)+m(6)),dp)/xnsq


c--- sub-subleading color
      ampsq_1gam3g=ampsq_1gam3g
     & +real(tempm0*conjg(tempm0),dp)*(xnsq+one)/xnsq**2

      enddo
      enddo
      enddo
      enddo
      enddo
      
!      write(6,*) 'ampsq_1gam3g',ampsq_1gam3g
!      pause
      
      return
      end
      
