! This routine returns the one-loop amplitudes for g+g -> ga+ga+g,
! with light quarks circulating in the loop
!
! It is based on the routine ggtogagag, which computes the amplitude squared
      subroutine get_gggaga_g_amps(p,za,zb,amp)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      real(dp):: fac,Qsumsq,p(mxpart,4)
      complex(dp):: Amp(2,2,2,2,2),
     & A1(2,2,2,2,2),A2(2,2,2,2,2),A3(2,2,2,2,2),
     & A4(2,2,2,2,2),A5(2,2,2,2,2),A6(2,2,2,2,2),
     & A7(2,2,2,2,2),A8(2,2,2,2,2),A9(2,2,2,2,2),
     & A10(2,2,2,2,2),A11(2,2,2,2,2),A12(2,2,2,2,2)
      integer:: h1,h2,h3,h4,h5
      
      call Aboxfill(1,2,3,4,5,za,zb,A1)
      call Aboxfill(1,2,4,3,5,za,zb,A2)
      call Aboxfill(1,3,2,4,5,za,zb,A3)
      call Aboxfill(3,1,2,4,5,za,zb,A4)

      call Aboxfill(1,4,2,3,5,za,zb,A5)
      call Aboxfill(3,1,4,2,5,za,zb,A6)
      call Aboxfill(1,3,4,2,5,za,zb,A7)
      call Aboxfill(1,4,3,2,5,za,zb,A8)

      call Aboxfill(4,1,2,3,5,za,zb,A9)
      call Aboxfill(4,1,3,2,5,za,zb,A10)
      call Aboxfill(3,4,1,2,5,za,zb,A11)
      call Aboxfill(4,3,1,2,5,za,zb,A12)

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do h5=1,2
      amp(h1,h2,h3,h4,h5)=
     &  A1(h1,h2,h3,h4,h5)+A2(h1,h2,h3,h4,h5)+A3(h1,h2,h3,h4,h5)
     & +A4(h1,h2,h3,h4,h5)+A5(h1,h2,h3,h4,h5)+A6(h1,h2,h3,h4,h5)
     & +A7(h1,h2,h3,h4,h5)+A8(h1,h2,h3,h4,h5)+A9(h1,h2,h3,h4,h5)
     & +A10(h1,h2,h3,h4,h5)+A11(h1,h2,h3,h4,h5)+A12(h1,h2,h3,h4,h5)
      enddo
      enddo
      enddo
      enddo
      enddo

      Qsumsq=Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2+Q(5)**2
! fac^2 below is taken from ggtogagag.f
!      fac^2=[Qsumsq*esq*sqrt(gsq)*ason2pi/2]^2 *xn/V
!           =[Qsumsq*esq*sqrt(gsq)*ason2pi*rt2]^2 *xn*V*avegg*statfac
      
      fac=Qsumsq*ason2pi*esq*sqrt(gsq)*rt2

c--- apply overall factor (note that Aboxfill gives -i*A)
      amp(:,:,:,:,:)=fac*amp(:,:,:,:,:)*(im)
 
      return
      end
