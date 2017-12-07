      subroutine ampsq_1gam1g2q(j1,j2,j3,j4,j5,j6,za,zb,
     & ampsq,ampsqid)
      implicit none
      include 'types.f'
c--- Returns matrix element squared for the process
c---     0  -->  qb(j1) + q(j2) + qbar(j3) + q(j4) + g(j5) + gam(j6)
c---
c--- returns arrays of matrix elements for non-identical quarks (ampsq)
c--- and for identical quarks (ampsqid), indexed by charges of
c--- (j1,j2) and (j3,j4) quark lines
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
      integer:: j1,j2,j3,j4,j5,j6,h1,h2,h3,h4,h5,j12,j34
      complex(dp):: amp(2,2,2,2,2,2),amp_swap(2,2,2,2,2,2),
     & A10(2,2,2,2,2,2),A01(2,2,2,2,2,2),
     & B10(2,2,2,2,2,2),B01(2,2,2,2,2,2),
     & A10_swap(2,2,2,2,2,2),A01_swap(2,2,2,2,2,2),
     & B10_swap(2,2,2,2,2,2),B01_swap(2,2,2,2,2,2)
      real(dp):: Q12,Q34,ampsq(2,2),ampsqid(2),bit,sign
              
c--- ordering of labels in amp is as follows:
c---  (charge of quark line 12, charge of quark line 34,
c---   helicity of antiquark j1, antiquark j3, gluon j5, photon j6;
c---   quark helicity j2=3-j1, quark j4=3-j3)

c--- loop over possible charges of quark lines
      do j12=1,2
      do j34=1,2

      if (j12 == 1) then
        Q12=-1._dp/3._dp
      else
        Q12=+2._dp/3._dp
      endif

      if (j34 == 1) then
        Q34=-1._dp/3._dp
      else
        Q34=+2._dp/3._dp
      endif

c--- basic amplitudes
      call amp_1gam1g2q_mpmppp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34,
     & A10(j12,j34,1,1,2,2),A01(j12,j34,1,1,2,2),
     & B10(j12,j34,1,1,2,2),B01(j12,j34,1,1,2,2))
     
      call amp_1gam1g2q_mppmpp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34,
     & A10(j12,j34,1,2,2,2),A01(j12,j34,1,2,2,2),
     & B10(j12,j34,1,2,2,2),B01(j12,j34,1,2,2,2))
     
      call amp_1gam1g2q_pmmppp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34,
     & A10(j12,j34,2,1,2,2),A01(j12,j34,2,1,2,2),
     & B10(j12,j34,2,1,2,2),B01(j12,j34,2,1,2,2))
     
      call amp_1gam1g2q_pmpmpp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34,
     & A10(j12,j34,2,2,2,2),A01(j12,j34,2,2,2,2),
     & B10(j12,j34,2,2,2,2),B01(j12,j34,2,2,2,2))
     

      call amp_1gam1g2q_pmpmmp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34,
     & A10(j12,j34,2,2,1,2),A01(j12,j34,2,2,1,2),
     & B10(j12,j34,2,2,1,2),B01(j12,j34,2,2,1,2))
     
      call amp_1gam1g2q_pmmpmp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34,
     & A10(j12,j34,2,1,1,2),A01(j12,j34,2,1,1,2),
     & B10(j12,j34,2,1,1,2),B01(j12,j34,2,1,1,2))
     

c--- simple symmetry
!      call amp_1gam1g2q_pmpmmp(j1,j2,j3,j4,j6,j5,za,zb,Q12,Q34,
!     & A10(j12,j34,2,2,2,1),A01(j12,j34,2,2,2,1),
!     & B10(j12,j34,2,2,2,1),B01(j12,j34,2,2,2,1))
     
!      call amp_1gam1g2q_pmmpmp(j1,j2,j3,j4,j6,j5,za,zb,Q12,Q34,
!     & A10(j12,j34,2,1,2,1),A01(j12,j34,2,1,2,1),
!     & B10(j12,j34,2,1,2,1),B01(j12,j34,2,1,2,1))

c--- new method
      call amp_1gam1g2q_mpmpmp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34,
     & A10(j12,j34,1,1,1,2),A01(j12,j34,1,1,1,2),
     & B10(j12,j34,1,1,1,2),B01(j12,j34,1,1,1,2))
      
!      call amp_1gam1g2q_mppmmp(j1,j2,j3,j4,j5,j6,za,zb,Q12,Q34,
!     & A10(j12,j34,1,2,1,2),A01(j12,j34,1,2,1,2),
!     & B10(j12,j34,1,2,1,2),B01(j12,j34,1,2,1,2))
      
!      call amp_1gam1g2q_pmpmmp(j2,j1,j4,j3,j5,j6,za,zb,Q12,Q34,
!     & A10(j12,j34,1,1,1,2),A01(j12,j34,1,1,1,2),
!     & B10(j12,j34,1,1,1,2),B01(j12,j34,1,1,1,2))
      
      call amp_1gam1g2q_pmmpmp(j2,j1,j4,j3,j5,j6,za,zb,Q12,Q34,
     & A10(j12,j34,1,2,1,2),A01(j12,j34,1,2,1,2),
     & B10(j12,j34,1,2,1,2),B01(j12,j34,1,2,1,2))
      

c--- remaining amplitudes by complex conjugation
      do h1=1,2
      do h2=1,2
      A10(j12,j34,h1,h2,1,1)=conjg(A10(j12,j34,3-h1,3-h2,2,2))
      A01(j12,j34,h1,h2,1,1)=conjg(A01(j12,j34,3-h1,3-h2,2,2))
      B10(j12,j34,h1,h2,1,1)=conjg(B10(j12,j34,3-h1,3-h2,2,2))
      B01(j12,j34,h1,h2,1,1)=conjg(B01(j12,j34,3-h1,3-h2,2,2))
      enddo
      enddo
      A10(j12,j34,1,1,2,1)=conjg(A10(j12,j34,2,2,1,2))
      A10(j12,j34,1,2,2,1)=conjg(A10(j12,j34,2,1,1,2))
      A10(j12,j34,2,2,2,1)=conjg(A10(j12,j34,1,1,1,2))
      A10(j12,j34,2,1,2,1)=conjg(A10(j12,j34,1,2,1,2))
          
      A01(j12,j34,1,1,2,1)=conjg(A01(j12,j34,2,2,1,2))
      A01(j12,j34,1,2,2,1)=conjg(A01(j12,j34,2,1,1,2))
      A01(j12,j34,2,2,2,1)=conjg(A01(j12,j34,1,1,1,2))
      A01(j12,j34,2,1,2,1)=conjg(A01(j12,j34,1,2,1,2))
          
      B10(j12,j34,1,1,2,1)=conjg(B10(j12,j34,2,2,1,2))
      B10(j12,j34,1,2,2,1)=conjg(B10(j12,j34,2,1,1,2))
      B10(j12,j34,2,2,2,1)=conjg(B10(j12,j34,1,1,1,2))
      B10(j12,j34,2,1,2,1)=conjg(B10(j12,j34,1,2,1,2))
          
      B01(j12,j34,1,1,2,1)=conjg(B01(j12,j34,2,2,1,2))
      B01(j12,j34,1,2,2,1)=conjg(B01(j12,j34,2,1,1,2))
      B01(j12,j34,2,2,2,1)=conjg(B01(j12,j34,1,1,1,2))
      B01(j12,j34,2,1,2,1)=conjg(B01(j12,j34,1,2,1,2))

c--- additional amplitudes needed for identical quarks
c--- (obtained by interchanging j2 and j4)
      if (j12. eq. j34) then
c--- basic amplitudes
      call amp_1gam1g2q_mpmppp(j1,j4,j3,j2,j5,j6,za,zb,Q12,Q34,
     & A10_swap(j12,j34,1,1,2,2),A01_swap(j12,j34,1,1,2,2),
     & B10_swap(j12,j34,1,1,2,2),B01_swap(j12,j34,1,1,2,2))
     
      call amp_1gam1g2q_mppmpp(j1,j4,j3,j2,j5,j6,za,zb,Q12,Q34,
     & A10_swap(j12,j34,1,2,2,2),A01_swap(j12,j34,1,2,2,2),
     & B10_swap(j12,j34,1,2,2,2),B01_swap(j12,j34,1,2,2,2))
     
      call amp_1gam1g2q_pmmppp(j1,j4,j3,j2,j5,j6,za,zb,Q12,Q34,
     & A10_swap(j12,j34,2,1,2,2),A01_swap(j12,j34,2,1,2,2),
     & B10_swap(j12,j34,2,1,2,2),B01_swap(j12,j34,2,1,2,2))
     
      call amp_1gam1g2q_pmpmpp(j1,j4,j3,j2,j5,j6,za,zb,Q12,Q34,
     & A10_swap(j12,j34,2,2,2,2),A01_swap(j12,j34,2,2,2,2),
     & B10_swap(j12,j34,2,2,2,2),B01_swap(j12,j34,2,2,2,2))
     

      call amp_1gam1g2q_pmpmmp(j1,j4,j3,j2,j5,j6,za,zb,Q12,Q34,
     & A10_swap(j12,j34,2,2,1,2),A01_swap(j12,j34,2,2,1,2),
     & B10_swap(j12,j34,2,2,1,2),B01_swap(j12,j34,2,2,1,2))
     
      call amp_1gam1g2q_pmmpmp(j1,j4,j3,j2,j5,j6,za,zb,Q12,Q34,
     & A10_swap(j12,j34,2,1,1,2),A01_swap(j12,j34,2,1,1,2),
     & B10_swap(j12,j34,2,1,1,2),B01_swap(j12,j34,2,1,1,2))
     

c--- simple symmetry
!      call amp_1gam1g2q_pmpmmp(j1,j4,j3,j2,j6,j5,za,zb,Q12,Q34,
!     & A10_swap(j12,j34,2,2,2,1),A01_swap(j12,j34,2,2,2,1),
!     & B10_swap(j12,j34,2,2,2,1),B01_swap(j12,j34,2,2,2,1))
     
!      call amp_1gam1g2q_pmmpmp(j1,j4,j3,j2,j6,j5,za,zb,Q12,Q34,
!     & A10_swap(j12,j34,2,1,2,1),A01_swap(j12,j34,2,1,2,1),
!     & B10_swap(j12,j34,2,1,2,1),B01_swap(j12,j34,2,1,2,1))

! new method
      call amp_1gam1g2q_mpmpmp(j1,j4,j3,j2,j5,j6,za,zb,Q12,Q34,
     & A10_swap(j12,j34,1,1,1,2),A01_swap(j12,j34,1,1,1,2),
     & B10_swap(j12,j34,1,1,1,2),B01_swap(j12,j34,1,1,1,2))
      
!      call amp_1gam1g2q_pmpmmp(j4,j1,j2,j3,j5,j6,za,zb,Q12,Q34,
!     & A10_swap(j12,j34,1,1,1,2),A01_swap(j12,j34,1,1,1,2),
!     & B10_swap(j12,j34,1,1,1,2),B01_swap(j12,j34,1,1,1,2))

      call amp_1gam1g2q_pmmpmp(j4,j1,j2,j3,j5,j6,za,zb,Q12,Q34,
     & A10_swap(j12,j34,1,2,1,2),A01_swap(j12,j34,1,2,1,2),
     & B10_swap(j12,j34,1,2,1,2),B01_swap(j12,j34,1,2,1,2))
     
!      A10_swap(:,:,1,1,1,2)=-A10_swap(:,:,1,1,1,2)
!      A10_swap(:,:,2,2,2,2)=-A10_swap(:,:,2,2,2,2)
!      A01_swap(:,:,1,1,1,2)=-A01_swap(:,:,1,1,1,2)
!      A01_swap(:,:,2,2,2,2)=-A01_swap(:,:,2,2,2,2)
!      B10_swap(:,:,1,1,1,2)=-B10_swap(:,:,1,1,1,2)
!      B10_swap(:,:,2,2,2,2)=-B10_swap(:,:,2,2,2,2)
!      B01_swap(:,:,1,1,1,2)=-B01_swap(:,:,1,1,1,2)
!      B01_swap(:,:,2,2,2,2)=-B01_swap(:,:,2,2,2,2)
      
      
c--- remaining amplitudes by complex conjugation
      do h1=1,2
      do h2=1,2
      A10_swap(j12,j34,h1,h2,1,1)=conjg(A10_swap(j12,j34,3-h1,3-h2,2,2))
      A01_swap(j12,j34,h1,h2,1,1)=conjg(A01_swap(j12,j34,3-h1,3-h2,2,2))
      B10_swap(j12,j34,h1,h2,1,1)=conjg(B10_swap(j12,j34,3-h1,3-h2,2,2))
      B01_swap(j12,j34,h1,h2,1,1)=conjg(B01_swap(j12,j34,3-h1,3-h2,2,2))
      enddo
      enddo
      A10_swap(j12,j34,1,1,2,1)=conjg(A10_swap(j12,j34,2,2,1,2))
      A10_swap(j12,j34,1,2,2,1)=conjg(A10_swap(j12,j34,2,1,1,2))
      A10_swap(j12,j34,2,2,2,1)=conjg(A10_swap(j12,j34,1,1,1,2))
      A10_swap(j12,j34,2,1,2,1)=conjg(A10_swap(j12,j34,1,2,1,2))
          
      A01_swap(j12,j34,1,1,2,1)=conjg(A01_swap(j12,j34,2,2,1,2))
      A01_swap(j12,j34,1,2,2,1)=conjg(A01_swap(j12,j34,2,1,1,2))
      A01_swap(j12,j34,2,2,2,1)=conjg(A01_swap(j12,j34,1,1,1,2))
      A01_swap(j12,j34,2,1,2,1)=conjg(A01_swap(j12,j34,1,2,1,2))
          
      B10_swap(j12,j34,1,1,2,1)=conjg(B10_swap(j12,j34,2,2,1,2))
      B10_swap(j12,j34,1,2,2,1)=conjg(B10_swap(j12,j34,2,1,1,2))
      B10_swap(j12,j34,2,2,2,1)=conjg(B10_swap(j12,j34,1,1,1,2))
      B10_swap(j12,j34,2,1,2,1)=conjg(B10_swap(j12,j34,1,2,1,2))
          
      B01_swap(j12,j34,1,1,2,1)=conjg(B01_swap(j12,j34,2,2,1,2))
      B01_swap(j12,j34,1,2,2,1)=conjg(B01_swap(j12,j34,2,1,1,2))
      B01_swap(j12,j34,2,2,2,1)=conjg(B01_swap(j12,j34,1,1,1,2))
      B01_swap(j12,j34,2,1,2,1)=conjg(B01_swap(j12,j34,1,2,1,2))

      endif
      
      enddo
      enddo
 
c--- square amplitudes
      do j12=1,2
      ampsqid(j12)=0._dp
      do j34=1,2
      ampsq(j12,j34)=0._dp
      
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      bit=real(
     & +A10(j12,j34,h1,h2,h3,h4)*conjg(A10(j12,j34,h1,h2,h3,h4))
     & +A01(j12,j34,h1,h2,h3,h4)*conjg(A01(j12,j34,h1,h2,h3,h4))
     & +B10(j12,j34,h1,h2,h3,h4)*conjg(B10(j12,j34,h1,h2,h3,h4))/xnsq
     & +B01(j12,j34,h1,h2,h3,h4)*conjg(B01(j12,j34,h1,h2,h3,h4))/xnsq
     & -two/xnsq*(A10(j12,j34,h1,h2,h3,h4)+A01(j12,j34,h1,h2,h3,h4))
     &     *conjg(B10(j12,j34,h1,h2,h3,h4)+B01(j12,j34,h1,h2,h3,h4)),dp)
      ampsq(j12,j34)=ampsq(j12,j34)+bit
      if (j12 == j34) then
      bit=bit+real(
     & +A10_swap(j12,j34,h1,h2,h3,h4)*conjg(A10_swap(j12,j34,h1,h2,h3,h4))
     & +A01_swap(j12,j34,h1,h2,h3,h4)*conjg(A01_swap(j12,j34,h1,h2,h3,h4))
     & +B10_swap(j12,j34,h1,h2,h3,h4)*conjg(B10_swap(j12,j34,h1,h2,h3,h4))/xnsq
     & +B01_swap(j12,j34,h1,h2,h3,h4)*conjg(B01_swap(j12,j34,h1,h2,h3,h4))/xnsq
     & -two/xnsq*(A10_swap(j12,j34,h1,h2,h3,h4)+A01_swap(j12,j34,h1,h2,h3,h4))
     &     *conjg(B10_swap(j12,j34,h1,h2,h3,h4)+B01_swap(j12,j34,h1,h2,h3,h4)),dp)
        if (h1 == h2) then
! Need to add the corresponding bit here (flip sign of first two to get others right)
! do not understand this sign
        if (h3 == h4) then
          sign=-one
        else
          sign=+one
        endif
        bit=bit+real(
     & -(sign)*two/xn*(A10(j12,j34,h1,h2,h3,h4)+A01(j12,j34,h1,h2,h3,h4))
     &   *conjg(A10_swap(j12,j34,h1,h2,h3,h4)+A01_swap(j12,j34,h1,h2,h3,h4))
     & +(sign)*two/xn*(A10(j12,j34,h1,h2,h3,h4)*conjg(B10_swap(j12,j34,h1,h2,h3,h4))
     &         +A01(j12,j34,h1,h2,h3,h4)*conjg(B01_swap(j12,j34,h1,h2,h3,h4))
     &         +A10_swap(j12,j34,h1,h2,h3,h4)*conjg(B10(j12,j34,h1,h2,h3,h4))
     &         +A01_swap(j12,j34,h1,h2,h3,h4)*conjg(B01(j12,j34,h1,h2,h3,h4)))
     & -two/xn**3*(B10(j12,j34,h1,h2,h3,h4)+B01(j12,j34,h1,h2,h3,h4))
     &      *conjg(B10_swap(j12,j34,h1,h2,h3,h4)+B01_swap(j12,j34,h1,h2,h3,h4)),dp)
        endif
        ampsqid(j12)=ampsqid(j12)+bit
      endif
 !     if ((j12 == 2) .and. (j34 == 2) .and. (h1 .eq. h2)) then
 !       if (bit .ne. 0d0) write(6,*) h1,h2,h3,h4,'   ',half*aveqq*bit*4d0*esq*gsq**3*xn**2*CF
 !     elseif ((j12 == 2) .and. (j34 == 2)) then
 !       write(6,*) h1,h2,h3,h4,'(a)',half*aveqq*4d0*esq*gsq**3*xn**2*CF*real(
 !    & +A10(j12,j34,h1,h2,h3,h4)*conjg(A10(j12,j34,h1,h2,h3,h4))
 !    & +A01(j12,j34,h1,h2,h3,h4)*conjg(A01(j12,j34,h1,h2,h3,h4))
 !    & +B10(j12,j34,h1,h2,h3,h4)*conjg(B10(j12,j34,h1,h2,h3,h4))/xnsq
 !    & +B01(j12,j34,h1,h2,h3,h4)*conjg(B01(j12,j34,h1,h2,h3,h4))/xnsq
 !    & -two/xnsq*(A10(j12,j34,h1,h2,h3,h4)+A01(j12,j34,h1,h2,h3,h4))
 !    &     *conjg(B10(j12,j34,h1,h2,h3,h4)+B01(j12,j34,h1,h2,h3,h4)),dp)
 !       write(6,*) h1,h2,h3,h4,'(b)',half*aveqq*4d0*esq*gsq**3*xn**2*CF*real(
 !    & +A10_swap(j12,j34,h1,h2,h3,h4)*conjg(A10_swap(j12,j34,h1,h2,h3,h4))
 !    & +A01_swap(j12,j34,h1,h2,h3,h4)*conjg(A01_swap(j12,j34,h1,h2,h3,h4))
 !    & +B10_swap(j12,j34,h1,h2,h3,h4)*conjg(B10_swap(j12,j34,h1,h2,h3,h4))/xnsq
 !    & +B01_swap(j12,j34,h1,h2,h3,h4)*conjg(B01_swap(j12,j34,h1,h2,h3,h4))/xnsq
 !    & -two/xnsq*(A10_swap(j12,j34,h1,h2,h3,h4)+A01_swap(j12,j34,h1,h2,h3,h4))
 !    &     *conjg(B10_swap(j12,j34,h1,h2,h3,h4)+B01_swap(j12,j34,h1,h2,h3,h4)),dp)
 !     endif
      enddo
      enddo
      enddo
      enddo
      
      enddo
      enddo
      
      return
      end
      
