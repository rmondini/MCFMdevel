!---- routine for
!-----q(i1)+qb(i2)+Q(i4)+QB(i5)+gamma(3)
!-----amplitude A is when gamma couples to i1,i2 line
!-----amplitude B is when gamma couples to i4,i5 line

!----- in MCFM notation gamma = 3
      subroutine ga_qqbQQb(i1,i2,i4,i5,za,zb,qqb_A,QQb_B)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'nf.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      integer i1,i2,i3,i4,i5
      complex(dp):: qqb_A(2,2,2),QQB_b(2,2,2)
      complex(dp):: qqbQQB_gaamp
      integer h1,h2,h3
      real(dp):: ampsq,qa,qb,Bigagam
      complex(dp):: temp
      integer i,j

      i=5
      j=-4

      qa=Q(i)
      qb=Q(j)
!     
!----- photon is 3 
      i3=3
      qqb_A(1,1,2)=qqbQQB_gaamp(i1,i2,i4,i5,i3,za,zb)
      qqb_B(1,1,2)=qqbQQB_gaamp(i4,i5,i1,i2,i3,za,zb)
!
      qqb_A(2,2,1)=-qqbQQB_gaamp(i1,i2,i4,i5,i3,zb,za)
      qqb_B(2,2,1)=-qqbQQB_gaamp(i4,i5,i1,i2,i3,zb,za)
!      
      qqb_A(1,2,2)=qqbQQB_gaamp(i1,i2,i5,i4,i3,za,zb)
      qqb_B(1,2,2)=-qqbQQB_gaamp(i5,i4,i1,i2,i3,za,zb)
!     
      qqb_A(2,1,1)=-qqbQQB_gaamp(i1,i2,i5,i4,i3,zb,za)
      qqb_B(2,1,1)=qqbQQB_gaamp(i5,i4,i1,i2,i3,zb,za)
!     
      qqb_A(2,2,2)=-qqbQQB_gaamp(i2,i1,i5,i4,i3,za,zb)
      qqb_B(2,2,2)=-qqbQQB_gaamp(i5,i4,i2,i1,i3,za,zb)
!
      qqb_A(1,1,1)=qqbQQB_gaamp(i2,i1,i5,i4,i3,zb,za)
      qqb_B(1,1,1)=qqbQQB_gaamp(i5,i4,i2,i1,i3,zb,za)
!
      qqb_A(2,1,2)=-qqbQQB_gaamp(i2,i1,i4,i5,i3,za,zb)
      qqb_B(2,1,2)=qqbQQB_gaamp(i4,i5,i2,i1,i3,za,zb)
!
      qqb_A(1,2,1)=qqbQQB_gaamp(i2,i1,i4,i5,i3,zb,za)
      qqb_B(1,2,1)=-qqbQQB_gaamp(i4,i5,i2,i1,i3,zb,za)

      return
      ampsq=zip
      write(6,*) i1,i2,i3,i4,i5
      do h1=1,2
         do h2=1,2
            do h3=1,2
               temp=qa*qqb_A(h1,h2,h3)+qb*qqb_B(h1,h2,h3)
!               write(6,*) h1,h2,h3,temp
!               write(6,*) h1,h2,h3,qqb_B(h1,h2,h3)
               ampsq=ampsq+real(temp*conjg(temp),dp)
            enddo
         enddo
      enddo
      
      write(6,*) ampsq,Bigagam(i1,i4,i2,i5,i3,i,j)
      write(6,*) ampsq/Bigagam(i1,i4,i2,i5,i3,i,j)*8._dp*xn*CF*gsq**2*esq
 !     pause
      
      return
      end

      
!---- i1^- i2^+ i3^- i4^+ ga + (ga couples to i1,i2) line
      function qqbQQB_gaamp(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: qqbQQB_gaamp
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer i1,i2,i3,i4,i5
   
      qqbQQB_gaamp=-(za(i1,i3)**2/(za(i1,i5)*za(i2,i5)*za(i3,i4)))
      return
      end
