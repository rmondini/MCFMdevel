      function A51mpppp(j1,j2,j3,j4,j5,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: A51mpppp
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer:: j1,j2,j3,j4,j5
C----Eq.(4) of hep-ph/9302280v1 of BDK multiplied by 16*pi^2*(i)
c----in order to match the normalization of A51mmppp.f and A51mpmpp.f
c---- (note additional minus sign from the relation A[1/2]=-A[1])
C--- It therefore gives (16*pi^2)*(-i)*A^{[1/2]}_{5;1}
      A51mpppp=-((s(j2,j3)+s(j3,j4)+s(j4,j5))*zb(j2,j5)**2
     & -zb(j2,j4)*za(j4,j3)*zb(j3,j5)*zb(j2,j5)
     & -zb(j1,j2)*zb(j1,j5)/(za(j1,j2)*za(j1,j5))
     & *(za(j1,j2)**2*za(j1,j3)**2*zb(j2,j3)/za(j2,j3)
     &  +za(j1,j3)**2*za(j1,j4)**2*zb(j3,j4)/za(j3,j4)
     &  +za(j1,j4)**2*za(j1,j5)**2*zb(j4,j5)/za(j4,j5)))
     & /(3._dp*zb(j1,j2)*za(j2,j3)*za(j3,j4)*za(j4,j5)*zb(j5,j1))

! JC: I have inserted this sign by hand in order to reproduce
! the correct singular limits 
      A51mpppp=-A51mpppp
      return
      end 
