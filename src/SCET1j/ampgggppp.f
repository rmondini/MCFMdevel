      function ampgggppp(p1,p2,p3,p5,p6,beta,betb,betc,za,zb)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
!     This routine calculates the amplitude for a right-handed quark
!     0--> g^+(1)+g^+(2)+g^+(3)+l^+(5)+lb^-(6)
!     according to Eq.(37) of 1309.3245v3
      integer:: p1,p2,p3,p5,p6
      complex(dp)::ampgggppp
      complex(dp)::beta,betb,betc
      ampgggppp=(
     & +beta*zb(p1,p3)*za(p1,p5)*zb(p1,p6)/(za(p1,p2)*za(p2,p3))
     & +betb*zb(p2,p3)*za(p2,p5)*zb(p2,p6)/(za(p1,p2)*za(p1,p3))
     & +betc*zb(p2,p3)*za(p2,p5)*zb(p1,p6)/(za(p1,p2)*za(p2,p3)))/rt2
      return
      end
