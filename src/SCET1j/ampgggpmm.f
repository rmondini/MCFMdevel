      function ampgggpmm(p1,p2,p3,p5,p6,alpa,alpb,alpc,za,zb)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
!     This routine calculates the amplitude for a right-handed quark
!     0--> g^+(1)+g^-(2)+g^-(3)+l^+(5)+lb^-(6)
!     according to Eq.(36) of 1309.3245v3
      integer:: p1,p2,p3,p5,p6
      complex(dp)::ampgggpmm
      complex(dp)::alpa,alpb,alpc
      ampgggpmm=za(p2,p3)/(za(p1,p2)*za(p1,p3)*zb(p2,p3))*(
     & +alpa*za(p2,p5)*za(p3,p5)*zb(p5,p6)
     & -alpb*za(p2,p3)*za(p2,p5)*zb(p2,p6)
     & -alpc*za(p2,p3)*za(p3,p5)*zb(p3,p6))/rt2
     
      return
      end
