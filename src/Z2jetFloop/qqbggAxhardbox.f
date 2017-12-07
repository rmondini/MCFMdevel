      function qqbggAxhardbox(p1,p2,p3,p4,p5,p6,za,zb)
C     returns the m^0 coefficient of the 1+,2-,3-,4+,5-,6+ hardbox
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6
      real(dp):: s3
      complex(dp):: qqbggAxhardbox,zab2
C--- begin statement function
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      s3(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)
C--- end statement functions
      qqbggAxhardbox=s(p3,p4)*s3(p1,p2,p4)
     & *((zab2(p2,p1,p4,p3)*zab2(p4,p1,p2,p6))**2
     & -(za(p2,p4)*zb(p3,p6)*s3(p1,p2,p4))**2)
     & /(2._dp*za(p1,p2)*zb(p5,p6)*zab2(p4,p1,p2,p3)**4)
      return
      end
