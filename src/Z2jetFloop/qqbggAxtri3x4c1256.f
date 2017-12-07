      subroutine qqbggAxtri3x4x1256(p1,p2,p3,p4,p5,p6,za,zb,coeff)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer p1,p2,p3,p4,p5,p6
      real(dp):: s3,mtsq,s12,s123,s124,s1234,s34
      complex(dp):: coeff(2,2),iza,izb,zab2,zba2,izab2,izba2,
     & qqbggAxhardbox

C--- begin statement function
      s3(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)
C--- end statement functions

      coeff(1,1) =czip
      coeff(2,2) =czip
      coeff(1,2)=
     & -qqbggAxhardbox(p1,p2,p3,p4,p5,p6,za,zb)/s3(p1,p2,p4)
     & -qqbggAxhardbox(p2,p1,p4,p3,p6,p5,zb,za)/s3(p1,p2,p3)
      coeff(2,1)=
     & +qqbggAxhardbox(p2,p1,p3,p4,p6,p5,zb,za)/s3(p1,p2,p4)
     & +qqbggAxhardbox(p1,p2,p4,p3,p5,p6,za,zb)/s3(p1,p2,p3)

      return
      end




   
