      subroutine namp_qqbggga(i1,i2,i3,i4,i5,za,zb,zab,zba,amp)
!--- CW June 16
!---  amplitudes for
!=====q(i1) + qb(i2) + g(i3) + g(i4) + gamma(i5) with
!-------- g(i4) contracted with vector n 
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer i1,i2,i3,i4,i5
      complex(dp):: zab(mxpart,mxpart),zba(mxpart,mxpart),amp(2,2,2)
      complex(dp) :: n_qqbggga_pp,n_qqbggga_pm,n_qqbggga_mm,n_qqbggga_mp
 
      amp(:,:,:)=czip
      
!----helicity array i1,i3,i5 
      amp(1,2,1)=n_qqbggga_pm(i1,i2,i3,i4,i5,za,zb,zab)
      amp(1,2,2)=n_qqbggga_pp(i1,i2,i3,i4,i5,za,zb,zab)
      amp(1,1,1)=n_qqbggga_mm(i1,i2,i3,i4,i5,za,zb,zab)
      amp(1,1,2)=n_qqbggga_mp(i1,i2,i3,i4,i5,za,zb,zab)

      
      amp(2,2,1)=-n_qqbggga_mp(i1,i2,i3,i4,i5,zb,za,zba)
      amp(2,1,2)=-n_qqbggga_pm(i1,i2,i3,i4,i5,zb,za,zba)
      amp(2,1,1)=-n_qqbggga_pp(i1,i2,i3,i4,i5,zb,za,zba)
      amp(2,2,2)=-n_qqbggga_mm(i1,i2,i3,i4,i5,zb,za,zba)
          
      return
     
      end

      function n_qqbggga_mp(i1,i2,i3,i4,i5,za,zb,zab)
      implicit none
      include 'types.f'
      complex(dp) :: n_qqbggga_mp
      include 'mxpart.f'
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer i1,i2,i3,i4,i5
      complex(dp):: zab(mxpart,mxpart)

      n_qqbggga_mp=
     &  (za(i1,i2)*za(i1,i4)*zb(i4,i1)*
     -     (-(za(i1,i3)*zab(i3,i2)*zb(i3,i2)) + 
     -       (-2*za(i3,i4)*zab(i1,i2) + za(i1,i4)*zab(i3,i2))*zb(i4,i2))
     -     - za(i1,i5)*(2*za(i3,i4)*
     -        (-(za(i1,i3)*zab(i1,i1)) + za(i3,i4)*zab(i1,i4))*zb(i4,i3)*
     -        zb(i5,i2) + za(i1,i4)*zb(i4,i1)*
     -        (2*za(i3,i4)*zab(i1,i5)*zb(i4,i2) + 
     -          za(i1,i3)*(2*zab(i3,i3)*zb(i5,i2) - 
     -             zab(i3,i2)*zb(i5,i3)) + za(i1,i4)*zab(i3,i2)*zb(i5,i4))
     -       ))/
     -  (2.*za(i1,i4)*za(i1,i5)*za(i2,i5)*za(i3,i4)*zb(i3,i2)*zb(i4,i1)*
     -    zb(i4,i3))

      n_qqbggga_mp=n_qqbggga_mp/rt2

      return
      end
      
      function n_qqbggga_mm(i1,i2,i3,i4,i5,za,zb,zab)
      implicit none
      include 'types.f'
      complex(dp) :: n_qqbggga_mm
      include 'mxpart.f'
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer i1,i2,i3,i4,i5
      complex(dp):: zab(mxpart,mxpart)

      n_qqbggga_mm=     (za(i3,i4)*zb(i4,i2)*(-2*zab(i1,i2)*zb(i2,i1) + 
     -       2*zab(i5,i2)*zb(i5,i2)) + 
     -    zab(i3,i2)*(-(za(i1,i3)*zb(i2,i1)*zb(i3,i2)) + 
     -       za(i1,i4)*zb(i2,i1)*zb(i4,i2) + 
     -       (-(za(i3,i5)*zb(i3,i2)) + za(i4,i5)*zb(i4,i2))*zb(i5,i2)))/
     -  (2.*za(i3,i4)*zb(i3,i2)*zb(i4,i3)*zb(i5,i1)*zb(i5,i2))

      n_qqbggga_mm=n_qqbggga_mm/rt2

      return
      end
      
      function n_qqbggga_pp(i1,i2,i3,i4,i5,za,zb,zab)
      implicit none
      include 'types.f'
      complex(dp) :: n_qqbggga_pp
      include 'mxpart.f'
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer i1,i2,i3,i4,i5
      complex(dp):: zab(mxpart,mxpart)

      n_qqbggga_pp=(2._dp*za(i1,i2)**2*zab(i3,i3)*zb(i3,i2) - 
     -    za(i1,i2)*(za(i1,i3)*zab(i2,i3)*zb(i3,i2) - 
     -       za(i1,i4)*zab(i2,i3)*zb(i4,i2) + 
     -       2._dp*za(i2,i4)*zab(i1,i2)*zb(i4,i3) + 
     -       2._dp*za(i3,i4)*zab(i1,i3)*zb(i4,i3) + 
     -       2._dp*za(i1,i5)*zab(i3,i3)*zb(i5,i3)) + 
     -    za(i1,i5)*(-2._dp*za(i2,i4)*zab(i1,i5)*zb(i4,i3) + 
     -       zab(i2,i3)*(za(i1,i3)*zb(i5,i3) - za(i1,i4)*zb(i5,i4))))/
     -  (2._dp*za(i1,i5)*za(i2,i3)*za(i2,i5)*za(i3,i4)*zb(i4,i3))

      n_qqbggga_pp=n_qqbggga_pp/rt2

      return
      end
      
      function n_qqbggga_pm(i1,i2,i3,i4,i5,za,zb,zab)
      implicit none
      include 'types.f'
      complex(dp) :: n_qqbggga_pm
      include 'mxpart.f'
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer i1,i2,i3,i4,i5
      complex(dp):: zab(mxpart,mxpart)

      n_qqbggga_pm=
     &     (-2*za(i1,i4)**2*za(i2,i3)*zab(i3,i3)*zb(i2,i1)*zb(i3,i2)*
     -     zb(i4,i1) + 2*za(i3,i4)**2*
     -     (za(i1,i5)*zab(i1,i1) + za(i4,i5)*zab(i1,i4))*zb(i3,i2)*
     -     zb(i4,i3)*zb(i5,i1) + 
     -    za(i1,i4)*zb(i4,i1)*
     -     (za(i2,i3)*(za(i1,i3)*zab(i4,i3)*zb(i2,i1)*zb(i3,i2) - 
     -          za(i1,i4)*zab(i4,i3)*zb(i2,i1)*zb(i4,i2) + 
     -          (za(i3,i5)*zab(i4,i3)*zb(i3,i2) - 
     -             za(i4,i5)*(2*zab(i3,i3)*zb(i3,i2) + 
     -                zab(i4,i3)*zb(i4,i2)))*zb(i5,i2)) + 
     -       2*za(i3,i4)*zb(i4,i3)*
     -        (za(i2,i4)*(-(zab(i1,i2)*zb(i2,i1)) + 
     -             zab(i5,i2)*zb(i5,i2)) + 
     -          za(i3,i4)*(-(zab(i1,i3)*zb(i2,i1)) + 
     -             zab(i5,i3)*zb(i5,i2)))))/
     -  (2.*za(i1,i4)*za(i2,i3)*za(i3,i4)**2*zb(i4,i1)*zb(i4,i3)*
     -    zb(i5,i1)*zb(i5,i2))

      n_qqbggga_pm=n_qqbggga_pm/rt2

      return
      end
