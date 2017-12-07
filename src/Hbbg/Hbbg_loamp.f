!-----amplitudes for H->b(i1,h1)+bbar(i2,h1)+g(i3,h2)
      subroutine Hbbg_loamp(i1,i2,i3,za,zb,amps)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f' 
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer i1,i2,i3
      complex(dp) amps(2,2),s123
      s123=s(i1,i2)+s(i2,i3)+s(i1,i3)
!---- ---
      amps(1,1)=s123/(zb(i1,i3)*zb(i2,i3))
!     --+
      amps(1,2)=za(i1,i2)**2/(za(i1,i3)*za(i2,i3))
!     +++
      amps(2,2)=-s123/(za(i1,i3)*za(i2,i3))
!     ++-
      amps(2,1)=-zb(i1,i2)**2/(zb(i1,i3)*zb(i2,i3))
      return
      end
      
      
