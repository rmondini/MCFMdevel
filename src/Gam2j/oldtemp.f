
      subroutine ga_qqbQQb(i1,i2,i4,i5,za,zb,amp_a,amp_b)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer i1,i2,i3,i4,i5
      complex(dp) :: amp_a(2,2,2),amp_b(2,2,2)
!---- amp_A couples to i1,i2
!---- amp_B couples to i4,i5
      complex(dp) ga_qqbQQb_amp
      
      amp_a(1,1,2)=ga_qqbQQb_amp(i1,i2,i4,i5,3,za,zb)
      amp_b(1,1,2)=ga_qqbQQb_amp(i4,i5,i1,i2,3,za,zb)
      
      return
      end

    
      function ga_qqbQQb_amp(i1,i2,i3,i4,i5,za,zb)
!-------amplitude for q(i1)qb(i2),Qb(i3)QB(i4),ga(i5), ga couples to i1,i2
      implicit none
      include 'types.f'
      complex(dp):: ga_qqbQQb_amp
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      
      ga_qqbQQb_amp=-(za(i1,i3)**2/(za(i1,i5)*za(i2,i5)*za(i3,i4)))

      return
      end
