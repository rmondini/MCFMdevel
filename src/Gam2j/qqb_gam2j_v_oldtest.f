
      subroutine qqb_gam2j_v(p,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'epinv.f'
      real(dp) :: msq(-nf:nf,-nf:nf),p(mxpart,4)
      complex(dp):: amp_qqbggga_slc_gaMHV,test
      complex(dp):: amp_qqbggga_slc_gMHVadj
      complex(dp):: amp_qqbggga_lc_gMHVadj
      complex(dp):: amp_qqbggga_lc_gMHV,amp_qqbggga_lc_gaMHV
      complex(dp):: amp_qqbggga_slc_gMHV
      complex(dp):: amp_qqbggga_del34_gaMHV
      complex(dp):: amp_qqbggga_del34_gMHV,amp_qqbQQbga_del1_lc
      complex(dp)::amp_del1(2,2,2),amp_del2(2,2,2)
      complex(dp):: amp_34(2,2,2,2),amp_43(2,2,2,2),amp_del(2,2,2,2)
!=====setup KC
!-----DEBUG
      epinv=0._dp
      musq=1._dp
      include 'kinpoint_qqbggga_kc.f'
!      include 'kinpoint_qqbQQga_kc.f'

!-----END DEBUG
!      call writeout(p) 
      call spinorz(5,p,za,zb)
!      call amp_qqbQQbga_v_assemble(1,2,3,4,5,za,zb,amp_del1,amp_del2)
      call amp_qqbggga_v_del34_assemble(1,2,3,4,5,za,zb,amp_34
     &     ,amp_43,amp_del)
!      test=amp_qqbQQbga_del1_lc(1,2,3,4,5,za,zb)
      stop
      return
      end
