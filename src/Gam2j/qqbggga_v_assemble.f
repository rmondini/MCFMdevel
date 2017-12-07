!---- CW June 16, routine to build the
!-----q(i1), qb(i2) g(i3) g(i4), gamma(i5), interfered with the LO
      
      subroutine qqbggga_v_assemble(p,i1,i2,i3,i4,i5,msqv)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'qcdouble.f'
      real(dp):: p(mxaprt,4),msqv
      integer i1,i2,i3,i4,i5
      real(dp) :: faclo,fac_virt,fac
      complex(dp) :: amp_lo(2,2,2,2),amp_virt(2,2,2,2)
      complex(dp) :: amp_nf(2,2,2,2)
      
!      fac_virt=rt2*sqrt(esq)
      
      fac=2._dp*esq*gsq

!------call lo

      

      
