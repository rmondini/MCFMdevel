!---  subtraction routine for H=>b(i1)+bb(i2)+f(i3)
!----- where f= g 
!-----dipoles are all final-final with following structue

      subroutine Hbb_decay_gs(p,i1,i2,i3,msqc)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'incldip.f'
      include 'masses.f'
      include 'breit.f'
      include 'hbbparams.f'
      integer i1,i2,i3,i4
      real(dp):: msqc(maxd),p(mxpart,4)
      real(dp):: sub13_2(4),sub13_2v
      real(dp):: sub23_1(4),sub23_1v
      real(dp) :: msq13_2,msq13_2v
      real(dp) :: msq23_1,msq23_1v
      integer:: j,k,nd
      real(dp), parameter:: tiny=1e-6_dp
      external Hbb_decay_ex
      integer isp_dipmin,isp_dipnpart
      common/isp_dip/isp_dipmin,isp_dipnpart
!$omp threadprivate(/isp_dip/)
      external Hbb_donothing_gvec
      
      msqc(:)=0._dp
               
      ndmax=2

      isp_dipmin=min(i1,i2,i3)
      isp_dipnpart=3

    
!------qg type dipoles
      call dips_HD(1,p,i1,i3,i2,sub13_2,sub13_2v,msq13_2,msq13_2v,
     &     Hbb_decay_ex,Hbb_donothing_gvec)
      call dips_HD(2,p,i2,i3,i1,sub23_1,sub23_1v,msq23_1,msq23_1v,
     &     Hbb_decay_ex,Hbb_donothing_gvec)

!---- build Matrix element

      msqc(1) = two*cf*msq13_2*sub13_2(qq)
      msqc(2) = two*cf*msq23_1*sub23_1(qq)

!      write(6,*) msqc(1),msqc(2)
!      pause
      return
      end

      subroutine Hbb_decay_ex(p,msq)
      implicit none
      include 'types.f'
!-----explicit routine for H=>b(1)+b(2)+g(3) in which
!---- positions of bs and g is fixed for dipole.
      include 'constants.f'
      include 'mxpart.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'hbbparams.f'
      include 'sprods_com.f'
      real(dp) :: p(mxpart,4),msq,s56
      
      s56=2._dp*(p(1,4)*p(2,4)-p(1,3)*p(2,3)-p(1,2)*p(2,2)-p(1,1)*p(2,1))
!----- lo
      msq=xn*gwsq*mb_eff**2/(4._dp*wmass**2)*2._dp*s56
!      write(6,*) 'ex = ',msq
!      pause
      return
      end


      subroutine Hbb_donothing_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      integer:: j,k,in
      real(dp):: msq,p(mxpart,4),n(4)
     
   
      msq=0._dp
      return
      end
