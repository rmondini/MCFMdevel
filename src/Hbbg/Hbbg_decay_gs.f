!---  subtraction routine for H=>b(i1)+bb(i2)+f(i3)+f(i4)
!----- where f= gg or QQb
!-----dipoles are all final-final with following structue

!-----q(i1) g(i3/i4)  qb(i2)  qg type
!----
!---  (i1,i3) i2
!---  (i1,i4) i2
!---- (i2,i3) i1
!---- (i2,i4) i1
      
!-----g(i3) g(i4) q(i1) gg (or QQb) type

!----- (i7,i8) i1 
!----- (i8,i7) i1
!----- (i7,i8) i2
!----- (i8,i7) i2


!---- for a total of 8 dipoles.

      subroutine Hbbg_decay_gs(p,i1,i2,i3,i4,msqc)
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
      real(dp):: sub13_2(4),sub13_2v,sub14_2(4),sub14_2v,
     &  sub24_1(4),sub24_1v,sub23_1(4),sub23_1v,   
     &  sub34_1(4),sub34_1v,sub43_1(4),sub43_1v
      real(dp) :: msq13_2,msq13_2v,msq14_2,msq14_2v,msq34_1,msq34_1v
      real(dp) :: msq24_1,msq24_1v,msq23_1,msq23_1v,msq43_1,msq43_1v
      real(dp) :: sub13_4(4),sub13_4v,msq13_4,msq13_4v
      real(dp) :: sub14_3(4),sub14_3v,msq14_3,msq14_3v
      real(dp) :: sub24_3(4),sub24_3v,msq24_3,msq24_3v
      real(dp) :: sub23_4(4),sub23_4v,msq23_4,msq23_4v
      real(dp) :: sub34_2(4),sub34_2v,msq34_2,msq34_2v
      real(dp) :: sub43_2(4),sub43_2v,msq43_2,msq43_2v
      integer:: j,k,nd
      real(dp), parameter:: tiny=1e-6_dp
      external Hbbg_decay_ex,Hbbg_decay_gvec_ex
      integer isp_dipmin,isp_dipnpart
      common/isp_dip/isp_dipmin,isp_dipnpart
!$omp threadprivate(/isp_dip/)
      external Hbbg_donothing_gvec
      logical nfonly
      common/nfonly/nfonly
!$omp threadprivate(/nfonly/)

      msqc(:)=0._dp
               
      ndmax=12

      isp_dipmin=min(i1,i2,i3,i4)
      isp_dipnpart=3

    
!------qg type dipoles
      call dips_HD(1,p,i1,i3,i2,sub13_2,sub13_2v,msq13_2,msq13_2v,
     &     Hbbg_decay_ex,Hbbg_donothing_gvec)
      call dips_HD(7,p,i1,i3,i4,sub13_4,sub13_4v,msq13_4,msq13_4v,
     &     Hbbg_decay_ex,Hbbg_donothing_gvec)

      call dips_HD(2,p,i1,i4,i2,sub14_2,sub14_2v,msq14_2,msq14_2v,
     & Hbbg_decay_ex,Hbbg_donothing_gvec)
      call dips_HD(8,p,i1,i4,i3,sub14_3,sub14_3v,msq14_3,msq14_3v,
     & Hbbg_decay_ex,Hbbg_donothing_gvec)

      call dips_HD(3,p,i2,i3,i1,sub23_1,sub23_1v,msq23_1,msq23_1v,
     & Hbbg_decay_ex,Hbbg_donothing_gvec)
      call dips_HD(9,p,i2,i3,i4,sub23_4,sub23_4v,msq23_4,msq23_4v,
     & Hbbg_decay_ex,Hbbg_donothing_gvec)

      call dips_HD(4,p,i2,i4,i1,sub24_1,sub24_1v,msq24_1,msq24_1v,
     & Hbbg_decay_ex,Hbbg_donothing_gvec)
      call dips_HD(10,p,i2,i4,i3,sub24_3,sub24_3v,msq24_3,msq24_3v,
     & Hbbg_decay_ex,Hbbg_donothing_gvec)
      
!---- QQ and gg types dipoles
      call dips_HD(5,p,i3,i4,i1,sub34_1,sub34_1v,msq34_1,msq34_1v,
     &     Hbbg_decay_ex,Hbbg_decay_gvec_ex)
      call dips_HD(6,p,i4,i3,i1,sub43_1,sub43_1v,msq43_1,msq43_1v,
     &     Hbbg_decay_ex,Hbbg_decay_gvec_ex)

!---  
      call dips_HD(11,p,i3,i4,i2,sub34_2,sub34_2v,msq34_2,msq34_2v,
     & Hbbg_decay_ex,Hbbg_decay_gvec_ex)
      call dips_HD(12,p,i4,i3,i2,sub43_2,sub43_2v,msq43_2,msq43_2v,
     &     Hbbg_decay_ex,Hbbg_decay_gvec_ex)

!---- build Matrix element

!---- work out pre-factors!
      
      msqc(1) = -half*msq13_2*sub13_2(qq)*(ca-2._dp*cf)
      msqc(2) = -half*msq14_2*sub14_2(qq)*(ca-2._dp*cf)
      msqc(3) = -half*msq23_1*sub23_1(qq)*(ca-2._dp*cf)
      msqc(4) = -half*msq24_1*sub24_1(qq)*(ca-2._dp*cf)
      msqc(7) = half*msq13_4*sub13_4(qq)*(ca)
      msqc(8) = half*msq14_3*sub14_3(qq)*(ca)
      msqc(9) = half*msq23_4*sub23_4(qq)*(ca)
      msqc(10) =half*msq24_3*sub24_3(qq)*(ca)

      
      msqc(5)= (ca)*(msq34_1*sub34_1(gg)+msq34_1v*sub34_1v)/four
      msqc(6)= (ca)*(msq43_1*sub43_1(gg)+msq43_1v*sub43_1v)/four
      msqc(11)= (ca)*(msq34_2*sub34_2(gg)+msq34_2v*sub34_2v)/four
      msqc(12)= (ca)*(msq43_2*sub43_2(gg)+msq43_2v*sub43_2v)/four

!------debug-just nf pieces
!      if(nfonly) then 
!         msqc(:)=zip
!      endif

!-----no nf pieces
!      return 

!---- and QQb contributions
      msqc(5) =  msqc(5)+tr*nf*(msq34_1*sub34_1(gq)-msq34_1v*sub34_1v)*half
      msqc(6) =  msqc(6)+tr*nf*(msq43_1*sub43_1(gq)-msq43_1v*sub43_1v)*half
      msqc(11) =  msqc(11)+tr*nf*(msq34_2*sub34_2(gq)-msq34_2v*sub34_2v)*half
      msqc(12) =  msqc(12)+tr*nf*(msq43_2*sub43_2(gq)-msq43_2v*sub43_2v)*half

      
      
      return
      end


      subroutine Hbbg_decay_ex(p,msq)
      implicit none
      include 'types.f'
!-----explicit routine for H=>b(1)+b(2)+g(3) in which
!---- positions of bs and g is fixed for dipole.
      include 'constants.f'
      include 'mxpart.f' 
      real(dp) :: p(mxpart,4),msq

      call Hbbg_decay(p,1,2,3,msq)
      return
      end


      subroutine Hbbg_donothing_gvec(p,n,in,msq)
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

