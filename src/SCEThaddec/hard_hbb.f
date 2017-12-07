      subroutine hard_hbb(Qsq,musq,hard)
      implicit none
!    Hard function for H=> bb (Qsq > 0) units of as/2/pi *LO 
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'
      real(dp),intent(in)::Qsq,musq
      real(dp),intent(out)::hard(2)
      complex(dp)::coeff(2)
      real(dp) :: tag,t1,t2
      logical nfonly
      common/nfonly/nfonly
!$omp threadprivate(/nfonly/)

      tag=one
      t1=one
      t2=one
      
!====== nlo coefficient
      hard(1)=-2*CF + (7*CF*pisq)/6. - CF*Log(musq/Qsq)**2

!===== nnlo coefficient  !--NF pieces
      hard(2)=  (CF*NF*(400 + pisq*(171 - 336*t2) + 18*(-1 + 3*t2)*zeta3 + 
     -      3*Log(musq/Qsq)*(-68 + pisq*(9 - 33*t2) + 180*t2 + 
     -     6*Log(musq/Qsq)*(-9 + 19*t2 + 2*t2*Log(musq/Qsq)))))/324.
      
!---- uncomment out this bit for nf only
!      hard(1)=zip
!      return 
!==== nnlo coefficient  !--- CA*CF pieces
      hard(2)=hard(2)+ (CA*CF*(9*Pi**4*(93 - 125*t2) + 30*pisq*(-681 + 1382*t2) - 
     -      20*(934 + 9*(-335 + 33*t2)*zeta3) + 
     -      30*(3*pisq*(-69 + 157*t2) + 4*(361 - 603*t2 + 702*zeta3))*Log(musq/Qsq) + 
     -     180*(99 + (-233 + 6*pisq)*t2)*Log(musq/Qsq)**2 - 3960*t2*Log(musq/Qsq)**3))/6480.
      
!=====nnlo coefficient ! CF^2 pieces

      hard(2)=hard(2)+ (CF**2*(60*pisq*(-11 + 40*t1 - 27*t2 - 14*tag) + Pi**4*(-88 + 245*t2 + 245*tag) + 
     -      720*(3 + tag - 15*zeta3) + 
     -      60*Log(musq/Qsq)*(-36 + 36*t1 + 3*pisq*(4 + 19*t1 - 19*t2) - 144*zeta3 + 
     -         Log(musq/Qsq)*(27 - 42*t1 + (27 - 19*pisq)*t2 + 12*tag + 5*pisq*tag + 
     -            3*Log(musq/Qsq)*(-6*t1 + 6*t2 + (t2 + tag)*Log(musq/Qsq))))))/720.
    
!---- No nf pieces- uncomment below
!      hard(2)=hard(2)- (CF*NF*(400 - 165*pisq + 36*zeta3 + 
!     -      12*Log(musq/Qsq)*(28 - 6*pisq + 3*Log(musq/Qsq)*(5 + Log(musq/Qsq)))))/324._dp

!      hard(2)=hard(2)/2._dp
      return
      
      return
      end
