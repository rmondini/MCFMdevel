      subroutine lumxmsq_h_bbg(p,order,xmsq)
      implicit none
      include 'types.f'
c---- Matrix element for H -> bbg decay NNLO production
c---- takes in 
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'facscale.f'
      include 'ewcharge.f'
      include 'scet_const.f'
      include 'hbbparams.f'
      integer:: order
      real(dp),intent(out):: xmsq
      real(dp):: p(mxpart,4),
     & soft1(-1:1),soft2(-1:3),hard(2),
     & bit,lobit,
     & assemble_dec_3j
      real(dp) :: xvar,yvar,zvar,y12,y13,y23
      real(dp) :: jeta1(-1:1),jeta2(-1:3),jetb1(-1:1),jetb2(-1:3),jetc1(-1:1),jetc2(-1:3)

c---------------------------------
!---- hard-code order and 3 massless momenta that sum to (hmass,0,0,0)

      order=2

      p(1,4)=hmass*0.32_dp
      p(1,1)=hmass*0.12_dp
      p(1,2)=hmass*(-0.02_dp)
      p(1,3)=hmass*(-0.295973_dp)

      p(2,4)=hmass*0.46_dp
      p(2,1)=hmass*0.01_dp
      p(2,2)=hmass*0.1165978_dp
      p(2,3)=hmass*0.444865_dp

      p(3,4)=hmass*0.22_dp
      p(3,1)=hmass*(-0.13_dp)
      p(3,2)=hmass*(-0.0965978_dp)
      p(3,3)=hmass*(-0.148892_dp)
c---------------------------------

!---- compute invariants
      call dotem(3,p,s)

c---------------------------------
      write(*,*) 'p(1,4)=',p(1,4)
      write(*,*) 'p(1,1)=',p(1,1)
      write(*,*) 'p(1,2)=',p(1,2)
      write(*,*) 'p(1,3)=',p(1,3)
      write(*,*)
      write(*,*) 'p(2,4)=',p(2,4)
      write(*,*) 'p(2,1)=',p(2,1)
      write(*,*) 'p(2,2)=',p(2,2)
      write(*,*) 'p(2,3)=',p(2,3)
      write(*,*)
      write(*,*) 'p(3,4)=',p(3,4)
      write(*,*) 'p(3,1)=',p(3,1)
      write(*,*) 'p(3,2)=',p(3,2)
      write(*,*) 'p(3,3)=',p(3,3)
c---------------------------------

!---- mandelstam invariants for H -> b(p1) bbar(p2) g(p3)
      xvar=s(1,2)/hmass**2
      yvar=s(1,3)/hmass**2
      zvar=s(2,3)/hmass**2

c---------------------------------
      write(*,*) 'xvar=',xvar
      write(*,*) 'yvar=',yvar
      write(*,*) 'zvar=',zvar
c---------------------------------

!---- LO H->bbg msq
      lobit=32._dp*pi**2*gwsq*mb_eff**2/(four*wmass**2)*(xvar**2+one)/(yvar*zvar)
      lobit=ason2pi*lobit*xn*cf

c---------------------------------
      write(*,*) 'lobit=',lobit
c---------------------------------

!---- soft function
      y12=s(1,2)/p(1,4)/p(2,4)/four
      y13=s(1,3)/p(1,4)/p(3,4)/four
      y23=s(2,3)/p(2,4)/p(3,4)/four

!      y13=0.12793846717065338_dp
!      y23=0.59619034077312028_dp
!      y12=0.88815670191868712_dp

!      y13=0.28667674908342355     
!      y23=0.34652200970874686
!      y12=0.0776031649478127750

      call soft_dec_3j(order,y12,y13,y23,soft1,soft2)

c---------------------------------
      write(*,*) 'y12=',y12
      write(*,*) 'y13=',y13
      write(*,*) 'y23=',y23
      write(*,*) 'soft1=',soft1
      write(*,*) 'soft2=',soft2
c---------------------------------

c---------------------------------
      musq=15000._dp
      write(*,*) 'musq=',musq
c---------------------------------

!---- hard function 
      call hard_hbbg(order,yvar,zvar,musq,hard)

c---------------------------------
      write(*,*) 'hard(1)=',hard(1)
      write(*,*) 'hard(2)=',hard(2)
c---------------------------------

!---- jet functions
      call jetq(order,two*p(1,4),jeta1,jeta2)
      call jetq(order,two*p(2,4),jetb1,jetb2)
      call jetg(order,two*p(3,4),jetc1,jetc2)
      
c---------------------------------
      write(*,*) 'jeta1=',jeta1
      write(*,*) 'jeta2=',jeta2
      write(*,*) 'jetb1=',jetb1
      write(*,*) 'jetb2=',jetb2
      write(*,*) 'jetc1=',jetc1
      write(*,*) 'jetc2=',jetc2
c---------------------------------

      xmsq=zip

      bit=assemble_dec_3j(order,
     & jeta1,jetb1,jetc1,jeta2,jetb2,jetc2,
     & soft1,soft2,hard)

      bit=bit*lobit
      xmsq=xmsq+bit

c---------------------------------
      write(*,*) 'xmsq=',xmsq
      stop
c---------------------------------

      return
      end

