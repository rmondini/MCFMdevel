      subroutine AGTY_dirgam_assemble_schan(p,i1,i2,i3,F)
!======C.Williams May 2016
!======routine to reproduce eq. 4.5 of AGTY, 
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f' 
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'scet_const.f'
      real(dp), intent(in) ::p(mxpart,4)
      real(dp):: ss,tt,uu
      integer i1,i2,i3
      real(dp):: ddilog,Li3,Li4
      real(dp), intent(out) :: F(3)
!=======basis functions from AGTY
      real(dp) :: AGTYAs,AGTYBs,AGTYCs,AGTYD1s,AGTYE1s
      real(dp) :: AGTYE2s,AGTYF1s,AGTYG1s
      real(dp) :: AGTYX1s,AGTYX2s,AGTYX3s,AGTYG2s,AGTYG3s
      real(dp) :: Bigx,Bigy,Bigs,BigT,BigU,x,y,z
      real(dp) Li2x,Li3x,Li4x,Li2y,Li3y,Li4y,Li2z,Li3z,Li4z,Li4zinv
      integer, parameter:: inite1x1=1,inite0x2=2,initeNf=3
      call dotem(4,p,s)

      !=====Mandelstam invariants 
      ss=s(i1,i2)
      tt=s(i2,i3)
      uu=s(i1,i3)

!----- define various pieces 
      BigX=log(-tt/ss)
      BigY=log(-uu/ss) 
      BigS=log(ss/musq) 
      BigU=log(-uu/musq)
      BigT=log(-tt/musq)
      
      x=-tt/ss
      y=-uu/ss 
      z=-uu/tt

      Li2x=ddilog(x)
      Li3x=Li3(x)
      Li4x=Li4(x)

      Li2y=ddilog(y)
      Li3y=Li3(y)
      Li4y=Li4(y)

      Li2z=ddilog(z)
      Li3z=Li3(z)
      Li4z=Li4(z)

      Li4zinv=Li4(1._dp/z)

!======== now build the functions

!====== this bit is special and does not scale with the LO quark charge. 
      F(initeNf)=AGTYAs(ss,tt,uu,BigX,BigY,Li2x,Li2y,Li3x,Li3y,
     &     Li4x,Li4y,Li4z,Li4zinv)
!------normalization of F(initeNf)
      F(initeNf) = 2._dp*xn*cf*(2*tr*cf-3*tr*ca/4._dp)*F(initeNf)

!======= B pieces (CF**2)
      F(inite0x2)=Cf**2*AGTYBs(ss,tt,uu,BigX,BigY,BigS,Li2x,Li3x,Li4x,
     & Li2y,Li3y,Li4y,Li4z,Li4zinv)
!      write(6,*) 'B = ',F(inite0x2)

!======= C pieces (CA**2)
      F(inite0x2)=F(inite0x2)+CA**2*AGTYCs(ss,tt,uu,BigX,BigY,BigS,Li2x,Li3x,Li4x,
     & Li2y,Li3y,Li4y,Li4z,Li4zinv)
!      write(6,*) 'C = ',F(inite0x2)

!======= D1 pieces (CF*CA)
      F(inite0x2)=F(inite0x2)+CF*CA*AGTYD1s(ss,tt,uu,BigX,BigY,BigS,Li2x,Li3x,Li4x,
     & Li2y,Li3y,Li4y,Li4z,Li4zinv)
!      write(6,*) 'D = ',F(inite0x2)

!======= E1 pieces (NF*CF)
      F(inite0x2)=F(inite0x2)+nf*Cf*AGTYE1s(tt,uu,BigX,BigY,BigS)
!      write(6,*) 'E1 = ',F(inite0x2)

!=======E2 pieces (NF*CA)
      F(inite0x2)=F(inite0x2)+nf*Ca*AGTYE2s(ss,tt,uu,BigX,BigY,BigS,Li2x,Li2y,Li3x,Li3y)
!      write(6,*) 'E2 = ',F(inite0x2)

!=======F1 pieces (NF**2)
      F(inite0x2)=F(inite0x2)+nf**2*AGTYF1s(tt,uu,BigX,BigY,BigS)
!      write(6,*) 'F1 = ',F(inite0x2)
        
!======overall factor outside square brackets in Eq. (4.5)
      F(inite0x2)=F(inite0x2)*2._dp*xn*cf

!************************************************
!==== one-loop squared pieces


!======G1s (CF**2)
      F(inite1x1)=CF**2*AGTYG1s(tt,uu,BigX,BigY)
!======G2s (CF*CA)  
      F(inite1x1)=F(inite1x1)+CF*CA*AGTYG2s(tt,uu,BigX,BigY,BigS)
!======G3s (CA*CA)  
      F(inite1x1)=F(inite1x1)+CA**2*AGTYG3s(tt,uu,BigX,BigY,BigS)
!======X1s (CF*NF)  
      F(inite1x1)=F(inite1x1)+Cf*Nf*AGTYX1s(tt,uu,BigX,BigY,BigS)
!======X1s (CA*NF)  
      F(inite1x1)=F(inite1x1)+CA*Nf*AGTYX2s(tt,uu,BigX,BigY,BigS)
!======X1s (NF**2)  
      F(inite1x1)=F(inite1x1)+Nf**2*AGTYX3s(tt,uu,BigX,BigY,BigS)

!=====overall factor outside square brackets in Eq. (5.2)
      F(inite1x1)=xn*cf*F(inite1x1)
      
!********************************************** 
      return
      end 
