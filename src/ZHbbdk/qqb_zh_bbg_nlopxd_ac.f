      subroutine qqb_zh_bbg_nlopxd_ac(P,msq)
      implicit none
      include 'types.f'
!------ this is the double real contribution for NLO X NLO ZH production 
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'hbbparams.f'
      include 'hdecaymode.f'
      integer:: j,k
      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: radiLL_Hbbg
      real(dp):: qqbZHgL,qbqZHgL,qgZHqL,gqZHqL,gqbZHqbL,qbgZHqbL
      real(dp):: qqbZHgR,qbqZHgR,qgZHqR,gqZHqR,gqbZHqbR,qbgZHqbR
      integer ig 
!--- ig is position of emitted Initial State gluon in phase-space
      ig=8 
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call dotem(ig,p,s)

      qqbZHgL=aveqq*radiLL_Hbbg(p,1,2,ig,5,6,7,3,4)
      qqbZHgR=aveqq*radiLL_Hbbg(p,1,2,ig,5,6,7,4,3)
      qbqZHgL=aveqq*radiLL_Hbbg(p,2,1,ig,5,6,7,3,4)
      qbqZHgR=aveqq*radiLL_Hbbg(p,2,1,ig,5,6,7,4,3)

      qgZHqL=-radiLL_Hbbg(p,1,ig,2,5,6,7,3,4)*aveqg
      qgZHqR=-radiLL_Hbbg(p,1,ig,2,5,6,7,4,3)*aveqg
      gqZHqL=-radiLL_Hbbg(p,2,ig,1,5,6,7,3,4)*aveqg
      gqZHqR=-radiLL_Hbbg(p,2,ig,1,5,6,7,4,3)*aveqg

      gqbZHqbL=-radiLL_Hbbg(p,ig,2,1,5,6,7,3,4)*aveqg
      gqbZHqbR=-radiLL_Hbbg(p,ig,2,1,5,6,7,4,3)*aveqg

      qbgZHqbL=-radiLL_Hbbg(p,ig,1,2,5,6,7,3,4)*aveqg
      qbgZHqbR=-radiLL_Hbbg(p,ig,1,2,5,6,7,4,3)*aveqg

      do j=-nf,nf
      do k=-nf,nf

      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 40
      if ((j == 0) .and. (k == 0)) goto 40

      if     ((j > 0) .and. (k < 0)) then
          msq(j,k)=((L(j)*l1)**2+(R(j)*r1)**2)*qqbZHgL
     &            +((R(j)*l1)**2+(L(j)*r1)**2)*qqbZHgR
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=((L(k)*l1)**2+(R(k)*r1)**2)*qbqZHgL
     &            +((R(k)*l1)**2+(L(k)*r1)**2)*qbqZHgR
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=((L(j)*l1)**2+(R(j)*r1)**2)*qgZHqL
     &            +((R(j)*l1)**2+(L(j)*r1)**2)*qgZHqR
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=((L(-j)*l1)**2+(R(-j)*r1)**2)*qbgZHqbL
     &            +((R(-j)*l1)**2+(L(-j)*r1)**2)*qbgZHqbR
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=((L(k)*l1)**2+(R(k)*r1)**2)*gqZHqL
     &            +((R(k)*l1)**2+(L(k)*r1)**2)*gqZHqR
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=((L(-k)*l1)**2+(R(-k)*r1)**2)*gqbZHqbL
     &            +((R(-k)*l1)**2+(L(-k)*r1)**2)*gqbZHqbR
      endif
 40   continue
      enddo
      enddo

      
      return
      end


      function radiLL_Hbbg(p,j1,j2,j3,j4,j5,j8,j6,j7)
      implicit none
      include 'types.f'
      real(dp):: radiLL_Hbbg
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'hdecaymode.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      integer:: j1,j2,j3,j4,j5,j6,j7,j8
      real(dp):: s458,s12,s13,s23,s123,prop
      real(dp):: fac,hdecay,msqhbb,msqhtautau,msqhgamgam
      real(dp):: p(mxpart,4)
      
      s12=s(j1,j2)
      s13=s(j1,j3)
      s23=s(j2,j3)
      s123=s12+s13+s23
c---calculate the 2 Z propagators
      prop=         ((s123-zmass**2)**2+(zmass*zwidth)**2)
      prop=prop*((s(j6,j7)-zmass**2)**2+(zmass*zwidth)**2)
      fac=four*V*gsq*esq**2
      fac=fac*gwsq*wmass**2/(one-xw)**2/prop

!---- Higgs decay => Hbbg
      call Hbbg_decay(p,j4,j5,j8,hdecay)
      s458=s(j4,j5)+s(j4,j8)+s(j5,j8)
      
      hdecay=hdecay/((s458-hmass**2)**2+(hmass*hwidth)**2)
      fac=fac*hdecay

      radiLL_Hbbg=s12/s13/s23
     & *(2._dp*s(j1,j7)*s(j2,j6)+s(j1,j7)*s(j3,j6)+s(j2,j6)*s(j3,j7))
     & +(s(j1,j7)*s(j2,j6)+s(j2,j6)*s(j3,j7)-s(j1,j6)*s(j1,j7))/s13
     & +(s(j1,j7)*s(j2,j6)+s(j1,j7)*s(j3,j6)-s(j2,j6)*s(j2,j7))/s23

      radiLL_Hbbg=fac*radiLL_Hbbg
      return
      end
