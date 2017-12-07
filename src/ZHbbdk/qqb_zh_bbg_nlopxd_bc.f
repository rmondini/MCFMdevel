!---  CW Sept 17
!---  this is a piece of the NLO x NLO  qqb => Z (H=>bb) process
!---  It corresponds to : The LO qqb => ZH+g  x the below cut H=>bb (i.e. real virtual,)

      
      subroutine qqb_zh_bbg_nlopxd_bc(p,msq)
      implicit none 
      include 'types.f'
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  H  + Z +g(p7)
c                           |    |
c                           |    --> fermion(p3)+antifermion(p4)
c                           |
c                           ---> b(p5)+b(p6)
c   for the moment --- radiation only from initial line
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'hbbparams.f'
      include 'hdecaymode.f'
      integer:: j,k
      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: radiLL_bc
      real(dp):: qqbZHgL,qbqZHgL,qgZHqL,gqZHqL,gqbZHqbL,qbgZHqbL
      real(dp):: qqbZHgR,qbqZHgR,qgZHqR,gqZHqR,gqbZHqbR,qbgZHqbR


      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call dotem(7,p,s)

      qqbZHgL=aveqq*radiLL_bc(1,2,7,5,6,3,4,p(5,4),p(6,4))
      qqbZHgR=aveqq*radiLL_bc(1,2,7,5,6,4,3,p(5,4),p(6,4))
      qbqZHgL=aveqq*radiLL_bc(2,1,7,5,6,3,4,p(5,4),p(6,4))
      qbqZHgR=aveqq*radiLL_bc(2,1,7,5,6,4,3,p(5,4),p(6,4))

      qgZHqL=-radiLL_bc(1,7,2,5,6,3,4,p(5,4),p(6,4))*aveqg
      qgZHqR=-radiLL_bc(1,7,2,5,6,4,3,p(5,4),p(6,4))*aveqg
      gqZHqL=-radiLL_bc(2,7,1,5,6,3,4,p(5,4),p(6,4))*aveqg
      gqZHqR=-radiLL_bc(2,7,1,5,6,4,3,p(5,4),p(6,4))*aveqg

      gqbZHqbL=-radiLL_bc(7,2,1,5,6,3,4,p(5,4),p(6,4))*aveqg
      gqbZHqbR=-radiLL_bc(7,2,1,5,6,4,3,p(5,4),p(6,4))*aveqg

      qbgZHqbL=-radiLL_bc(7,1,2,5,6,3,4,p(5,4),p(6,4))*aveqg
      qbgZHqbR=-radiLL_bc(7,1,2,5,6,4,3,p(5,4),p(6,4))*aveqg

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

      function radiLL_bc(j1,j2,j3,j4,j5,j6,j7,E1,E2)
      implicit none
      include 'types.f'
      real(dp):: radiLL_bc
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'hdecaymode.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'hbbparams.f'
      integer:: j1,j2,j3,j4,j5,j6,j7
      real(dp):: s45,s12,s13,s23,s123,prop,E1,E2
      real(dp):: fac,hdecay,msqhbb,msqhtautau,msqhgamgam

      s12=s(j1,j2)
      s13=s(j1,j3)
      s23=s(j2,j3)
      s123=s12+s13+s23
c---calculate the 2 Z propagators
      prop=         ((s123-zmass**2)**2+(zmass*zwidth)**2)
      prop=prop*((s(j6,j7)-zmass**2)**2+(zmass*zwidth)**2)
      fac=four*V*gsq*esq**2
      fac=fac*gwsq*wmass**2/(one-xw)**2/prop
      s45=s(j4,j5)
C   Deal with Higgs decay
      call Hbb_bc_bit(s45,E1,E2,1,hdecay)
      hdecay=hdecay/((s45-hmass**2)**2+(hmass*hwidth)**2)
      fac=fac*hdecay
     
     
      
      radiLL_bc=s12/s13/s23
     & *(2._dp*s(j1,j7)*s(j2,j6)+s(j1,j7)*s(j3,j6)+s(j2,j6)*s(j3,j7))
     & +(s(j1,j7)*s(j2,j6)+s(j2,j6)*s(j3,j7)-s(j1,j6)*s(j1,j7))/s13
     & +(s(j1,j7)*s(j2,j6)+s(j1,j7)*s(j3,j6)-s(j2,j6)*s(j2,j7))/s23

      radiLL_bc=fac*radiLL_bc
      return
      end
