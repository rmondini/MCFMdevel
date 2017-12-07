      subroutine gg_hgbis(p,msq)
      implicit none
      include 'types.f'
c---- matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H(-->  b(p3)+b~(p4))+g(p5)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'hdecaymode.f'
      integer:: j,k,iglue
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: ss,tt,uu,mhsq,hdecay
      real(dp):: s34,Asq,fac,msqgamgam
      real(dp):: Hampgggsq,gg
      parameter(iglue=5)

      msq(:,:)=zip

c--- note: important to use this statement (rather than dot product) in
c--- order to retain the possibility of massive decay products
      s34=(p(3,4)+p(4,4))**2
     &   -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2

C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqgamgam(hmass)
      else
      write(6,*) 'Unimplemented process in gg_hg'
      stop
      endif
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

      call dotem(iglue,p,s)
c--- first parameter in H1ampgggsq is order
      gg=Hampgggsq(0,1,2,5,3,4)
      
      Asq=(as/(3._dp*pi))**2/vevsq
      fac=Asq*gsq*hdecay

      msq(0,0)=avegg*fac*V*xn*gg
      msq(1,-1)=+aveqq*fac*V/2._dp*(tt**2+uu**2)/ss
      msq(0,+1)=-aveqg*fac*V/2._dp*(ss**2+tt**2)/uu
      msq(+1,0)=-aveqg*fac*V/2._dp*(ss**2+uu**2)/tt

      do j=-nf,nf
      do k=-nf,nf
      if ((k == -j) .and. (j .ne. 0)) then
      msq(j,k)=msq(1,-1)
      elseif ((j == 0) .and. (k .ne. 0)) then
      msq(j,k)=msq(0,1)
      elseif ((j .ne. 0) .and. (k == 0)) then
      msq(j,k)=msq(1,0)
      endif
      enddo
      enddo
      
      return
      end
