      subroutine gg_hg_vbis(p,msq,order)
      implicit none
      include 'types.f'
c---- Matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H(-->  b(p3)+b~(p4))+g(p5)
c
c---- Computed using amplitudes extracted from PeTer
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'hdecaymode.f'
      integer:: j,k,iglue,order
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: ss,tt,uu,mhsq,hdecay
      real(dp):: s34,Asq,fac,msqgamgam
      real(dp):: gg,qqb,qg,gqb
      real(dp):: Hgg(2),Hqqb(2),Hqg(2),Hgqb(2)
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
c--- compute matrix elements
      call Hampgggsq(order,5,1,2,3,4,gg,Hgg(1),Hgg(2),qg,Hqg(1),Hqg(2))
      call Hampgggsq(order,2,5,1,3,4,gg,Hgg(1),Hgg(2),gqb,Hgqb(1),Hgqb(2))
      call Hampgggsq(order,2,1,5,3,4,gg,Hgg(1),Hgg(2),qqb,Hqqb(1),Hqqb(2))
      
c--- pick out higher-order coefficients
      if (order == 1) then
        gg =ason2pi*Hgg(1)*gg
        qqb=ason2pi*Hqqb(1)*qqb
        qg =ason2pi*Hqg(1)*qg
        gqb=ason2pi*Hgqb(1)*gqb
      endif
      if (order == 2) then
        gg =ason2pi**2*Hgg(2)*gg
        qqb=ason2pi**2*Hqqb(2)*qqb
        qg =ason2pi**2*Hqg(2)*qg
        gqb=ason2pi**2*Hgqb(2)*gqb
      endif
            
      Asq=(as/(3._dp*pi))**2/vevsq
      fac=Asq*gsq*hdecay

      msq(0,0)=avegg*fac*V*xn*gg
      msq(1,-1)=+aveqq*fac*V/two*qqb
      msq(0,-1)=-aveqg*fac*V/two*gqb
      msq(+1,0)=-aveqg*fac*V/two*qg

      do j=-nf,nf
      do k=-nf,nf
      if ((k == -j) .and. (j .ne. 0)) then
      msq(j,k)=msq(1,-1)
      elseif ((j == 0) .and. (k .ne. 0)) then
      msq(j,k)=msq(0,-1)
      elseif ((j .ne. 0) .and. (k == 0)) then
      msq(j,k)=msq(1,0)
      endif
      enddo
      enddo
      
      return
      end
