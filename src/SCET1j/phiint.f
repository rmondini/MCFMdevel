      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp)::nDn(3,3),soft,soft123,n1Dn2,n1Dn3,n2Dn3,tauC,
     & soft121
      
      tauC=10.2123232356_dp
      n1Dn2=0.112345_dp
      n1Dn3=0.2345_dp
      n2Dn3=0.45112345_dp

      nDn(:,:)=zip
      nDn(1,2)=n1Dn2
      nDn(1,3)=n1Dn3
      nDn(2,3)=n2Dn3
      nDn(2,1)=nDn(1,2)
      nDn(3,1)=nDn(1,3)
      nDn(3,2)=nDn(2,3)

      soft123=soft(1,2,3,nDn,tauC)
      soft121=soft(1,2,1,nDn,tauC)
      write(6,*) 'soft121',soft121
      stop
      end

      function soft(i,j,k,nDn,tauC)
!     calculate soft function for jettiness
!     eikonal between i and j 
!     jetiness minimum wrt to line k
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp)::soft
      real(dp)::nDn(3,3),ndotn(3,3),tauC,tauCC,egauss,
     & phiintA,phiintB
      integer::i,j,k,ii,jj,kk
      common/ijk/ii,jj,kk
      common/naDnb/ndotn
      common/tauC/tauCC
      external phiintA,phiintB
      tauCC=tauC
      ndotn(:,:)=nDn
      ii=i
      jj=j
      kk=k
      if ((k .eq. i) .or. (k.eq.j)) then
      soft=egauss(phiintA,0._dp,twopi,0.00001_dp)
      else
      soft=egauss(phiintB,0._dp,twopi,0.00001_dp)
      endif
      return
      end
    
c      function I1real(n1Dn2,n)
c      implicit none
c      include 'types.f'
c      real(dp)::I1real
c      real(dp)::egauss,n1Dn2,n1Dn3,n2Dn3,n1n2,n1n3,n2n3,phiint
c      common/ijk/ii,jj,kk
c      common/naDnb/n1n2,n1n3,n2n3
c      external phiint
c      n1n2=n1Dn2
c      n1n3=n1Dn3
c     n2n3=n2Dn3
c      I1real=egauss(phiint,0._dp,1._dp,0.00001_dp)
c      return
c      end


      function phiintA(x)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp) x,phiintA,phi,phip,dgauss,sintA,sinphi,tauC
      external sintA
      real(dp)::nDn(3,3)
      integer::i,j,k
      common/ijk/i,j,k
      common/naDnb/nDn
      common/phi/phip
      common/tauC/tauC
      phi=twopi*x
      phip=phi
      sinphi=sin(phi)
      phiintA=dgauss(sintA,0._dp,1._dp,0.00001_dp)

C       + ds * (
C          + 1/2*theta(s,Aijk)*ln(s)*s^-1
C          - theta(s,Aijk)*ln(sinphi)*s^-1
C          - theta(s,Aijk)*ln(tauc)*s^-1
C          + 1/2*theta(s,Aijk)*ln(nDn(i,j))*s^-1
C          )

      phiintA=phiintA
     &  - log(sinphi)**2
     &     - 2._dp*log(sinphi)*log(tauc)
     &     + log(sinphi)*log(nDn(i,j))
     &     - log(tauc)**2
     &     + log(tauc)*log(nDn(i,j))
     &     - 0.25_dp*log(nDn(i,j))**2

      return
      end




      function sintA(s)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp) s,phi,sintA,slim,rts,root
      real(dp)::nDn(3,3),tauC,cosphi
      integer::i,j,k
      common/ijk/i,j,k
      common/naDnb/nDn
      common/tauC/tauC
      common/phi/phi

      cosphi=cos(phi)
      root=(cosphi**2*nDn(i,j)-nDn(i,k))*nDn(j,k)+nDn(i,j)*nDn(i,k)
      write(6,*) 'cosphi,root',cosphi,root
      if (root > 0) then
      rts =(cosphi*sqrt(nDn(i,j)*nDn(j,k))+sqrt(root))
     & /(nDn(j,k)-nDn(i,j))
      write(6,*) 'rts',rts
      rts =(cosphi*sqrt(nDn(i,j)*nDn(j,k))-sqrt(root))
     & /(nDn(j,k)-nDn(i,j))
      write(6,*) 'rts',rts
       
      slim=rts**2
      if (s .gt. slim) then
      sintA=(half*log(s)+half*log(nDn(i,j))-log(sin(phi))-log(tauC))/s
      else 
      sintA=0._dp
      endif
      else
      sintA=0._dp
      endif
!      call newton(i,j,k,nDn,phi,slim)
!       + ds * (
!          + 1/2*theta(s,Aijk)*ln(s)*s^-1
!          - theta(s,Aijk)*ln(sinphi)*s^-1
!          - theta(s,Aijk)*ln(tauc)*s^-1
!          + 1/2*theta(s,Aijk)*ln(nDn(i,j))*s^-1
!          )
      return
      end

      function phiintB(x)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp) x,phiintB,phi,phip,dgauss,sintB,slim,rts,cosphi,root
      real(dp)::nDn(3,3)
      integer::i,j,k
      common/ijk/i,j,k
      common/naDnb/nDn
      common/phi/phip
      external sintB

      phi=twopi*x
      phip=phi
      cosphi=cos(phi)
      write(6,*) 'cosphi',cos(phi)
      root=(cosphi**2*nDn(k,i)-nDn(k,j))*nDn(i,j)+nDn(k,i)*nDn(k,j)
      if (root > 0) then
      rts =(cosphi*sqrt(nDn(k,i)*nDn(i,j))+sqrt(root))
     & /(nDn(i,j)-nDn(k,i))
      write(6,*) 'rts',rts
      rts =(cosphi*sqrt(nDn(k,i)*nDn(i,j))-sqrt(root))
     & /(nDn(i,j)-nDn(k,i))
      write(6,*) 'rts',rts
      slim=rts**2
c      call newton(k,i,j,nDn,phi,slim)
      phiintB=dgauss(sintB,0._dp,slim,0.00001_dp)
      else
      phiintB=0
      endif     
      return
      end

      function sintB(s)
      implicit none
      include 'types.f'
      real(dp)::sintB,s
      real(dp)::A,nDpj,phi,nDn(3,3),tauC
      integer::i,j,k
      common/phi/phi 
      common/ijk/i,j,k
      common/naDnb/nDn
      common/tauC/tauC
      A(i,j,k,s)=(nDn(i,k)+s*nDn(j,k)
     & -2*cos(phi)*sqrt(s*nDn(i,k)*nDn(j,k)))/nDn(i,j)  
      nDpj=A(k,i,j,s)
      if (nDpj .gt. s) then
      sintB=nDn(i,j)/(nDn(i,k)*nDpj)
     & *(log(sin(phi)*tauc)-0.5_dp*log(s*nDn(i,k)))
      else
      sintB=0._dp
      endif
      return
      end

