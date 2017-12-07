      subroutine newton(i,j,k,nDn,phi,s)
!  find limit on s at fixed phi
!  solves for the condition A(i,j,k)=s
      implicit none
      include 'types.f'
      integer::i,j,k
      real(dp)::nDn(3,3),phi,s,f,fp,shift
      real(dp),save::preci
      f(s)=(nDn(i,k)+s*nDn(j,k)-2*cos(phi)*sqrt(s*nDn(i,k)*nDn(j,k)))
     & /nDn(i,j)-s  
      fp(s)=(nDn(j,k)-cos(phi)*sqrt(nDn(i,k)*nDn(j,k)/s))/nDn(i,j)-1._dp 
      preci=1.e-6_dp
      s=1._dp
 20   continue
      shift=-f(s)/fp(s)
      s=s+shift
      write(6,*) 's',s
      if (f(s) .gt. preci) go to 20
      if (s .gt. 1._dp) s=1._dp
      if (s .lt. 0._dp) s=0._dp
c      write(6,*) 's,f(s)',s,f(s)
c      pause
      return
      end
