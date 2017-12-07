      subroutine gen6vh(r,q,wt6,*)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'mxdim.f'
      include 'debug.f'
      include 'phasemin.f'
      include 'limits.f'
      include 'kprocess.f'
      include 'masses.f'
      include 'jetcuts.f'
      include 'x1x2.f'
      include 'energy.f'
      integer:: nu
      real(dp):: r(mxdim)
      real(dp):: wt6,q(mxpart,4)
      real(dp):: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4),p8(4)
      real(dp) :: Qsqmin,Qsqmax
      real(dp):: y,pswt,xjac,tau,vs,vsqmax,vsqmin,
     & s34,rtshat,ymax,yave
      real(dp), parameter:: Qsqmincut = 0.1_dp
      real(dp) :: lntaum,x1mx2,surd,muk,rho
      real(dp) :: xsoft,theta,phi
      
      wt6=0._dp
      Qsqmin=max(m3456min**2,Qsqmincut)
      Qsqmax=m3456max**2
      
c--- generic process

      lntaum=log(taumin)
      tau=exp(lntaum*(one-r(9)))
      xjac=-lntaum*tau

c      tau=(one-taumin)*r(9)**2+taumin
c      xjac=2*r(9)*(one-taumin)

      x1mx2=two*r(10)-one
      surd=sqrt(x1mx2**2+four*tau) 
           
      xx(1)=half*(+x1mx2+surd)
      xx(2)=half*(-x1mx2+surd)
!      write(*,*) r(1),r(9),r(10)
!      write(*,*) xx(1),xx(2)

      xjac=xjac*two/surd


!      write(6,*) xmin,xx(1)
!      write(6,*) xmin,xx(2)
!      pause
c---  phase space volume only checked for x1=x2=1
      if ((kcase==kvlchwg) .or. (kcase==kvlchwh)) then
        xx(1)=1._dp
        xx(2)=1._dp
        xjac=1._dp
      endif

c---if x's out of normal range alternative return
      if   ((xx(1) > 1._dp)
     & .or. (xx(2) > 1._dp)
     & .or. (xx(1) < xmin)
     & .or. (xx(2) < xmin)) return 1

      p1(4)=-xx(1)*sqrts*half
      p1(1)=zip
      p1(2)=zip
      p1(3)=-xx(1)*sqrts*half

      p2(4)=-xx(2)*sqrts*half
      p2(1)=zip
      p2(2)=zip
      p2(3)=+xx(2)*sqrts*half

      if     ((kcase==kW_twdk) .or. (kcase==kWtbwdk)
     &   .or. (kcase==kW_cwdk) .or. (kcase==kvlchwh)) then
c--- W+t process, radiation in production
        call phase6a(r,p1,p2,p3,p4,p5,p6,p7,p8,pswt,*999)
      elseif ((kcase==kWtdkay) .or. (kcase==kvlchwg)) then
c--- W+t process, radiation in decay
         call phase6b(r,p1,p2,p3,p4,p5,p6,p7,p8,pswt,*999)
      elseif (kcase==kZHbbjt) then
         call phase6h(r,p1,p2,p3,p4,p5,p6,p7,p8,pswt,*999)
         goto 1113
       else
c--- generic case 
        call phase6(r,p1,p2,p3,p4,p5,p6,p7,p8,pswt,*999)
      endif

      goto 1112
!---  DEBUG Generate soft phase space point for checking
      xsoft=1.e-7_dp
      do nu=1,4
      p7(nu)=xsoft*p7(nu)
!-----adjust p8 to conserve momenutm
      p8(nu)=p8(nu)+(one-xsoft)*p7(nu)/xsoft
      enddo
!---- make p8 massless
      muk = 2._dp*(p8(4)*p5(4)-p8(3)*p5(3)-p8(2)*p5(2)-p8(1)*p5(1))
      rho = p8(4)**2-p8(3)**2-p8(2)**2-p8(1)**2
      p8(:)=-rho/muk*p5(:)+p8(:)
      p5(:)=(one+rho/muk)*p5(:)
 1112 continue

      goto 1113
!-----DEBUG generate a collinear phase space point for checking
      theta=r(9)*twopi
      phi=r(10)*twopi
      xsoft=1.e-4_dp
      p7(2)=+p7(4)*sin(theta)
      p7(3)=+p7(4)*cos(theta)*sin(phi)
      p7(1)=+p7(4)*cos(theta)*cos(phi)
      p8(2)=+p8(4)*sin(theta+xsoft)
      p8(3)=+p8(4)*cos(theta+xsoft)*sin(phi)
      p8(1)=+p8(4)*cos(theta+xsoft)*cos(phi)
 !      goto 1113

!      write(6,*) 'in p5',p5
!-----adjust p8 to conserve momentum 
      do nu=1,4
         p5(nu)=-p1(nu)-p2(nu)-p3(nu)-p4(nu)-p6(nu)-p7(nu)-p8(nu)
      enddo
     
 !     goto 1113
!----- now make p8 massless 
      muk = 2._dp*(p5(4)*p1(4)-p5(3)*p1(3)-p5(2)*p1(2)-p5(1)*p1(1))
      rho = p5(4)**2-p5(3)**2-p5(2)**2-p5(1)**2
      p5(:)=-rho/muk*p1(:)+p5(:)
      p1(:)=(one+rho/muk)*p1(:)

      
!      write(6,*) 's(58) =',p8(4)*p5(4)-p8(3)*p5(3)-p8(2)*p5(2)-p8(1)*p5(1)
!      stop
 1113 continue
      
       
      
      
      do nu=1,4
        q(1,nu)=p1(nu)
        q(2,nu)=p2(nu)
        q(3,nu)=p3(nu)
        q(4,nu)=p4(nu)
        q(5,nu)=p5(nu)
        q(6,nu)=p6(nu)
        q(7,nu)=p7(nu)
        q(8,nu)=p8(nu)
      enddo 
 !     call writeout(q)
 !     write(6,*) xjac,pswt
 !     pause
      wt6=xjac*pswt

!      if (debug) write(6,*) 'wt6 in gen6',wt6
      return
 999  q(:,:)=0._dp
      return 1
      
      end

