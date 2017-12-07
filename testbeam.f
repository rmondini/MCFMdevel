      program testbeam
! Attempt to replicate beam-function plots in Fig.3 of arXiv:1405.1044
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'pdlabel.f'
      include 'couple.f'
      include 'scale.f'
      include 'facscale.f'
      include 'selectpdfs.f'
      include 'masses.f'
      integer iflav,ix
      real(dp):: muB,x,eps,dgauss,B0,B1,B2,fx(-5:5),alphas,ason2pi
      common/setupvars/iflav,muB,x
      external integrand0,integrand1,integrand2

! mass parameters for running routines
      zmass=91.1876_dp

! PDF set MSTW2008 NNLO
      pdlabel='mstw8nn'
      call pdfwrap
      
! muB == sqrt(tmax)
      muB=30._dp

! scales
      facscale=muB
      scale=muB
      musq=scale**2

! accuracy of integration
      eps=5.e-3_dp

! include all parton flavors
      selectpdfs(:)=99
      
! parton flavor
      iflav=0
      
! alpha_s(muB) with 2-loop running
      ason2pi=alphas(muB,amz,3)/two/pi
      write(6,*) 'parton flavor',iflav
      write(6,*) 'alpha_s(MZ)=',amz
      write(6,*) 'alpha_s(muB)=',ason2pi*two*pi

      do ix=-45,-1
      x=10._dp**(ix/15._dp)
      
      selectpdfs(1)=99 ! select all pdfs
      call fdist_select(1,1,x,muB,fx)
      B0=fx(iflav)
      selectpdfs(1)=-1  ! isolate a->b contributions here
! Note: change of integral from [0,1] to [x,1] takes place inside xbeam
      B1=ason2pi*dgauss(integrand1,0._dp,1._dp,eps)
      B2=ason2pi**2*dgauss(integrand2,0._dp,1._dp,eps)
!      write(6,*) 'x,B0,B1,B2',x,B0,B1,B2
      write(6,*) 'x, B1/B0, B2/B0 [%]',x,B1/B0*100._dp,B2/B0*100._dp
      
      enddo
      
      
      stop
      end
 
! This is the 1-loop function to integrate
      function integrand1(z)
      implicit none
      include 'types.f'
      integer iflav
      real(dp):: integrand1,z,muB,x,btau(-5:5,-1:1)
      common/setupvars/iflav,muB,x

      call xbeam1bis(1,1,z,x,muB,btau)
      
      integrand1=btau(iflav,-1)
      
      return
      end
      
! This is the 2-loop function to integrate
      function integrand2(z)
      implicit none
      include 'types.f'
      integer iflav
      real(dp):: integrand2,z,muB,x,btau(-5:5,-1:3)
      common/setupvars/iflav,muB,x

      call xbeam2bis(1,1,z,x,muB,btau)
      
      integrand2=btau(iflav,-1)
      
      return
      end
      
