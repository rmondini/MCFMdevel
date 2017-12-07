      subroutine scaleset_Msqpt3sq(p,mu0)
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c---  sqrt(M^2+pt3^2), where M is mass2
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'kprocess.f'
      include 'breit.f'
      real(dp):: p(mxpart,4),mu0,pt

      if((kcase==ktt_tot) .or.
     &   (kcase==kbb_tot) .or.
     &   (kcase==kcc_tot)) then
        mu0=mass2**2+pt(3,p)**2
        mu0=sqrt(abs(mu0))
      else
        write(6,*) 'dynamicscale sqrt(M^2+pt3^2)'//
     &             ' not supported for this process.'
        stop
      endif
      
      return
      end
      

