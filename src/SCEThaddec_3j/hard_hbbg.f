      subroutine hard_hbbg(order,yvar,zvar,musq,hard)
      implicit none
!---- hard function for H->bbg, units of (as/2/pi)**n * LOmsq
      include 'types.f'
      include 'constants.f'
      include 'masses.f'
      include 'nf.f'
      include 'scet_const.f'
      integer,intent(in)::order
      real(dp),intent(in)::yvar,zvar,musq
      real(dp),intent(out)::hard(2)
      real(kind=dp) logmh2omu2
      real(kind=dp) nlohe0,nnlohe0,nloxnlohe0
      logical inclnnlohard

!---------------------
!==== whether to include time-consuming NNLO piece
      inclnnlohard=.true.
!---------------------

      logmh2omu2 = log(hmass**2/musq)
      hard(:)=zip

      if (order < 1) return

!==== nlo coefficient
      call nloheval(yvar,zvar,logmh2omu2,nlohe0)
      hard(1)=nlohe0

      if (order < 2) return

!==== nlo^2+nnlo coefficient
      if(inclnnlohard) then
!         call nnloevalfinite(yvar,zvar,logmh2omu2,nnlohe0)
         call nnloevalfinitetdhpl(yvar,zvar,musq,nnlohe0)
      else
         nnlohe0=0._dp
      endif
      call nloxnlofeval(yvar,zvar,musq,nloxnlohe0)
      hard(2) = nnlohe0+nloxnlohe0
      
      return
      
      return
      end

