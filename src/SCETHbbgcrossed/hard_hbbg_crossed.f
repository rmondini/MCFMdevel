      subroutine hard_hbbg_crossed(order,nproc,uvar,vvar,musq,hard)
      implicit none
!---- hard function for crossed processes of H->b(p1)bbar(p2)g(p3), units of (as/2/pi)**n * LOmsq
!---- nproc=1: bbar(p1)b(p2)->H+g(p3)
!---- nproc=2: b(p2)g(p3)->H+b(p1)
      include 'types.f'
      include 'constants.f'
      include 'masses.f'
      include 'nf.f'
      include 'scet_const.f'
      integer,intent(in)::order,nproc
      real(dp),intent(in)::uvar,vvar,musq
      real(dp),intent(out)::hard(2)
      real(kind=dp) logmh2omu2
      real(kind=dp) nlohe0,nnlohe0,nloxnlohe0

      logmh2omu2 = log(hmass**2/musq)
      hard(:)=zip

      if (order < 1) return

!==== nlo coefficient

      if (nproc.eq.1) then
         call nloevalfinite_proc1(uvar,vvar,musq,nlohe0)
      else if (nproc.eq.2) then
         call nloevalfinite_proc2(uvar,vvar,musq,nlohe0)
      endif

      hard(1)=nlohe0

      if (order < 2) return

!==== nlo^2+nnlo coefficient

      if (nproc.eq.1) then
         call nnloevalfinite_proc1(uvar,vvar,musq,nnlohe0)
      else if (nproc.eq.2) then
         call nnloevalfinite_proc2(uvar,vvar,musq,nnlohe0)
      endif

      if (nproc.eq.1) then
         call nloxnloevalfinite_proc1(uvar,vvar,musq,nloxnlohe0)
      else if (nproc.eq.2) then
         call nloxnloevalfinite_proc2(uvar,vvar,musq,nloxnlohe0)
      endif

      hard(2) = nnlohe0+nloxnlohe0
      
      return
      end

