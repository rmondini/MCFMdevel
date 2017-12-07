      function getdynamictau(p)
      implicit none
      include 'types.f'
      include 'kprocess.f'
      include 'mxpart.f'
      include 'taucut.f'
      real(dp):: getdynamictau,p(mxpart,4),pt

! Return normal value of tau if not(dynamictau)
      if (dynamictau .eqv. .false.) then
        getdynamictau=taucut
        return
      endif
      
      if ( (kcase == kdirgam) .or. (kcase == kgamjet)
     & .or.(kcase == kgam_2j)) then
        getdynamictau=pt(3,p)
      else
        write(6,*) 'Invalid process for dynamic taucut!'
        stop
      endif
      
      getdynamictau=getdynamictau*taucut
      
      return
      end
      
