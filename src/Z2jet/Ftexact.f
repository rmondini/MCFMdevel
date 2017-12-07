      function Ftexact(s23,mt2)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'qlfirst.f'
      complex(dp)::Ftexact
      complex(dp)::qlI3,qlI2
      real(dp)::s23,mt2,musq

      if (qlfirst) then
        call qlinit
        qlfirst=.false.
      endif

c--- musq is irrelevant, set it to some value
      musq=abs(s23)
      
      Ftexact=
     & -(one+six*mt2*qlI3(s23,0._dp,0._dp,mt2,mt2,mt2,musq,0)
     &   +12._dp*mt2/s23*(qlI2(s23,mt2,mt2,musq,0)-qlI2(0._dp,mt2,mt2,musq,0)))
      
      return
      end
