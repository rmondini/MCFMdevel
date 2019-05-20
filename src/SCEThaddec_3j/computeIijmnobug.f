      subroutine computeIijmnobug(y12,y23,y31,Iijm)
      implicit none
c--- One-time computation of integrals required for one-jet soft functions
c--- RM Nov 18: we fix the bug n=1,5->n=1,6 and Iijm(5)->Iijm(6) and use this file
c---            in soft_dec_3j.f. The files in the SCET1j folder have not been
c---            modified so they still contain the typo.
      include 'types.f'
      include 'constants.f'
      real(dp)::y12,y23,y31,y(3,3),Iijm(6),al,be,I0JSTW,I1JSTW,I0,I1
      integer::n
      integer,parameter::i(6)=(/1,2,2,3,3,1/)
      integer,parameter::j(6)=(/2,1,3,2,1,3/)
      integer,parameter::m(6)=(/3,3,1,1,2,2/)

      y(:,:)=zip

      y(1,2)=y12
      y(2,1)=y12
      y(2,3)=y23
      y(3,2)=y23
      y(1,3)=y31
      y(3,1)=y31

      do n=1,6
      al=y(j(n),m(n))/y(i(n),j(n))
      be=y(i(n),m(n))/y(i(n),j(n))
      I0=I0JSTW(al,be)
      I1=I1JSTW(al,be)
      Iijm(n)=I0*log(al)+I1
      enddo

      return
      end
      
      
