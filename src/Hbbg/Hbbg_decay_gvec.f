!--- CW March 17, Routine to compute the H=> bb g virtual matrix element squared

      subroutine Hbbg_decay_gvec_ex(p,n,ipt,msqv)
      implicit none
      include 'types.f'
!---- wrapper routine to gvec 
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'hbbparams.f'
      include 'epinv.f'
      include 'scale.f'
      real(dp) :: msqv,p(mxpart,4),n(4)
      integer i1,i2,i3,ipt
      complex(dp) :: amp_p(2)
      integer h1,h2,imax
      complex(dp) :: zab(mxpart,mxpart),zba(mxpart,mxpart)
      real(dp) :: fac
      
      if(ipt.ne.3) then
         write(6,*) 'Error in Higgs decay gvec'
         write(6,*) 'ipt != 3',ipt
         stop
      endif

      call Hbbg_decay_gvec(p,1,2,3,n,msqv)
      return
      end
      
      subroutine Hbbg_decay_gvec(p,i1,i2,i3,n,msqgv)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'hbbparams.f'
      include 'epinv.f'
      include 'scale.f'
      real(dp) :: msqgv,p(mxpart,4),n(4)
      integer i1,i2,i3
      complex(dp) :: amp_p(2)
      integer h1,h2,imax
      complex(dp) :: zab(mxpart,mxpart),zba(mxpart,mxpart)
      real(dp) :: fac
      
      msqgv=zip
!---- LO factor 
      fac=V*mb_eff**2*gsq*gwsq/wmass**2/8._dp
      call spinoru(i3,p,za,zb)
      call spinork(i3,p,zab,zba,n) 
      call checkndotp(p,n,i3)


      amp_p(1)=(za(i1,i2)*zab(i1,i1))/s(i1,i3) - 
     -  (za(i2,i3)*zab(i1,i3))/s(i1,i3) - 
     -  (za(i1,i2)*zab(i2,i2))/s(i2,i3)
     &     - (za(i1,i3)*zab(i2,i3))/s(i2,i3)
      amp_p(2)=-((zb(i2,i1)*zba(i1,i1))/s(i1,i3)) + 
     -  (zb(i3,i2)*zba(i1,i3))/s(i1,i3) + 
     -  (zb(i2,i1)*zba(i2,i2))/s(i2,i3) 
     &     + (zb(i3,i1)*zba(i2,i3))/s(i2,i3)


      do h1=1,2
         msqgv=msqgv+real(amp_p(h1)*conjg(amp_p(h1)),dp)
      enddo

      msqgv=msqgv*fac
      
      return
      end
