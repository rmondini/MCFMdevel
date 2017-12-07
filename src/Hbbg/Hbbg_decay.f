!---  routine for assembling the LO Hbbg_decay routines (from helicity amplitudes to aid NLO calculation)

!---- CW March 17
!-----rotuine is for H=>(b(i1)+bbar(i2)+g(i3))
!---  Note that helicity conservation requires (b^(h1), bbar^(h1))
!---- b quarks are massless
      subroutine Hbbg_decay(p,i1,i2,i3,hdecay)
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
      integer i1,i2,i3,imax
      real(dp) :: p(mxpart,4)
      complex(dp) :: amp(2,2),amp_11,amp_12
      real(dp) :: fac,hdecay,s123
      integer h1,h2
      real(dp) :: babis,s12,s23,s13
      real(dp) :: babisHbbg_msq
      amp(:,:)=czip
      imax=max(i1,i2,i3)
!---- overall factor
      
      fac=V*mb_eff**2*gsq*gwsq/wmass**2/four

    
      call spinoru(imax,p,za,zb) 

      s123=s(i1,i2)+s(i2,i3)+s(i1,i3)
      s12=s(i1,i2)
      s23=s(i2,i3)
      s13=s(i1,i3)
      amp_11 =s123/(zb(i1,i3)*zb(i2,i3))
      amp_12 =za(i1,i2)**2/(za(i1,i3)*za(i2,i3))

      amp(1,1)=amp_11
      amp(2,2)=s123/(za(i1,i3)*za(i2,i3))

      amp(1,2)=amp_12
      amp(2,1)=zb(i1,i2)**2/(zb(i1,i3)*zb(i2,i3))

      hdecay=zip

      do h1=1,2
         do h2=1,2
            hdecay=hdecay+abs(amp(h1,h2))**2
         enddo
      enddo
!      babis=mb_eff**2*(4*CF*gsq*(s123**2 + s12**2)*(s13 + s23)*xn)/
!     &     ((s123 - s12)*s13*s23)/vevsq*2._dp

      
      hdecay=hdecay*fac
!      babis= babisHbbg_msq(i1,i2,i3)
!      write(6,*) 'b check'
!      write(6,*) babis,hdecay,babis/hdecay
!      write(6,*) '****'
!      pause
      return
      end
      
