!---- CW May 16
!---- master formulas for amplitudes for
!------g(i1)+g(i2)+g(i3)+gamma(i4)

      subroutine fill_gg_ga(p,i1,i2,i3,i4,msqgg)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'nf.f'
      include 'ewcharge.f'
      include 'qcdcouple.f'
      real(dp):: p(mxpart,4),msqgg
      integer i1,i2,i3,i4
      complex(dp):: gg_gga(2,2,2,2)
      real(dp):: fac
      complex(dp) :: gg_gga_mmpp,gg_gga_mppp,gg_gga_pppp
      complex(dp) :: gg_gga_pmpm,gg_gga_mppm
      integer h1,h2,h3,h4,j
      real(dp):: qsum
      
      msqgg=zip
      gg_gga(:,:,:,:)=czip
      qsum=zip
      do j=1,nf
         qsum=qsum+Q(j)
      enddo
      
!-----sqaure of factor taken out of KC + (1/16pi^2) from loop
      fac=qsum**2*esq*gsq*2._dp*(ason4pi)**2
      
!---- color factor 
      fac=fac*xn*V
      
!-----debug KC check
!      include 'kc_check_gg_gga.f'
!      call spinorz(4,p,za,zb)
!-----end debug
      call spinoru(4,p,za,zb)
!---- basic amplitudes
      gg_gga(2,2,2,2)=gg_gga_pppp(i1,i2,i3,i4,za,zb)
      gg_gga(1,2,2,2)=gg_gga_mppp(i1,i2,i3,i4,za,zb)
      gg_gga(1,1,2,2)=gg_gga_mmpp(i1,i2,i3,i4,za,zb)
      gg_gga(2,1,2,1)=gg_gga_pmpm(i1,i2,i3,i4,za,zb)
      gg_gga(1,2,2,1)=gg_gga_mppm(i1,i2,i3,i4,za,zb)
  
!----- swaps
      gg_gga(2,1,2,2)=gg_gga_mppp(i2,i3,i1,i4,za,zb)
      gg_gga(2,2,1,2)=gg_gga_mppp(i3,i1,i2,i4,za,zb)
      gg_gga(2,2,2,1)=gg_gga_mppp(i4,i1,i2,i3,za,zb)

      gg_gga(2,1,1,2)=gg_gga_mmpp(i2,i3,i1,i4,za,zb)
      gg_gga(2,2,1,1)=gg_gga_mmpp(i3,i4,i1,i2,za,zb)

!----- conjugates
      gg_gga(1,1,1,1)=gg_gga_pppp(i1,i2,i3,i4,zb,za)
      gg_gga(2,1,1,1)=gg_gga_mppp(i1,i2,i3,i4,zb,za)
      gg_gga(2,2,1,1)=gg_gga_mmpp(i1,i2,i3,i4,zb,za)

      gg_gga(1,2,1,1)=gg_gga_mppp(i2,i3,i1,i4,zb,za)
      gg_gga(1,1,2,1)=gg_gga_mppp(i3,i1,i2,i4,zb,za)
      gg_gga(1,1,1,2)=gg_gga_mppp(i4,i1,i2,i3,zb,za)

      gg_gga(1,2,1,2)=gg_gga_mmpp(i1,i3,i2,i4,zb,za)
      gg_gga(1,2,1,2)=gg_gga_pmpm(i1,i2,i3,i4,zb,za)
      gg_gga(2,2,1,1)=gg_gga_mmpp(i1,i2,i3,i4,zb,za)
      gg_gga(2,1,1,2)=gg_gga_mppm(i1,i2,i3,i4,zb,za)
  
      msqgg=zip
      do h1=1,2
         do h2=1,2
            do h3=1,2
               do h4=1,2
                  msqgg=msqgg
     &        +real(gg_gga(h1,h2,h3,h4)*conjg(gg_gga(h1,h2,h3,h4)),dp)             
!     write(6,*) h1,h2,h3,h4,gg_gga(h1,h2,h3,h4)
               enddo
            enddo
         enddo
      enddo

      msqgg=msqgg*fac
      return
      end
      
      
      function gg_gga_mmpp(i1,i2,i3,i4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: gg_gga_mmpp
      include 'mxpart.f'
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'

      integer i1,i2,i3,i4
      complex(dp):: Ls0,lnrat,boxcoeff,boxint,bub,rat,box

      boxcoeff= -(((za(i1,i4)*zb(i3,i2)*zb(i4,i1)**2 - 
     -        za(i2,i4)*zb(i3,i1)*zb(i4,i2)**2)*zb(i4,i3))/
     -    (2._dp*za(i1,i2)*zb(i2,i1)**4))
      boxint=-lnrat(-s(i1,i3),-s(i1,i4))**2-pi**2

      box=boxcoeff*boxint
      bub= -(lnrat(-s(i1,i3),-s(i1,i4))*za(i1,i3)**2*zb(i2,i1)*
     -    (za(i2,i3)*zb(i3,i2) - za(i2,i4)*zb(i4,i2)))/
     -     (za(i1,i2)**3*zb(i4,i2)**2)
      rat=zb(i4,i3)**2/zb(i2,i1)**2

      gg_gga_mmpp=box+bub+rat
      return
      end

       function gg_gga_mppm(i1,i2,i3,i4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: gg_gga_mppm
      include 'mxpart.f'
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'

      integer i1,i2,i3,i4
      complex(dp):: Ls0,lnrat,boxcoeff,boxint,bub,rat,box

      boxcoeff=   (zb(i3,i2)*(za(i2,i4)*zb(i3,i1)*zb(i4,i2)**2 - 
     -      za(i1,i2)*zb(i2,i1)**2*zb(i4,i3)))/
     -  (2.*za(i1,i4)*zb(i4,i1)**4)
      boxint=-lnrat(-s(i1,i3),-s(i1,i2))**2-pi**2

      box=boxcoeff*boxint
      bub= -(lnrat(-s(i1,i3),-s(i1,i2))*za(i1,i4)*
     -    zb(i3,i1)**2*(-(za(i2,i4)*zb(i4,i2)) + 
     -      za(i3,i4)*zb(i4,i3)))/
     -  (za(i2,i4)**2*zb(i4,i1)**3)

      rat=zb(i3,i2)**2/zb(i4,i1)**2
      
      gg_gga_mppm=box+bub+rat
!      write(6,*) 'tota',gg_gga_mppm
!      write(6,*) 'box',box
!      write(6,*) 'bub',bub
!      write(6,*) 'rat',rat
      return
      end


      function gg_gga_pmpm(i1,i2,i3,i4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: gg_gga_pmpm
      include 'mxpart.f'
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'

      integer i1,i2,i3,i4
      complex(dp):: Ls0,lnrat,boxcoeff,boxint,bub,rat,box
      boxcoeff=   -(za(i2,i4)*(-(za(i1,i4)**2*za(i2,i3)*
     -          zb(i4,i1)) - 
     -       za(i1,i2)*za(i3,i4)**2*zb(i4,i3)))/
     -  (2.*za(i1,i3)**4*zb(i3,i1))
      boxint=-lnrat(-s(i1,i2),-s(i1,i4))**2-pi**2
      
      box=boxcoeff*boxint
      
      bub=  (lnrat(-s(i1,i4),-s(i1,i2))*zb(i1,i3)*
     -    za(i2,i1)**2*(za(i2,i3)*zb(i3,i2) - 
     -      za(i3,i4)*zb(i4,i3)))/
     -  (zb(i3,i4)**2*za(i3,i1)**3)

      rat=za(i4,i2)**2/za(i3,i1)**2

      gg_gga_pmpm=box+bub+rat
      
      return
      end


      
      function gg_gga_mppp(i1,i2,i3,i4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: gg_gga_mppp
      include 'mxpart.f'
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'

      integer i1,i2,i3,i4

      gg_gga_mppp=-(zb(i3,i2)*zb(i4,i2)*zb(i4,i3)**2)/
     -  (s(i1,i2)*zb(i3,i1)*zb(i4,i1))
      return
      end


      
      function gg_gga_pppp(i1,i2,i3,i4,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: gg_gga_pppp
      include 'mxpart.f'
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'

      integer i1,i2,i3,i4

      gg_gga_pppp=-(zb(i1,i2)*zb(i3,i4))/(za(i1,i2)*za(i3,i4))

      return
      end
