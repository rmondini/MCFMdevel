!---- assembly routine for Ql pieces
!---- fills amplitudes and returns lc and slc amplitudes
!---- CW Jan 17 

!-----Order is q(i1)^h1 qb(i2)^(-h1) + l(i3)^h2 +lb(i4)^(-h2)+gamma(i5)^h3 +g(i6)^h4
!---- returns amp(h1,h2,h3,h4)

!---- KC.sty prefactor rt2/2/ee/gw^2*costhw/gs^3


!---- combo is xn * amp_lc + (1/xn)* amp_slc
      subroutine amp_zaj_ql_assemble(i1,i2,i3,i4,i5,i6,za,zb,
     &     amp_lc,amp_slc)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'zprods_decl.f'
      integer i1,i2,i3,i4,i5,i6
      complex(dp) :: amp_lc(2,2,2,2),amp_slc(2,2,2,2)
      complex(dp) :: amp_zaj_ql_lc_MHV,amp_zaj_ql_slc_MHV
     &     ,amp_zaj_ql_lc_NMHV,amp_zaj_ql_slc_NMHV
      integer h1,h2,h3,h4

      
!---- basic ordering MHV and NMV
      amp_lc(2,1,2,2)  = amp_zaj_ql_lc_MHV(i1,i2,i3,i4,i5,i6,za,zb)
      amp_slc(2,1,2,2) = amp_zaj_ql_slc_MHV(i1,i2,i3,i4,i5,i6,za,zb)
      amp_lc(2,1,1,2)  = amp_zaj_ql_lc_NMHV(i1,i2,i3,i4,i5,i6,za,zb)
      amp_slc(2,1,1,2) = amp_zaj_ql_slc_NMHV(i1,i2,i3,i4,i5,i6,za,zb)
      
!----- conjugation 
      amp_lc(1,2,1,1)  = amp_zaj_ql_lc_MHV(i1,i2,i3,i4,i5,i6,zb,za)
      amp_slc(1,2,1,1) = amp_zaj_ql_slc_MHV(i1,i2,i3,i4,i5,i6,zb,za)
      amp_lc(1,2,2,1)  = amp_zaj_ql_lc_NMHV(i1,i2,i3,i4,i5,i6,zb,za)
      amp_slc(1,2,2,1) = amp_zaj_ql_slc_NMHV(i1,i2,i3,i4,i5,i6,zb,za)

!----  single line reversals
      amp_lc(1,1,2,2)  = -amp_zaj_ql_lc_MHV(i2,i1,i3,i4,i5,i6,za,zb)
      amp_slc(1,1,2,2) = -amp_zaj_ql_slc_MHV(i2,i1,i3,i4,i5,i6,za,zb)
      amp_lc(1,1,1,2)  = -amp_zaj_ql_lc_NMHV(i2,i1,i3,i4,i5,i6,za,zb)
      amp_slc(1,1,1,2) = -amp_zaj_ql_slc_NMHV(i2,i1,i3,i4,i5,i6,za,zb)

      amp_lc(2,2,2,2)  = -amp_zaj_ql_lc_MHV(i1,i2,i4,i3,i5,i6,za,zb)
      amp_slc(2,2,2,2) = -amp_zaj_ql_slc_MHV(i1,i2,i4,i3,i5,i6,za,zb)
      amp_lc(2,2,1,2)  = -amp_zaj_ql_lc_NMHV(i1,i2,i4,i3,i5,i6,za,zb)
      amp_slc(2,2,1,2) = -amp_zaj_ql_slc_NMHV(i1,i2,i4,i3,i5,i6,za,zb)

!----- double line reverasls 
      amp_lc(1,2,2,2)  = amp_zaj_ql_lc_MHV(i2,i1,i4,i3,i5,i6,za,zb)
      amp_slc(1,2,2,2) = amp_zaj_ql_slc_MHV(i2,i1,i4,i3,i5,i6,za,zb)
      amp_lc(1,2,1,2)  = amp_zaj_ql_lc_NMHV(i2,i1,i4,i3,i5,i6,za,zb)
      amp_slc(1,2,1,2) = amp_zaj_ql_slc_NMHV(i2,i1,i4,i3,i5,i6,za,zb)


!---- single line reversals + conj
      amp_lc(2,2,1,1)  = -amp_zaj_ql_lc_MHV(i2,i1,i3,i4,i5,i6,zb,za)
      amp_slc(2,2,1,1) = -amp_zaj_ql_slc_MHV(i2,i1,i3,i4,i5,i6,zb,za)
      amp_lc(2,2,2,1)  = -amp_zaj_ql_lc_NMHV(i2,i1,i3,i4,i5,i6,zb,za)
      amp_slc(2,2,2,1) = -amp_zaj_ql_slc_NMHV(i2,i1,i3,i4,i5,i6,zb,za)

      amp_lc(1,1,1,1)  = -amp_zaj_ql_lc_MHV(i1,i2,i4,i3,i5,i6,zb,za)
      amp_slc(1,1,1,1) = -amp_zaj_ql_slc_MHV(i1,i2,i4,i3,i5,i6,zb,za)
      amp_lc(1,1,2,1)  = -amp_zaj_ql_lc_NMHV(i1,i2,i4,i3,i5,i6,zb,za)
      amp_slc(1,1,2,1) = -amp_zaj_ql_slc_NMHV(i1,i2,i4,i3,i5,i6,zb,za)

      
!---- doble line reversals + conj 

      amp_lc(2,1,1,1)  = amp_zaj_ql_lc_MHV(i2,i1,i4,i3,i5,i6,zb,za)
      amp_slc(2,1,1,1) = amp_zaj_ql_slc_MHV(i2,i1,i4,i3,i5,i6,zb,za)
      amp_lc(2,1,2,1)  = amp_zaj_ql_lc_NMHV(i2,i1,i4,i3,i5,i6,zb,za)
      amp_slc(2,1,2,1) = amp_zaj_ql_slc_NMHV(i2,i1,i4,i3,i5,i6,zb,za)



!---- checking code
      if(1.eq.2) then
         write(6,*) '****** LC amplitudes **********' 
!----  writeout LC
         do h1 =1,2
            do h2=1,2
               do h3=1,2
                  do h4=1,2
                     write(6,*) h1,h2,h3,h4,amp_lc(h1,h2,h3,h4)
                  enddo
               enddo
            enddo
         enddo
         
         write(6,*) '******************************' 
         write(6,*) 
         write(6,*) '****** SLC amplitudes **********' 
         
!----  writeout SLC
         do h1 =1,2
            do h2=1,2
               do h3=1,2
                  do h4=1,2
                     write(6,*) h1,h2,h3,h4,amp_slc(h1,h2,h3,h4)
                  enddo
               enddo
            enddo
         enddo
         write(6,*) '******************************'
      endif
    
      
    
      return
      end
