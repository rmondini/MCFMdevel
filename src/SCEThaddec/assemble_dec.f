      function assemble_dec(order,beama0,beamb0,jeta1,jetb1,
     & jeta2,jetb2,soft1,soft2,hard)
      implicit none
!-----This routine is a modification of the standard assembly
!---- style routines, designed for hadronic decays, takes in LO beam functions
!-----and corrections to jet functions and soft functions
!---- CW Mar 17 
c---- Given beam, soft, jet and hard functions, calculates
c---- O(alphas^order) correction after all convolutions
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'scet_const.f'
      include 'taucut.f'
      include 'kpart.f'
      integer:: order
      real(dp):: assemble_dec
      real(dp):: soft1(-1:1),soft2(-1:3),hard(2),
     &     beama0,beamb0,
     & jeta1(-1:1),jeta2(-1:3),
     & jetb1(-1:1),jetb2(-1:3),
     & Lt1cut,full(-1:3)

      Lt1cut = log(taucut/scale)

      if(coeffonly)then
         assemble_dec=zip
      else
         assemble_dec=beama0*beamb0
      endif
      
      
      if((order==1).or.
     & ((order==2).and.(coeffonly.eqv..false.)))then
         full(1)=(beama0*beamb0*jeta1(1))/2._dp
     &        + (beama0*beamb0*jetb1(1))/2._dp + beama0*beamb0*soft1(1)
         
         full(0)=(beama0*beamb0*jeta1(0))/2._dp
     &        + (beama0*beamb0*jetb1(0))/2._dp + beama0*beamb0*soft1(0)
         
         full(-1)=  beama0*beamb0*hard(1)
     &  + (beama0*beamb0*jeta1(-1))/2._dp         
     &  +(beama0*beamb0*jetb1(-1))/2._dp + beama0*beamb0*soft1(-1)
         
      assemble_dec=assemble_dec+ason2pi*(
     & +      full(-1)
     & +full(0)*Lt1cut
     &        +full(1)*Lt1cut**2/two)
     
      endif

      if(order>1)then
         full(3)=jeta2(3)/4. + (jeta1(1)*jetb1(1))/4. + 
     -  jetb2(3)/4. + (jeta1(1)*soft1(1))/2. + 
     -  (jetb1(1)*soft1(1))/2. + soft2(3)
         
      full(2)=  jeta2(2)/4. + (3*jeta1(1)*jetb1(0))/8. + 
     -  (3*jeta1(0)*jetb1(1))/8. + jetb2(2)/4. + 
     -  (3*jeta1(1)*soft1(0))/4. + 
     -  (3*jetb1(1)*soft1(0))/4. + 
     -  (3*jeta1(0)*soft1(1))/4. + 
     -        (3*jetb1(0)*soft1(1))/4. + soft2(2)
      
      full(1)=(hard(1)*jeta1(1))/2. + jeta2(1)/4. + 
     -  (jeta1(1)*jetb1(-1))/4. + (jeta1(0)*jetb1(0))/2. + 
     -  (hard(1)*jetb1(1))/2. + (jeta1(-1)*jetb1(1))/4. - 
     -  (Pi**2*jeta1(1)*jetb1(1))/12. + jetb2(1)/4. + 
     -  (jeta1(1)*soft1(-1))/2. + 
     -  (jetb1(1)*soft1(-1))/2. + jeta1(0)*soft1(0) + 
     -  jetb1(0)*soft1(0) + hard(1)*soft1(1) + 
     -  (jeta1(-1)*soft1(1))/2. - 
     -  (Pi**2*jeta1(1)*soft1(1))/6. + 
     -  (jetb1(-1)*soft1(1))/2. - 
     -  (Pi**2*jetb1(1)*soft1(1))/6. + soft2(1)
      
      full(0)=  (hard(1)*jeta1(0))/2. + jeta2(0)/4. + 
     -  (jeta1(0)*jetb1(-1))/4. + (hard(1)*jetb1(0))/2. + 
     -  (jeta1(-1)*jetb1(0))/4. - 
     -  (Pi**2*jeta1(1)*jetb1(0))/24. - 
     -  (Pi**2*jeta1(0)*jetb1(1))/24. + 
     -  (zeta3*jeta1(1)*jetb1(1))/2. + jetb2(0)/4. + 
     -  (jeta1(0)*soft1(-1))/2. + 
     -  (jetb1(0)*soft1(-1))/2. + hard(1)*soft1(0) + 
     -  (jeta1(-1)*soft1(0))/2. - 
     -  (Pi**2*jeta1(1)*soft1(0))/12. + 
     -  (jetb1(-1)*soft1(0))/2. - 
     -  (Pi**2*jetb1(1)*soft1(0))/12. - 
     -  (Pi**2*jeta1(0)*soft1(1))/12. + 
     -  zeta3*jeta1(1)*soft1(1) - 
     -  (Pi**2*jetb1(0)*soft1(1))/12. + 
     -  zeta3*jetb1(1)*soft1(1) + soft2(0)
     
      
      full(-1)=  hard(2) + (hard(1)*jeta1(-1))/2. + jeta2(-1)/4. + 
     -  (hard(1)*jetb1(-1))/2. + 
     -  (jeta1(-1)*jetb1(-1))/4. - 
     -  (Pi**2*jeta1(0)*jetb1(0))/24. + 
     -  (zeta3*jeta1(1)*jetb1(0))/4. + 
     -  (zeta3*jeta1(0)*jetb1(1))/4. - 
     -  (Pi**4*jeta1(1)*jetb1(1))/1440. + jetb2(-1)/4. + 
     -  hard(1)*soft1(-1) + (jeta1(-1)*soft1(-1))/2. + 
     -  (jetb1(-1)*soft1(-1))/2. - 
     -  (Pi**2*jeta1(0)*soft1(0))/12. + 
     -  (zeta3*jeta1(1)*soft1(0))/2. - 
     -  (Pi**2*jetb1(0)*soft1(0))/12. + 
     -  (zeta3*jetb1(1)*soft1(0))/2. + 
     -  (zeta3*jeta1(0)*soft1(1))/2. - 
     -  (Pi**4*jeta1(1)*soft1(1))/720. + 
     -  (zeta3*jetb1(0)*soft1(1))/2. - 
     -  (Pi**4*jetb1(1)*soft1(1))/720. + soft2(-1)
      
       assemble_dec=assemble_dec+ason2pi**2*(
     & +full(-1)
     & +full(0)*Lt1cut
     & +full(1)*Lt1cut**2/two
     & +full(2)*Lt1cut**3/three
     &     +full(3)*Lt1cut**4/four)*beama0*beamb0

      endif

      return
      end
