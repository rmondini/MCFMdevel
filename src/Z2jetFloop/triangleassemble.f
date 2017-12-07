      subroutine triangleassemble(p1,p2,p3,p4,p5,p6,za,zb,mtsq,Ccoeff)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'qqbggintnames.f'
C     c12_34=1,c12_3=2,c12_4=3,c3_124=4,c4_123=5,c3_4=6
      integer p1,p2,p3,p4,p5,p6,b3,b4
      real(dp):: mtsq,mt2
      complex(dp):: Ccoeff(6,2,2),coeff(2,2),ans(2,2)
      mt2=mtsq


      ans(1,1)=dcmplx(-9.63230767649925768D-002,-0.48813970226006931d0     )
      ans(1,2)=dcmplx(-0.58034951522681344d0,-0.75898912763715420d0)
      ans(2,1)=dcmplx(0.28624655060222559d0,-0.60146328543455840d0)
      ans(2,2)=dcmplx(-0.25986822666137871d0,0.21478754583785398d0)

      call qqbggAxtri12x3x456(p1,p2,p3,p4,p5,p6,b3,b4,za,zb,mt2,coeff)
      Ccoeff(c12_3,1,1)=+coeff(1,1)
      Ccoeff(c12_3,1,2)=+coeff(1,2)
      Ccoeff(c12_3,2,1)=+coeff(2,1)
      Ccoeff(c12_3,2,2)=+coeff(2,2)
      write(6,*) '12_3,mm',Ccoeff(c12_3,1,1)/ans(1,1)
      write(6,*) '12_3,mp',Ccoeff(c12_3,1,2)/ans(1,2)
      write(6,*) '12_3,pm',Ccoeff(c12_3,2,1)/ans(2,1)
      write(6,*) '12_3,pp',Ccoeff(c12_3,2,2)/ans(2,2)
      write(6,*)
      stop

      call qqbggAxtri12x3x456(p1,p2,p4,p3,p5,p6,b3,b4,za,zb,mt2,coeff)
      Ccoeff(c12_4,1,1)=-coeff(1,1)
      Ccoeff(c12_4,1,2)=-coeff(2,1)
      Ccoeff(c12_4,2,1)=-coeff(1,2)
      Ccoeff(c12_4,2,2)=-coeff(2,2)
      write(6,*) '12_4,mm',Ccoeff(c12_4,1,1)
      write(6,*) '12_4,mp',Ccoeff(c12_4,1,2)
      write(6,*) '12_4,pm',Ccoeff(c12_4,2,1)
      write(6,*) '12_4,pp',Ccoeff(c12_4,2,2)
      write(6,*)

      call qqbggAxtri3x4x1256(p1,p2,p3,p4,p5,p6,za,zb,coeff)
      Ccoeff(c3_4,1,1)=+coeff(1,1)
      Ccoeff(c3_4,1,2)=+coeff(1,2)
      Ccoeff(c3_4,2,1)=+coeff(2,1)
      Ccoeff(c3_4,2,2)=+coeff(2,2)
      write(6,*) '3_4,mm',Ccoeff(c3_4,1,1)
      write(6,*) '3_4,mp',Ccoeff(c3_4,1,2)
      write(6,*) '3_4,pm',Ccoeff(c3_4,2,1)
      write(6,*) '3_4,pp',Ccoeff(c3_4,2,2)
      write(6,*)
      pause

c      ans(1,1)=dcmplx(0.18356436548127525d0,0.93438580157035200d0)
c      ans(1,2)=dcmplx(1.0699796935851986d0,1.5222343216637184d0)
c      ans(2,1)=dcmplx(-0.53212292719039400d0,1.1376175241037432d0)
c      ans(2,2)=dcmplx(0.42182580015634852d0,-0.38980581816940651d0)


      b3=6
      b4=5
      write(6,*) 'mt2',mt2
      write(6,*) 'b3,b4',b3,b4
      call qqbggAxtri123x4x56(p1,p2,p3,p4,p5,p6,b3,b4,za,zb,mt2,coeff)
      Ccoeff(c4_123,1,1)=+coeff(1,1)
      Ccoeff(c4_123,1,2)=+coeff(1,2)
      Ccoeff(c4_123,2,1)=+coeff(2,1)
      Ccoeff(c4_123,2,2)=+coeff(2,2)
      write(6,*) '4_123,mm',Ccoeff(c4_123,1,1)
      write(6,*) '4_123,mp',Ccoeff(c4_123,1,2)
      write(6,*) '4_123,pm',Ccoeff(c4_123,2,1)
      write(6,*) '4_123,pp',Ccoeff(c4_123,2,2)
      write(6,*)


      call qqbggAxtri123x4x56(p1,p2,p4,p3,p5,p6,b3,b4,za,zb,mtsq,coeff)
      Ccoeff(c3_124,1,1)=-coeff(1,1)
      Ccoeff(c3_124,1,2)=-coeff(2,1)
      Ccoeff(c3_124,2,1)=-coeff(1,2)
      Ccoeff(c3_124,2,2)=-coeff(2,2)
      write(6,*) '3_124,mm',Ccoeff(c3_124,1,1)
      write(6,*) '3_124,mp',Ccoeff(c3_124,1,2)
      write(6,*) '3_124,pm',Ccoeff(c3_124,2,1)
      write(6,*) '3_124,pp',Ccoeff(c3_124,2,2)
      write(6,*)

      return
      end
