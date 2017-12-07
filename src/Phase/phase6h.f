      subroutine phase6h(r,p1,p2,p3,p4,p5,p6,p7,p8,wt)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'breit.f'
      include 'kprocess.f'
c********* generate phase space for 2-->5 process
c********* r(mxdim),p1(4),p2(4) are inputs 
c--------- incoming p1 and p2 reversed in sign from physical values 
c---- i.e. phase space for -p1-p2 --> p3+p4+p5+p6+p7
c---- with all 2 pi's (ie 1/(2*pi)^11)

      real(dp):: r(mxdim),wt,wt345678,wt34,wt5678,wt67,tmp,wt56,wt78
      real(dp):: p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),p7(4),p8(4)
      real(dp):: p5678(4),p12(4),p67(4),p34(4),p56(4),p78(4),s3min
      integer:: oldn3,oldn2
      real(dp), parameter:: wt0=1._dp/twopi**4, tiny=1e-6_dp

      p12(:)=-p1(:)-p2(:)

      oldn3=n3
      n3=1
      call phi1_2(r(1),r(2),r(3),r(4),p12,p5678,p34,wt345678,*99)
      n3=0
      call phi3m0(r(5),r(6),p34,p3,p4,wt34,*99)
      oldn2=n2
      n2=0
      n3=0
      call phi1_2(r(7),r(8),r(11),r(12),p5678,p56,p78,wt5678,*99)
      call phi3m0(r(13),r(14),p56,p5,p6,wt56,*99)
      call phi3m0(r(15),r(16),p78,p7,p8,wt78,*99)
      
      
      if(p5(4).ne.p5(4)) goto 99
      if(p6(4).ne.p6(4)) goto 99
      if(p7(4).ne.p7(4)) goto 99
      if(p8(4).ne.p8(4)) goto 99

      
      n2=oldn2
      n3=oldn3
      wt=wt0*wt345678*wt5678*wt56*wt78*wt34
  
      
      return
 99   continue
      wt=0._dp
      n3=oldn3
      n2=oldn2
      return
      end

