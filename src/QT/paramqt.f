      double precision Pi,Zeta2,Zeta3,Zeta4,Zeta5
      parameter (Pi=3.14159265358979323846d0)
      parameter (Zeta2=1.6449340668482264d0)
      parameter (Zeta3=1.2020569031595942d0)
      parameter (Zeta4=1.0823232337111382d0)
      parameter (Zeta5=1.0369277551433699d0)
      double precision CA,CF,TF
      parameter (CA=3d0)
      parameter (CF=1.3333333333333333d0)
      parameter (TF=0.5d0)

      double precision nf,sumQsq
      double precision beta0,beta1,beta2
      double precision G0,G1,G2
      double precision d1,d2,d3
      double precision gH0qq,gH1qq,gH2qq,cH1qq,cH2qq,cH3qq
      double precision gH0qq2,gH1qq2,gH2qq2
      double precision gH0gg,gH1gg,gH2gg,cH1gg,cH2gg,cH3gg
      double precision gH0gg2,gH1gg2,gH2gg2
      common /qtparam/ nf,sumQsq,
     &                 beta0,beta1,beta2,
     &                 G0,G1,G2,
     &                 d1,d2,d3,
     &                 gH0qq,gH1qq,gH2qq,cH1qq,cH2qq,cH3qq,
     &                 gH0gg,gH1gg,gH2gg,cH1gg,cH2gg,cH3gg,
     &                 gH0qq2,gH1qq2,gH2qq2,gH0gg2,gH1gg2,gH2gg2
!$omp threadprivate(/qtparam/)

      double precision PDFbuf1(-5:5,0:1)
      double precision PDFbuf2(-5:5,0:1)
      common /pdfdata/ PDFbuf1,PDFbuf2
!$omp threadprivate(/pdfdata/)
      
      double precision Li2,Li3,Li4
      double precision PlusDist
      external Li2,Li3,PlusDist
      
c      double precision YLP0gg,YLP0qiqi,YLP0qg,YLP0gq
c      double precision YLP1qiqi,YLP1qiqj,YLP1qiqbi,YLP1qg,YLP1gg,YLP1gq
c      double precision I1gg,I1qiqi,I1qg,I1gq
c      double precision YLI2qiqi,YLI2qiqj,YLI2qiqbi,YLI2qg,YLI2gg,YLI2gq

c      double precision YLP0qiqiYLP0qiqi,YLP0qgYLP0gq,YLP0qgYLP0gg,YLP0qiqiYLP0qg
c      double precision YLP0ggYLP0gg,YLP0gqYLP0qg,YLP0gqYLP0qiqi,YLP0ggYLP0gq
c      double precision I1qiqiYLP0qiqi,I1qgYLP0gq,I1qgYLP0gg,I1qiqiYLP0qg
c      double precision I1ggYLP0gg,I1gqYLP0qg,I1gqYLP0qiqi,I1ggYLP0gq
c      double precision Ii1gq,Ii1gg
