C=NAME gpsur
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/gpsur.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:37 $
C=TYPE subroutine
C=PURPOSE general
C=KEYWORDS parametric_bicubic_splines
C=BLOCK ABSTRACT
C
C It computes the position vector on a bicubic surface
C
C=END ABSTRACT
C=BLOCK USAGE
C
C call gpsur(ndimn,u,v,r,apatch)
C
C  INPUT       DIMENSION      DESCRIPTION
c
c ndimn                     N. of dimensions
c u,v                       local parametric coordinates on the patch
c apatch      ndimn,4,4     coefficent matrix as computed by using EVAPA2
c
C OUTPUT      DIMENSION     DESCRIPTION
c
c r           3             position vector
c
c
c NOTES: This routine uses the matrix apatch containing the coefficents of the bicubic
c        patch, as it is produced by the routine EVAPA2. This array is NOT in the standard
c        format (the header is missing). If ndimn=3 it is suggested then to use the routine evsurg,
c        together with routine EVAPA3. EXAMPLE:
c
c evaluation of the position vector oa the point (0.5,0.6) on patch (i=2,j=4) of
c a surface formed by 15*10 patches (16*11 knots)
c
c with gpsur:
c
c
c   parameter(ndimn=3)
c   real apatch(ndimn,4,4)
c       ...
c   call evapa2(16,11,ndimn,tanret,choret,coor,apatch,2,4)
c   call gpsur(ndimn,0.5,0.6,r,apatch)
c       ....
c
c with evsurg (suggested mode if ndimn=3):
c
c   paramater(ndimn=3)
c   real apatch(ndimn,1+4*4)
c      ...
c   call evapa3(16,11,ndimn,tanret,choret,coor,apatch,2,4)
c   call evsurg(apatch,0.5,0.6,4,4,r,ru,rv,ruv,ruu,rvv,1)
c      ....
c
c
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine gpsur(ndimn,u,v,r,apatch)
c
c get a point of coordinate parameters u,v of the
c bicubic patch
c
      implicit none
      integer*8 id,ndimn
      real*8 apatch(ndimn,4,4),d1,d2,d3,d4,r(ndimn),u,v
c
      do 10 id=1,ndimn
         d1 = ((apatch(id,1,4) *v+
     1          apatch(id,1,3))*v+
     1          apatch(id,1,2))*v+apatch(id,1,1)
c
         d2 = ((apatch(id,2,4) *v+
     1          apatch(id,2,3))*v+
     1          apatch(id,2,2))*v+apatch(id,2,1)
c
         d3 = ((apatch(id,3,4) *v+
     1          apatch(id,3,3))*v+
     1          apatch(id,3,2))*v+apatch(id,3,1)
c
         d4 = ((apatch(id,4,4) *v+
     1          apatch(id,4,3))*v+
     1          apatch(id,4,2))*v+apatch(id,4,1)
c
         r(id) = ((d4*u+d3)*u+d2)*u+d1
10     continue
       return
       end
C
C=END SOURCE
