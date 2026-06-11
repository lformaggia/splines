C=NAME gders1
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/gders1.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:36 $
C=TYPE subroutine
C=PURPOSE general
C=KEYWORDS parametric_bicubic_spline
C=BLOCK ABSTRACT
C
C It computes the various 1st anbd 2nd derivatives at a point r = r(u,v)
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
C call gders1(ndimn,ru,rv,ruv,ruu,rvv,apatch,u,v)
C
C
C  INPUT       DIMENSION      DESCRIPTION
c
c ndimn                     n. of dimensions
c apatch      ndimn,4,4     array of coefficents of the bicubic expansion
c                             r = apatch(i,j)*(v**(j-1))*(u**(i-1)) (i,j=1,...,4)
c u,v                       parametric coordinates
c
C OUTPUT      DIMENSION     DESCRIPTION
c
c ru          ndimn         r,u
c rv          ndimn         r,v
c ruv         ndimn         r,uv
c rvv         ndimn         r,vv
c ruu         ndimn         r,uu
c
C=END USAGE
C
C=BLOCK SOURCE
C
       subroutine gders1(ndimn,ru,rv,ruv,ruu,rvv,apatch,u,v)
       implicit none
c
c this routine gets the various derivatives at the point u,v
c of the bicubic patch  r = uav
c
       integer*8 id,ndimn
       real*8 apatch(ndimn,4,4),d2,d3,d4,ru(*),ruu(*),ruv(*),rv(*)
       real*8 rvv(*),s1,s2,s3,s4,u,v,z1,z2,z3,z4
c
      do 10 id =1,ndimn
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
         s1 = (3*apatch(id,1,4) *v+
     1         2*apatch(id,1,3))*v+
     1           apatch(id,1,2)
c
         s2 = (3*apatch(id,2,4) *v+
     1         2*apatch(id,2,3))*v+
     1           apatch(id,2,2)
c
         s3 = (3*apatch(id,3,4) *v+
     1         2*apatch(id,3,3))*v+
     1           apatch(id,3,2)
c
         s4 = (3*apatch(id,4,4) *v+
     1         2*apatch(id,4,3))*v+
     1           apatch(id,4,2)
c
         z1 =  6*apatch(id,1,4)*v+
     1         2*apatch(id,1,3)
c
         z2 =  6*apatch(id,2,4)*v+
     1         2*apatch(id,2,3)
c
         z3 =  6*apatch(id,3,4)*v+
     1         2*apatch(id,3,3)
c
         z4 =  6*apatch(id,4,4)*v+
     1         2*apatch(id,4,3)
c
        rv(id) = ((s4*u+s3)*u+s2)*u+s1
        ru(id) = (3*d4*u+2*d3)*u + d2
        ruv(id)= (3*s4*u+2*s3)*u + s2
        ruu(id)=  6*d4*u +2*d3
        rvv(id) = ((z4*u+z3)*u+z2)*u+z1
10    continue
      return
      end
C
C=END SOURCE
