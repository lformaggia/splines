C=NAME curvps
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/curvps.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:35 $
C=TYPE real function
C=PURPOSE general
C
C=KEYWORDS parametric_surfaces
C=BLOCK ABSTRACT
C It evaluates the curvature on a surface along a given direction on the
C parametric plane
c
c   K = (vDv)/(vGv)
c
c where v is the direction on the parametric plane and G and D are the first and
C second fundamental matrices
C
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
C  xk-=curvps(ru,rv,ruv,ruu,rvv,xnor,v)
c
c
C  INPUT       DIMENSION      DESCRIPTION
C
C  ru         3             r,u
C  rv         3             r,v
c  ruu        3             r,uu
c  rvv        3             r,vv
c  ruv        3             r,uv
c
C OUTPUT      DIMENSION     DESCRIPTION
C
c  xnor       3             SURFACE NORMAL (NORMALISED)
c  curvps                   curvature
C
C=END USAGE
C
C=BLOCK SOURCE
C
       real*8 function curvps(ru,rv,ruv,ruu,rvv,xnor,v)
c
c the function evaluates the curvature of a parametric surface
c along the direction v. it return also the surface normal direction
c xnor. ru,rv,ruv,ruu,rvv are the surface 1st and 2nd derivatives at
c the given point
c
c       curvature = ( v d v) / (v g v)
c
c where d=2nd fundamental matrix and g= 1st fundamental matrix
c
      implicit none
      real*8 anor,d1,d2,d3,g1,g2,g3,ru(3),ruv(3),ruu(3),rv(3)
      real*8 rvv(3),v(2),v1,v2,xnor(3)
c
      xnor(1) = ru(2)*rv(3) - rv(2)*ru(3)
      xnor(2) = ru(3)*rv(1) - rv(3)*ru(1)
      xnor(3) = ru(1)*rv(2) - rv(1)*ru(2)
      anor = sqrt(xnor(1)*xnor(1)+xnor(2)*xnor(2)+xnor(3)*xnor(3))
      anor = 1.d0/anor
      xnor(1) = xnor(1)*anor
      xnor(2) = xnor(2)*anor
      xnor(3) = xnor(3)*anor
c
c     |d1  d2|       |g1 g2|
c  d =|      |    g =|     |
c     |d2  d3|       |g2 g3|
c
      d1 = ruu(1)*xnor(1)+ruu(2)*xnor(2)+ruu(3)*xnor(3)
      d2 = ruv(1)*xnor(1)+ruv(2)*xnor(2)+ruv(3)*xnor(3)
      d3 = rvv(1)*xnor(1)+rvv(2)*xnor(2)+rvv(3)*xnor(3)
      g1 = ru(1)*ru(1) +  ru(2)*ru(2) +  ru(3)*ru(3)
      g2 = ru(1)*rv(1) +  ru(2)*rv(2) +  ru(3)*rv(3)
      g3 = rv(1)*rv(1) +  rv(2)*rv(2) +  rv(3)*rv(3)
      v1 = v(1)
      v2 = v(2)
      curvps=(v1*v1*d1+2*v1*v2*d2+v2*v2*d3)/
     1       (v1*v1*g1+2*v1*v2*g2+v2*v2*g3)
      return
      end
C
C=END SOURCE
