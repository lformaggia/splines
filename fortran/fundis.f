C=NAME fundis
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/fundis.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:36 $
C=TYPE real function
C=PURPOSE local_to xmin1d
C=KEYWORDS Parametric_polynomial_curves 1D_Minimization
C=BLOCK ABSTRACT
c
c  The routine computes the distance between given point xp and the vector
c  position at a given point on a parametric curve, r(u)
C
C=END ABSTRACT
C
C=BLOCK USAGE
C
c     real  function fundis(u,xln1d,isn,v,dx,xp,r)
C
C
c   INPUT     DIMENSION       DESCRIPTION
c
c u                         parametric coordinate of point on curve
c n                         n. of knots on the curve (n. of arcs +1).
c xln1d       3,*           curve coefficent in standard format.
c isn         n-1           array of pointers to xln1d.
c xp          3             coordinates of given point.
c
c   OUTPUT    DIMENSION       DESCRIPTION
c
c fundis                    distance
c dx                        derivative d(||r(u)-xp||)/du
c r           3             position vector r(u)
c
C=END USAGE
C
C=BLOCK SOURCE
C
      real*8 function fundis(ug,xln1d,isn,n,xd,xp,r)
      implicit none
      integer*8 is,isn(*),ista,n
      real*8 h,r(3),rd(3),s,ug,u,v1,v2,v3,xd,xln1d(*),xp(3),xxx(3),zeps
      parameter(zeps=1.d-32)
c
c get arc and u value corresponding to ug
c
      is = min(max(1,int(ug)),n-1)
      u  = max(0.d0,min(ug-dble(is),1.d0))
      ista = isn(is)
      call evps1d(xln1d(ista),u,r,rd,xxx,xxx,s,2)
      v1 = r(1) - xp(1)
      v2 = r(2) - xp(2)
      v3 = r(3) - xp(3)
      h    = v1*v1 + v2*v2 + v3*v3
      fundis = sqrt(h)
      xd     = 2.d0*(v1*rd(1) + v2*rd(2) + v3*rd(3))
      xd     = xd/max(fundis,zeps)
      return
      end
C
C=END SOURCE
