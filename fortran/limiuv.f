C=NAME limiuv
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/limiuv.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:38 $
C=TYPE subroutine
C=PURPOSE Local_to xminli
C=KEYWORDS Parametric_polynomial_surfaces Line_search_Minimization
C=BLOCK ABSTRACT
c
c This routine is used to find the interval into which the minimum is searched.
c The 2D minimization is limited into the interval defined by u10,u20,v10,v20. Given
c the coordinate of a stating point (ugs,vgs) and a direction xi, the routine compute
c the value of bl1 and bl2 such that
c
c       u10<=ugs + x*xi(1)<=u20   and   v10<=vgs + x*xi(2)<=v20
c
c when                   bl1<=x<=bl2
c
C=END ABSTRACT
C=BLOCK USAGE
C
c call limiuv(ugs,vgs,xi,u10,u20,v10,v20,bl1,bl2)
C
C
c   INPUT     DIMENSION       DESCRIPTION
c
c ugs,vgs                   starting value of the global parametric coordinates
c xi          2             search direction
c u10,u20                   limits of interval (u-coords.)
c v10,v20                   limits of interval (v-coords.)
c
c   OUTPUT     DIMENSION       DESCRIPTION
c
c bl1,bl2                   interval of x so that
c                           u10<=ugs + x*xi(1)<=u20 , v10<=vgs + x*xi(2)<=v20
c                           when bl1<=x<=bl2
c
C=END USAGE
C
C=BLOCK SOURCE
C
       subroutine limiuv(ugs,vgs,xi,u10,u20,v10,v20,bl1,bl2)
      implicit none
      real*8 b1,b2,b3,b4,bl1,bl2,t,u10,u20,ugs,v10,v20,vgs,xi(2),zero
      parameter(zero=1.d-20)
c
c     find b1l, b2l such that  u10< ugs + b1l*xi(1) < u20
c     and                      v10<= vgs + b
c
      if(abs(xi(1)).gt.zero)then
        b1 = (u10 - ugs)/xi(1)
        b2 = (u20 - ugs)/xi(1)
        t  = max(b1,b2)
        b2 = b1+ b2 -t
        b1 = t
      else
        b1 =  1.d+10
        b2 = -1.d+10
      endif
c
      if(abs(xi(2)).gt.zero)then
        b3 = (v10 - vgs)/xi(2)
        b4 = (v20 - vgs)/xi(2)
        t  = max(b3,b4)
        b4 = b3+ b4 -t
        b3 = t
      else
        b3 =  1.d+10
        b4 = -1.d+10
      endif
c
      bl1 = max(b2,b4)
      bl2 = min(b1,b3)
      return
      end
C
C=END SOURCE
