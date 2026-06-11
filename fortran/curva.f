C=NAME curva
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/curva.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:35 $
C=TYPE subroutine
C=PURPOSE general
C=KEYWORDS parametric_cubic_splines
C=BLOCK ABSTRACT
C
C It gives an estimate of the curvature at the knots for a parametric
C cubic spline
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
C call CURVA(ndimn,n,q,t,cs,rad
C
C
C  INPUT       DIMENSION      DESCRIPTION
C
C  ndimn                    number of dimension
C  n                        number of knots
C  q          ndimn,n       knots coordinates
C  t          ndimn,n       tangent vector at knots
C  cs         ndimn,n-1     chord length
C
C OUTPUT      DIMENSION     DESCRIPTION
C
c rad         n-1           estimate of the curvature at the arcs
C
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine curva(ndimn,n,q,t,cs,rad)
c
c evaluates an estimate of the curvature at the knots
c
      implicit none
      integer*8 i,id,n,ndimn
      real*8 clen,cs(*),pp,q(ndimn,*),rad(*),t(ndimn,*),xk
c
      do 10 i=1,n-1
        xk=0.d0
        clen=cs(i)
        do 20 id=1,ndimn
           pp = 2.d0*(3.d0*(q(id,i+1)-q(id,i))-
     1               clen*(t(id,i+1)+2*t(id,i)))
           xk = xk + pp*pp
20      continue
        xk = sqrt(xk)
        rad(i) = xk/(clen*clen)
10    continue
      xk =0.d0
      clen=cs(n-1)
      do 30 id =1,ndimn
         pp = 6.d0*(q(id,n-1)-q(id,n))+
     1        2.d0*(t(id,n-1)+2.d0*t(id,n))*clen
         xk =xk +pp*pp
30    continue
      xk =sqrt(xk)
      rad(n) =xk/(clen*clen)
      return
      end
C
C=END SOURCE
