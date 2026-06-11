C=NAME getp2
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/getp2.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:37 $
C=TYPE subroutine
C=PURPOSE general
C=KEYWORDS parametric_cubic_splines
C=BLOCK ABSTRACT
C
C It computes the position vector and various derivatives of a cubic spline
c at a given point in the parameter space
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
C call getp2(ndimn,n,q,t,cs,is,u,x,xu,xuu,xuuu,s,in)
C
C
C  INPUT       DIMENSION      DESCRIPTION
c
c ndimn                     N. of dimensions (IT MUST BE equal to 2 or 3 !!)
c n                         n. of knots
c q           ndimn,n       Knots coordinates
c cs          n             Arc chord length
C t           ndimn,n       Tangent at knots
c is                        Number of arc ( 1<= is <=n-1)
c u                         Parametric coordinate, local to arc is (0.<=u<=1.)
c in                        Switch :
c                            =0 -> xu will contain the derivative x,u
c                           !=0 -> xu will contain the tangent versor
c                                  (||xu|| =1)
c
C OUTPUT      DIMENSION     DESCRIPTION
c
c x           ndimn         position vector at local parametric coordinate u
c                           (arc number =is)
c xu          ndimn         x,u 1st derivative with respect to parameter u at the
c                           point (or tangent versor if in !=0)
c s                         ||x,u|| modulus of x,u
c xuu         ndimn         x,uu 2nd derivative
c xuuu        ndimn         x,uuu 3rd derivative
c
c
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine getp2(ndimn,n,q,t,cs,is,u,x,xd,xdd,xddd,s,in)
c
c  gives the values of the coordinate and various derivatives
c  for segmentis of the cubic spline at location u
c
c  x    =  point location         xd  = 1st derivative = dx/du
c  xdd  =  d2x/du2                xddd= d3x/du3
c  s    = | xd |                  in = normalization switch:
c                                  /=0 -> xd = dx/ds
      implicit none
      integer*8 i,id,in,ir,is,n,ndimn
      real*8 a1,a2,a3,clen,cs(*),q(ndimn,*),r12,s,s1,t(ndimn,*)
      real*8 t1,t2,u,x(ndimn),xd(ndimn),xdd(ndimn),xddd(ndimn)
c
      s  = 0.d0
      ir = is+1
      clen = cs(is)
      do 10 id=1,ndimn
         r12 = q(id,ir-1) - q(id,ir )
         t1  = t(id,ir-1)
         t2  = t(id,ir  )
         a1  =  2.d0*r12 + clen*(t1+t2)
         a2  = -3.d0*r12 - clen*(2.d0*t1+t2)
         a3  =  clen*t1
         x(id)=
     1         ( ( a1                        *u +
     1             a2                       )*u +
     1             a3                       )*u +
     1              q(id,ir-1)
         xd(id)= ( 3.d0*a1                   *u +
     1             2.d0*a2                  )*u +
     1               a3
         xdd(id) =   6.d0*a1*u + 2.d0*a2
         xddd(id)=   6.d0*a1
         s       =   s + xd(id)**2
 10   continue
      s  = sqrt(s)
      if (in.ne.0)then
        s1=1.d0/s
        do 30 i=1,ndimn
30      xd(i)  = xd(i)*s1
      endif
      return
      end
C
C=END SOURCE
