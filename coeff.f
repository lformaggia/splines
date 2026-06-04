C=NAME coeff
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/coeff.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:35 $
C=TYPE SUBROUTINE
C=PURPOSE  local_to slen psplin
C=KEYWORDS Parametric_cubic_splines
C=BLOCK ABSTRACT
C
C A secondary routine, associeted to subroutine SLEN which computes the length of an
C arc of cubic spline. The routine evaluates the coefficents of the polinomial
C (4th order) P + ||x,u||*2 = scoef(1) + scoef(2)*u+ ... + scoef(5)*u**4.
C=END ABSTRACT
C
C=BLOCK USAGE
c
C
C
C call COEFF(ndimn,n,q,t,cs,is,scoef)
C
C  INPUT       DIMENSION      DESCRIPTION
C
C  ndimn                    number of dimension
C  n                        number of knots
C  q          ndimn,n       knots coordinates
C  t          ndimn,n       tangent vector at knots
C  cs         ndimn,n-1     chord length
C  is                       arc number
C
C OUTPUT      DIMENSION     DESCRIPTION
C
C scoef       5             coefficents for the polynomial expansion of
C                           ||x,u|| on the arc.
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine coeff(ndimn,n,q,t,cs,is,scoef)
c
c *** this sub. computes the coeficients of the polynomial |r'|**2
c     for segment n. is (from point is-1 to is)
c
      implicit none
      integer*8 i,id,is,n,ndimn
      real*8 a1,a2,a3,a4,a5,f,len,p,pp,q(ndimn,*),qq,r12
      real*8 scoef(5),s,t(ndimn,*),z,cs(*)
c
c     do 10 i=2,n
      i  =is+1
      a1 = 0.d0
      a2 = 0.d0
      a3 = 0.d0
      a4 = 0.d0
      a5 = 0.d0
c
      do 1000 id=1,ndimn
      r12 = q(id,i-1)-q(id,i)
      len = cs(i-1)
      p   = len*t(id,i-1)
      z   = len*t(id,i  )
      pp = p*p
      qq = z*z
      s = 3.d0*( 2.d0*r12+   p+z)
      f = 2.d0*(-3.d0*r12-2.d0*p-z)
      a1 = a1+pp
      a2 = a2+2.d0*p*f
      a3 = a3+f*f+2.d0*p*s
      a4 = a4+2.d0*f*s
      a5 = a5+s*s
 1000 continue
c
c     ic=i-1
      scoef(1)=a1
      scoef(2)=a2
      scoef(3)=a3
      scoef(4)=a4
      scoef(5)=a5
c
      return
      end
C
C=END SOURCE
