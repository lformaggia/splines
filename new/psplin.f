C=NAME psplin
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/new/psplin.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:42 $
C
C=TYPE subroutine
C
C
C=AUTHOR Luca Formaggia
C
C=PURPOSE general
C
C=KEYWORDS parametric_cubic_splines Parametric_curves
C
C=BLOCK ABSTRACT
C
c  It evaluates a cubic parametric spline passing through
c  the knots Q. The tangents are evaluated and stored in
c  T. Arrays a,c,b, are help array for the tridiagonal solver
C
C=END ABSTRACT
C
C
C=BLOCK USAGE
c
c  call psplin(ndimn,n,q,t,cs,len,ispty,a,b,c)
C
C   INPUT     DIMENSION     DESCRIPTION
c
c   ndimn                   n. of dimension
c   n                       n. of knots
c   q         ndimn,n       knots position vector
c   ispty     2             end knots condition (see VARIABLES DESCRIPTION below)
C   t         ndimn,n       in case of imposed tangent at the end points (ispty=2)
c                           t(*,1) or/and t(*,n) will contain the tangent value
c
C   OUTPUT    DIMENSION     DESCRIPTION
c
c   t         ndimn,n       tangent vector at knots (t  = dx/dc, where c is the
c                           'chord scaled' parameter.(see VARIABLES DESCRIPTION below)
c   cs        n-1           arcs chord length
c   len       n-1           arcs length
c
c   HELP ARRAYS
c
c   a         n             used by the tridiagonal solver
c   b         n               ..
c   c         n               ..
C
C   VARIABLES DESCRIPTION
c
c   ispty contains an indicator used to determine which condition to
c   apply at the spline curve end points. The allowed values are:
c
c    ispty =  0  natural spline (d2r/ds2 = 0)
c    ispty =  1  tangent imposed with direction given by the last two
c                spline point (see also routine GTIMP)
c    ispty =  2  tangent imposed, with associated value contained into t at the
c                corresponding position
c
c   The array t contains the 1st detivative at the knots with respect to a parameter
c   scaled with the chord length, t=x,c. It means that the derivative at the end of an arc,
c   with respect to the standard parameter u (whose variation between the two
c   end of the arc is unitary) can be computed as x,u = t*c , where c is the
c   arc chord length.
c   As a consequence x,u is not continous between arcs, while x,c it is.
c   In fact x,c is a better approximation of the tangent versor T = x,u/ ||x,u|| than x,u
c   utself.
c   For more details please  consult the USER'S GUIDE.
C
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine psplin(ndimn,n,q,t,cs,len,
     1                  ispty,a,b,c)
      parameter(eps=0.02,mdimn=8)
c
      real      q(ndimn,*),t(ndimn,*),cs(*)
      real      len(*),timp(2*mdimn)
      real      a(*),b(*),c(*)
      integer   ispty(2)
c
c evaluate chord length
c
      call cholen(ndimn,n,q,cs)
c
c get imposed tangents at ends
c
      if(ispty(1).eq.1)call gtimp(ndimn,timp,q,n,1)
      if(ispty(2).eq.1)call gtimp(ndimn,timp,q,n,2)
      if(ispty(1).eq.2)then
        call putimp(ndimn,timp,1,t,1)
      endif
      if(ispty(2).eq.2)then
        call putimp(ndimn,timp,2,t,n)
      endif
c
c evaluate tangents
c
      call evtan(ndimn,n,q,cs,a,b,c,t,ispty,timp)
c
c evaluate sectors length with accuracy eps*chord
c
      call slen(ndimn,n,q,t,cs,len,eps)
c
c evaluate curvature at knots
c
      return
      end
C
C=END SOURCE
