C=NAME pspl
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/pspl.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:39 $
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
c  It may use an arbitrary knot sequence
C
C=END ABSTRACT
C
C
C=BLOCK USAGE
c
c  call pspl(ndimn,n,q,t,us,ispty,a,b,c)
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
c   us        n-1           knot sequence: i.e. the interval between the value
c                           of the parameter between two adiacent knots
c                           It replaces cs
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
c    ispty =  0  natural spline (d2r/ds2 = 0).
c    ispty =  1  tangent imposed with direction given by the last two
c                spline point (see also routine GTIMP).
c    ispty =  2  tangent imposed, with associated value contained into t at the
c                corresponding position.
c    ispty =  3  Bessel end conditions: d2r/ds2 computed by a parabolic
c                interpolation.
c    ispty =  4  Quadratic end condition: d2r/ds2 at the end point put
c                equal to the adjacent one.
c
c   The array t contains the 1st detivative at the knots with respect to the parameter
c   whose knot sequence is specified into us. However, the routines which computes
c   the curve quantities refers to the standard parameter u with uniform knot sequence.
c   It means that the derivative at the end of an arc,
c   with respect to the standard parameter u (whose variation between the two
c   end of the arc is unitary) can be computed as x,u = t*us , where c is the
c   knot sequence.
c   As a consequence x,u is not continous between arcs, while x,c is.
c   For more details please  consult the USER'S GUIDE, Release 1.1
C
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine pspl(ndimn,n,q,t,us,ispty,a,b,c)
      implicit none
      integer*8 mdimn,n,ndimn,ispty(2)
      real*8 eps
      parameter(eps=0.02d0,mdimn=8)
      real*8 q(ndimn,*),t(ndimn,*),us(*),timp(2*mdimn)
      real*8 a(*),b(*),c(*)
c
c get imposed tangents at ends
c
      if(ispty(1).eq.1)call gtimp(ndimn,timp,q,n,1)
      if(ispty(2).eq.1)call gtimp(ndimn,timp,q,n,2)
      if(ispty(1).eq.2)then
        call putimp(ndimn,timp,1_8,t,1_8)
      endif
      if(ispty(2).eq.2)then
        call putimp(ndimn,timp,2_8,t,n)
      endif
c
c evaluate tangents
c
      call evtan(ndimn,n,q,us,a,b,c,t,ispty,timp)
c
      return
      end
C
C=END SOURCE
