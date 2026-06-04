C=NAME pspli2
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/new/pspli2.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:41 $
C
C=TYPE subroutine
C
C
C=AUTHOR Luca Formaggia
C
C=PURPOSE general
C
C=KEYWORDS parametric_cubic_splines
C
C=BLOCK ABSTRACT
C
c  It evaluates the parameters for cubic spline interpolation
c  It differs from the similar routine psplin because now q contains the
c  value of the variable to be interpolated by a cubic spline.
C
C=END ABSTRACT
C
C
C=BLOCK USAGE
c
c  call psplin2(ndimn,n,q,t,cs,ispty,a,b,c,niv,rxy)
C
C   INPUT     DIMENSION     DESCRIPTION
c
c   ndimn                   n. of components of vector q
c   n                       n. of knots
c   q         ndimn,n       value of variables at knots
c   ispty     2             end knots condition (see VARIABLES DESCRIPTION below)
C   t         ndimn,n       in case of imposed tangent at the end points (ispty=2)
c                           t(*,1) or/and t(*,n) will contain the imposed value of the
c                           derivative at the two ends of the curve
c   niv                     number of spatial dimensions
c   rxy       niv,n         knots coordinates
c
C   OUTPUT    DIMENSION     DESCRIPTION
c
c   t         ndimn,n       derivative at knots (t  = dq/dc, where c is the
c                           'chord scaled' parameter.(see VARIABLES DESCRIPTION below)
c   cs        n-1           arcs chord length
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
c    ispty =  0  natural  (d2q/ds2 = 0)
c    ispty =  1  1st derivative imposed with direction given by the last two
c                value of the quantity q
c    ispty =  2  1st derivative imposed, with associated value contained into t at the
c                corresponding position
c
c   The array t contains the 1st detivative at the knots with respect to a parameter
c   scaled with the chord length, t=q,c. It means that the derivative at the end of an arc,
c   with respect to the standard parameter u (whose variation between the two
c   end of the arc is unitary) can be computed as q,u = t*c , where c is the
c   arc chord length.
c   As a consequence q,u is not continous between arcs, while q,c is.
c   In fact q,c is a better approximation of the derivative with respect to arc length
c   than q,u.
c
c   For more details please  consult the USER'S GUIDE.
C
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine pspli2(ndimn,n,q,t,cs,
     1                  ispty,a,b,c,niv,rxy)
      implicit none
      integer*8 mdimn,n,ndimn,niv,ispty(2)
      parameter(mdimn=8)
      real*8 q(ndimn,*),t(ndimn,*),cs(*),timp(2*mdimn)
      real*8 a(*),b(*),c(*),rxy(niv,*)
c
c  new version rxy are the point coordinates q is the
c  values to interpolate
c
c evaluate chord length
c
      call cholen(niv,n,rxy,cs)
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
      call evtan(ndimn,n,q,cs,a,b,c,t,ispty,timp)
c
c evaluate curvature at knots
c
      return
      end
C
C=END SOURCE
