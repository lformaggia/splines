c=NAME correctt
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/new/correctt.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:42 $
C=TYPE subroutine
C=PURPOSE psplin
C=KEYWORDS  parametric_cubic_splines
C=BLOCK ABSTRACT
C
C It corrects the tangent at the end of the spline
c It should be called only when the tangent value is imposed
c Its use is necessary in order to use a value for the imposed
c tangent which do not depend strongly from the chosen 
c parametrisation. This is accomplished by multiplying the
c tangent by the ration chord/delta
c
C=END ABSTRACT
C
C=BLOCK USAGE
C
C       call correctt(timp,ndimn,delta,chord,n,iend)
c
c INPUT
c
c VARIABLES    DIMENSION    DESCRIPTION
c
c n                         N. of knots
c ndimn                     N. of dimensions
c chord        n-1          chord length
c delta        n-1          forward differences of knot sequence
c                           chosen for the parametrisation
c iend                      End of the spline 1-> first 2-> second
c timp        ndimn,2       Unscaled tangents components
c
c OUTPUT
c
c VARIABLES    DIMENSION    DESCRIPTION
c 
c timp        ndimn,2       Scaled tangents components
c
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine correctt(timp,ndimn,delta,chord,n,iend)
c
      implicit none
      integer*8 i,igo,iend,n,ndimn
      real*8 chord(n),delta(n),den,scal,timp(ndimn,2)
      igo =iend
      if(iend.eq.2)igo=n-1
c
      do 10 i=1,ndimn
c
         den = abs(delta(igo))
c
         if(den.eq.0.d0)then
           scal=1.d0
         else
           scal = chord(igo)/den
         endif
c     
         timp(i,iend) = timp(i,iend)*scal
c
10    continue
      return
      end
C
C=END SOURCE
