c=NAME uniform
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/new/uniform.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:41 $
C=TYPE subroutine
C=PURPOSE psplin
C=KEYWORDS  parametric_cubic_splines
C=BLOCK ABSTRACT
C
CIt assemble into the array delta the delta's corresponding
c to uniform parametrisation
C
C=END ABSTRACT
C
C=BLOCK USAGE
C
C       call uniform(n,delta)
C
C INPUT
C
C VARIABLES       DIMENSION    DESCRIPTION
C
c  n                           N. of knots of the cardinal spline
c
C OUTPUT
C
C VARIABLES       DIMENSION    DESCRIPTION
c 
c delta            n-1         delta(i)=first forward difference
c                              of knot sequence. In the case of 
c                              uniform parametrization all deltas
c                              are unitary.
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine uniform(n,delta)
c
      implicit none
      integer*8 i,n
      real*8 delta(n)
c
      do 10 i=1,n-1
       delta(i) =1.d0
10    continue
      return
      end
C
C=END SOURCE
