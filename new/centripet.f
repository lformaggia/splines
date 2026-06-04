c=NAME centripet
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/new/centripet.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:42 $
C=TYPE subroutine
C=PURPOSE psplin
C=KEYWORDS  parametric_cubic_splines
C=BLOCK ABSTRACT
C
C It computes the deltas of the knot sequence corresponding
c to a centripetal parametrization.
C
C=END ABSTRACT
C
C=BLOCK USAGE
C
C       call centripet(n,cs,delta)
c
c INPUT
c
c VARIABLES    DIMENSION    DESCRIPTION
c
c n                         N. of knots
c cs           n-1          chord length
c
c OUTPUT
c
c VARIABLES    DIMENSION    DESCRIPTION
c 
C delta        n-1          delta's corresponding to centripetal
c                           parametrisation 
c
c=BLOCK BIBLIOGRAPHY 
c
c G.Farin, Curves and Surfaces for Computer Aided Geometric Design
c A practical Guide, Academic press, 1988, pag 110.
c
C=END BIBLIOGRAPHY     
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine centripet(n,cs,delta)
c
      implicit none
      integer*8 i,n
      real*8 cs(n),delta(n)
c
c evaluates centripetal parametrisation
c
      do 10 i=1,n-1
         delta(i)=sqrt(max(cs(i),0.d0))
10    continue
      return
      end
C
C=END SOURCE
