C=NAME getta3
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/getta3.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:37 $
C=TYPE SUBROUTINE
C
C=AUTHOR  L. Formaggia
C=PURPOSE General
C
C=KEYWORDS Parametric_bicubic_splines, parametric_polynomial_surfaces
C
C=BLOCK ABSTRACT
c
C The subroutine interpolates a surface spline formed by tensor-product
c bicubic patches on a lattice of points (knots). It uses the method
c described in the reference. It allows for different specification
c of tangent boundary conditions.
c
c Modification: 1/4/90 : it allows for degenerate patches, i.e. all points in
c a iso-parametric line at one of the border of the surface in the parametric
c plane may map into a single point in the 'physical' space.
c
c (A border on the parametric plane is a curve at u=constant=1 or u=constant=n etc.etc.)
c
c
c=END ABSTRACT
C
c=BLOCK REFERENCE 
c                 Faux and Pratt, Computational Geometry for Design and Manufacture,
c                 Ellis-Horwood Ltd., 1979, pp. 223-226.
c
c                 G.Farin, Curves and Surfaces for Computer Aided Geometric Design
c                 A practical Guide, Academic press, 1988, pp. 106-108.
c
C=END REFERENCE
C
C=BLOCK USAGE
C
C      call getta3(ndimn,n,m,mm,coor,tanret,choret,q,cs,a,b,c,t,ispt3)
c
c   INPUT     DIMENSION       DESCRIPTION
c
c ndimn                     n. of dimensions
c n                         n. of knots along u
c m                         n. of knots along v
c mm                        max n. of points along one direction ( >= max(n,m) )
c coor        ndimn,*       knot coordinates. The coordinate corresponding to
c                           knot (i,j) in the lattice is coor(*,i+n*(j-1))
c ispt3       3,2           tangent end condition marker              :
c                           (1,*) -> condition on r,u on  iso-v boundary curve
c                           (2,*) -> condition on r,v on  iso-u boundary curves
c                           (3,*) -> condition on r,uv on the surface corners.
c                           (*,1) -> first end
c                           (*,2) -> second end
c
c                            ispt3 = -1 : Derivative evaluated using an internal
c                                         procedure (the same as +1)
c                            ispt3 =  0 : Natural spline , second derivative =0.
c                            ispt3 =  1 : Derivative imposed with direction given by the last two
c                                         spline point (see also routine GTIMP).
c                            ispt3 =  2 : Derivative imposed. The value is contained in the
c                                         corresponding position of array  tanret.
c                            ispt3 =  3 : Bessel end condition: the derivative is set equal
c                                         to the one corresponding to a parabola through the
c                                         last 3 data points 
c                            ispt3 =  4 : Quadratic end condition: the second derivative
c                                         at the end point is teken equal to the one at the
c                                         previous point.
c                            ispt3 =  5 : Not-a-Knot condition. The last two polinomial expansions
c                                         merge into a single cubic
c
c tanret      ndimn,*,n,m   value of imposed derivative
c                           (only if the corresponding value of ispt3 is equal to 1)
c
c
c
c   OUTPUT    DIMENSION     DESCRIPTION
c
c tanret      ndimn,3,n*m   evaluated tangents at knots:
c                                         (*,1,*) -> cu * r,u
c                                         (*,2,*) -> cv * r,v
c                                         (*,3,*) -> cuv* r,uv
c
c                           where cu,cv,cuv are the scaling factors which
c                            take into account the change of parametrization.
c                           The parametrization used for inter-patch derivative
c                           continuity is in fact based on chord length. For more
c                           details consult the User's manual
c
c choret      2,n,m         chord length of the arcs of the iso-parametric curves
c                           departing from node (i,j),relative to patch (i,j):
c                           (1,*,*) -> chord lenght on v-constant isolines
c                           (2,*,*) -> chord length on u-constant isolines
c                           [The actual memory occupied by choret is 2*(n-1)*(m-1).]
c
C
C  HELP ARRAYS
C
c q           mm
c cs          mm
c a           mm
c b           mm
c t           mm
c
c    Patch description:
c
c               SEE USER'S MANUAL
c
c
c
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine getta3(ndimn,n,m,mm,coor,
     1                  tanret,choret,
     1                  q,cs,a,b,c,t,ispt3)
c
c m,n dimension of the lattice, mm=max(m,n) .        n. of
c points in the lattice= m*n
c
      real tanret(ndimn,3,n,m),choret(2,n,m)
      real q(ndimn,mm),cs(mm),a(mm),b(mm),c(mm),t(ndimn,mm)
      real coor(ndimn,*),timp(6)
c
      integer ispty(2),ispt3(3,2)
      iretp(i,j)=(j-1)*n +i
c
c get r,v along the u=cost lines
c
      ispty(1)=1
      ispty(2)=1
      do 100 i=1,n
         do 110 j=1,m
         do 110 id=1,ndimn
             q(id,j) = coor(id,iretp(i,j))
110      continue
         call cholen(ndimn,m,q,cs)
         if(abs(ispt3(2,1)).eq.1)then
            ispty(1)=1
            call gtimp(ndimn,timp,q,m,1)
         else if(ispt3(2,1).eq.0)then
            ispty(1)=0
         else if(ispt3(2,1).eq.2)then
            ispty(1)=1
            call movt(ndimn,timp,1,tanret,2,i,1,n,m)
         else
            ispty(1)=ispt3(2,1)
         endif
         if(abs(ispt3(2,2)).eq.1)then
            ispty(2)=1
            call gtimp(ndimn,timp,q,m,2)
         else if(ispt3(2,2).eq.0)then
            ispty(2)=0
         else if(ispt3(2,2).eq.2)then
            ispty(2)=1
            call movt(ndimn,timp,2,tanret,2,i,m,n,m)
         else
            ispty(2)=ispt3(2,2)
         endif
         call evtan(ndimn,m,q,cs,a,b,c,t,ispty,timp)
         do 120 j=1,m
         do 120 id=1,ndimn
             tanret(id,2,i,j)=t(id,j)
120       continue
          do 130 j=1,m-1
             choret(2,i,j)=cs(j)
130       continue
100   continue
c
c interpolate r,u along v=cost lines
c
      do 200 j=1,m
         do 210 i=1,n
         do 210 id=1,ndimn
             q(id,i) = coor(id,iretp(i,j))
210      continue
         call cholen(ndimn,n,q,cs)
         if(abs(ispt3(1,1)).eq.1)then
            ispty(1)=1
            call gtimp(ndimn,timp,q,n,1)
         else if(ispt3(1,1).eq.0)then
            ispty(1)=0
         else if(ispt3(1,1).eq.2)then
            ispty(1)=1
            call movt(ndimn,timp,1,tanret,1,1,j,n,m)
         else
            ispty(1)=ispt3(1,1)
         endif
         if(abs(ispt3(1,2)).eq.1)then
            ispty(2)=1
            call gtimp(ndimn,timp,q,n,2)
         else if(ispt3(1,2).eq.0)then
            ispty(2)=0
         else if(ispt3(1,2).eq.2)then
            ispty(2)=1
            call movt(ndimn,timp,2,tanret,1,n,j,n,m)
         else
            ispty(2)=ispt3(1,2)
         endif
         call evtan(ndimn,n,q,cs,a,b,c,t,ispty,timp)
c
         do 220 i=1,n
         do 220 id=1,ndimn
             tanret(id,1,i,j)=t(id,i)
220       continue
          do 230 i=1,n-1
             choret(1,i,j)=cs(i)
230       continue
200   continue
c
c evaluate r,vu along the u lines j=1 and j=m
c
c
      do 400 k=1,2
        if(k.eq.1)j=1
        if(k.eq.2)j=m
        do 410 i=1,n
        do 410 id=1,ndimn
           q(id,i) = tanret(id,2,i,j)
410     continue
         if(abs(ispt3(3,1)).eq.1)then
            ispty(1)=1
            call gtimp(ndimn,timp,q,n,1)
         else if(ispt3(3,1).eq.0)then
            ispty(1)=0
         else if(ispt3(3,1).eq.2)then
            ispty(1)=1
            call movt(ndimn,timp,1,tanret,3,1,j,n,m)
         else
            ispty(1)=ispt3(3,1)
         endif
c
         if(abs(ispt3(3,2)).eq.1)then
            ispty(2)=1
            call gtimp(ndimn,timp,q,n,2)
         else if(ispt3(3,2).eq.0)then
            ispty(2)=0
         else if(ispt3(3,2).eq.2)then
            ispty(2)=1
            call movt(ndimn,timp,2,tanret,3,n,j,n,m)
         else
            ispty(2)=ispt3(3,2)
         endif
c
        do 420 i=1,n-1
           cs(i) = choret(1,i,j)
420     continue
c
        call evtan(ndimn,n,q,cs,a,b,c,t,ispty,timp)
        do 430 i=1,n
        do 430 id=1,ndimn
           tanret(id,3,i,j)=t(id,i)
430     continue
400   continue
c
c again: interpolate r,uv vorking along the v lines and
c using as boundary condition the tangent evaluated in the
c previous step
c
      ispty(1) = 1
      ispty(2) = 1
c
      do 500 i=1,n
c
        do 510 j=1,m
        do 510 id=1,ndimn
           q(id,j) = tanret(id,1,i,j)
510     continue
        do 520 j=1,m-1
           cs(j)   = choret(2,i,j)
520     continue
c
        do 522 id=1,ndimn
          timp(id)       =tanret(id,3,i,1)
          timp(id+ndimn) =tanret(id,3,i,m)
522     continue
c
        call evtan(ndimn,m,q,cs,a,b,c,t,ispty,timp)
c
        do 540 j=1,m
        do 540 id=1,ndimn
           tanret(id,3,i,j) = t(id,j)
540     continue
500   continue
      return
      end
C
C=END SOURCE
