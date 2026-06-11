c=NAME trcoe
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/trcoe.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:40 $
C=TYPE subroutine
C=AUTHOR Luca Formaggia
C=PURPOSE local_to psplin psplin2
C=KEYWORDS  parametric_cubic_splines Parametric_curves
C=BLOCK ABSTRACT
C
c Get the tridiagonal system, and stores it into a,b and c
c
C=END ABSTRACT
C
C=BLOCK USAGE
c
c      call  TRCOE(ndimn,n,q,cs,a,b,c,t,ispty,timp)
c
c INPUT       DIMENSION     DESCRIPTi0N
c ndimn                     N. of dimensions
c n                         n. of knots
c q           ndimn,n       knots coordinates or value of the variables to
c                           be interpolated at knots
c cs          n             Arc chord length
c ispty       2             End conditions:
c                           0 -> natural spline 
c                           1,2 -> tangent imposed
c                           3 -> Bessel end conditions
c                           4 -> Quadratic end conditions
c                           5 -> Not-aKnot conditions
c timp        ndimn,2       value of imposed tangent (used only if
c                           corresponding ispty !=0)
c
c  OUTPUT     DIMENSION     DESCRIPTION
c
C t           ndimn,n       right hand side of tri-diagonal linear sistem
c a           n             coefficents used by tridiagonal solver (lower diagonal)
c b           n              coefficents used by tridiagonal solver (principal diagonal)
c c           n             coefficents used by tridiagonal solver (upper diagonal)
c
c  NOTES
c
c The tridiagonal sistem set up by this routine is:
c
c    |   b1  c1               |     t1
c    |   a2  b2  c2           |     t2
c    |       a3  b3  c3       | X =  .
c    |           ....         |      .
c    |                 an  bn |     tn
c
c where
c        a(i)    = cs(i)                       !
c        b(i)    = 2*(cs(i)+cs(i-1)            !
c        c(i)    = cs(i-1)                     ! for 2 <= i <= n
c        ti(*) = (3./(a(i)+c(i)))*             !
c                 (c(i)**2*(q(*,i+1)-q(*,i))+  !
c                  a(i)**2*(q(*,i)-q(*,i-1)))  !
c
c
c while the values at i=1 and i=n are determined by the boundary conditions imposed
c
c natural spline :
c
c  b(1) = 2*cs(1)    c(1)=cs(1  )  t(*,1) = 3*(q(*,2)-q(*,1)  ) (1st end)
c  b(n) = 2*cs(n-1)  a(n)=cs(n-1)  t(*,n) = 3*(q(*,n)-q(*,n-1)) (2nd end)
c
c tangent imposed :
c
c  b(1) = 1.0    c(1)= 0.0  t(*,1) = timp(*,1) (1st end)
c  b(n) = 1.0    a(n)= 0.0  t(*,n) = timp(*,n) (2nd end)
c
c Bessel : The tangent vector at the end is calcultated from a parabola interpolating
c          the last 3 points
c 
c Quadratic: The second derivative at the last 2 points is set equal
c
c Not-a-Knot: The last 2 polinomial segments mesrge into a single cubic
c
c
c For details on Bessel, not-a-knot and quadratic conditions, please refer to Reference
c
c
c   The routine is able to handle degenerate situations where some
c   arc chord lengths are zero.
c
c WARNING: The only degenerate situation succesfully tested so far is the one where
c          all the chord length are zero, i.e all the knots coincides with a point.
c
c
c=BLOCK REFERENCE 
c                 Faux and Pratt, Computational Geometry for Design and Manufacture,
c                 Ellis-Horwood Ltd., 1979, pp. 223-226.
c
c                 G.Farin, Curves and Surfaces for Computer Aided Geometric Design
c                 A practical Guide, Academic press, 1988, pp. 106-108.
c
C=END REFERENCE
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine trcoe(ndimn,n,q,cs,a,b,c,t,ispty,timp)
      implicit none
      integer*8 i,i0,i1,i2,id,ii1,ii2,ispty(2),n,ndimn
      real*8 a(*),a1,b(*),beta,c(*),c1,cs(*),q(ndimn,*),t(ndimn,*)
      real*8 timp(ndimn,*),zero
      parameter (zero=1.d-10)
c
      ii1   = ispty(1)
      ii2   = ispty(2)
c
      if(ii1.lt.0)ii1=0
      if(ii2.lt.0)ii2=0
c
      if(n.le.2.or.
     1   abs(cs(1)).lt.zero.or.abs(cs(n-1)).lt.zero)then
        if(ii1.gt.2)ii1=0
        if(ii2.gt.2)ii2=0
      endif
c
c
c Set tri-di system for interior points
c
      do 10 i=2,n-1
        a1   = cs(i  )
        c1   = cs(i-1)
        if (abs(a1).lt.zero)a1=zero
        if (abs(c1).lt.zero)c1=zero
        a(i) = a1
        b(i) = 2.d0*(a1+c1)
        c(i) = c1
        do 10 id=1,ndimn
         t(id,i) = (3.d0/(a1*c1))*
     1                  (c1*c1*(q(id,i+1)-q(id,i  ))+
     1                   a1*a1*(q(id,i  )-q(id,i-1)) )
10    continue
c
      if(ii1.eq.1.or.ii1.eq.2)then
c
c tangent imposed
c
         do 30 id = 1,ndimn
            t(id   ,1)=timp(id,1)
30       continue
         b(1)=1
         c(1)=0
         a(1)=0
      else if(ii1.eq.0)then
c
c natural spline
c
         do 40 id =1,ndimn
            t(id   ,1)=3.d0*(q(id   ,2) -q(id,1))
40       continue
         b(1)=2.d0*cs(1)
         c(1)=cs(1)
         a(1)=0.d0
         if(b(1).lt.zero)b(1)=zero
      else if(ii1.eq.3)then
c
c Bessel
c
         beta=b(2)/2.d0
         i0 = 1
         i1 = 2
         i2 = 3
         do 33 id=1,ndimn
            t(id,1)=
     1      -q(id,i0)*(2*cs(i0)+cs(i1))/(beta*cs(i0))+
     1      q(id,i1)*beta/(cs(i0)*cs(i1))-
     1      q(id,i2)*(cs(i0)/(cs(i1)*beta))
33       continue
c
         a(1)=0.d0
         b(1)=1.d0
         c(1)=0.d0
c
      else if(ii1.eq.4)then
c
c quadratic
c
         i0 = 1
         i1 = 2
         do 34 id=1,ndimn
            t(id,1)=2.d0*(q(id,i1)-q(id,i0))/cs(i0)
34       continue
         a(1)=0.d0
         b(1)=1.d0
         c(1)=1.d0
c
      else if(ii1.eq.5)then
c
c not-a-knot
c
         i0 = 1
         i1 = 2
         i2 = 3
         beta=b(2)
         do 35 id=1,ndimn
            t(id,1)=
     1      (q(id,i2)-q(id,i1))*(cs(i0)*cs(i0))/cs(i1)+
     1      (q(id,i1)-q(id,i0))*(cs(i1)/cs(i0))*
     1          (3*cs(i0)+2*cs(i1))
35       continue
         a(1)=0.d0
         b(1)=cs(i1)*beta
         c(1)=beta*beta
      endif
c
      if(ii2.eq.1.or.ii2.eq.2)then
c
c tangent imposed
c
         do 32 id = 1,ndimn
            t(id   ,n)=timp(id,2)
32       continue
         b(n)=1.d0
         a(n)=0.d0
         c(n)=0.d0
      else if(ii2.eq.0)then
c
c natural spline
c
         do 41 id =1,ndimn
            t(id   ,n)=3.d0*(q(id   ,n) -q(id,n-1) )
41       continue
         b(n)=2.d0*cs(n-1)
         a(n)=cs(n-1)
         c(n)=0.d0
         if(b(n).lt.zero)b(n)=zero
c
      else if(ii2.eq.3)then
c
c Bessel
c
         i0 = n
         i1 = n-1
         i2 = n-2
         beta=b(n-1)/2.d0
         do 43 id=1,ndimn
            t(id,n)=
     1      q(id,i0)*(2*cs(i1)+cs(i2))/(beta*cs(i1))-
     1      q(id,i1)*beta/(cs(i1)*cs(i2))+
     1      q(id,i2)*(cs(i1)/(cs(i2)*beta))
43       continue
c
         a(n)=0.d0
         b(n)=1.d0
         c(n)=0.d0
c
      else if(ii2.eq.4)then
c
c quadratic
c
         i0 = n-1
         i1 = n
         do 44 id=1,ndimn
            t(id,n)=2.d0*(q(id,i1)-q(id,i0))/cs(i0)
44       continue
         a(n)=1.d0
         b(n)=1.d0
         c(1)=0.d0
c
      else if(ii2.eq.5)then
c
c not-a-knot
c
         i0 = n
         i1 = n-1
         i2 = n-2
         beta=b(n-1)
         do 45 id=1,ndimn
            t(id,n)=
     1   (q(id,i1)-q(id,i2))*(cs(i1)*cs(i1))/cs(i2)+
     1   (q(id,i0)-q(id,i1))*(cs(i2)/cs(i1))*
     1   (3.d0*cs(i1)+2.d0*cs(i2))
45       continue
         a(n)=beta*beta
         b(n)=cs(i2)*beta
         c(n)=0.d0
      endif
      return
      end
C
C=END SOURCE
