C=NAME xmaxpa
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/xmaxpa.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:40 $
C=TYPE subroutine
C=PURPOSE local
C=KEYWORDS parametric_Polynimial_surfaces
C=BLOCK ABSTRACT
C
c It computes the coordinates of the knots (i.e. the points at the corner of each patch) \
c and the min-max dimensions of a surface
c
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
c call xmaxpa(n,m,keypa,aps,xmax,np,c3d)
c
C  INPUT       DIMENSION      DESCRIPTION
c
c n                         number of patches -1 along u direction (equal to the
c                           number of knots)
c m                         number of patches -1 along v direction (equal to the
c                           number of knots)
c keypa       *             array of pointers to aps
c aps         3,*           array of coefficients for the surface
c
C OUTPUT      DIMENSION     DESCRIPTION
c
c np                        number of knots ( n*m)
c c3d         3,np          knots coordinates
c xmax        3,2           (*,1) = min.  coordinates value
c                           (*,2) = max.  coorinates value
c
c NOTE
c
c The points are ordered along the u-coordinate. It means that if we assign the
c couple of numbers (i,j) to a knot according to its position in the lattice, where
c point (i,j) correpond to c3d(*,k) where k= (j-1)*n + i. Viceversa
c  i  = mod(k-1,n)+1  j = 1+(k-1)/n (INTEGER DIVISION!)
c
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine xmaxpa(n,m,keypa,aps,xmax,np,c3d)
c
c find max coordinates of the point on a surface
c
      real aps(3,*),xmax(3,2),xp(3),c3d(3,*)
      integer keypa(*)
c
      np =0
      do 10 id=1,3
        xmax(id,1)= 1.e+32
        xmax(id,2)= -1.e+32
10     continue
c
       do 500 j=1,m-1
         kkv=0
         if(j.eq.m-1)kkv=1
         do 530 kv=0,kkv
           v= float(kv)
           do 510 i=1,n-1
             kku =0
             if(i.eq.n-1)kku=1
             npa =(i-1)*(m-1)+j
             kpa =keypa(npa)
             nu = aps(1,kpa)
             nv = aps(2,kpa)
             do 510 ku=0,kku
               u=float(ku)
               k = k +1
               call evsurg(aps(1,kpa),u,v,nu,nv,xp,
     1                               xx,xx,xx,xx,xx,1)
               xmax(1,1)=min(xmax(1,1),xp(1))
               xmax(1,2)=max(xmax(1,2),xp(1))
               xmax(2,1)=min(xmax(2,1),xp(2))
               xmax(2,2)=max(xmax(2,2),xp(2))
               xmax(3,1)=min(xmax(3,1),xp(3))
               xmax(3,2)=max(xmax(3,2),xp(3))
               np = np+1
               c3d(1,np)=xp(1)
               c3d(2,np)=xp(2)
               c3d(3,np)=xp(3)
c
510        continue
530      continue
500   continue
      return
      end
C
C=END SOURCE
C
C=END DECK
c

