C=NAME getiso
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/getiso.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:37 $
C=TYPE subroutine
C=AUTHOR L.formaggia
C=PURPOSE General
C=KEYWORDS Parametric_polynomial_curves Parametric_polynomial_surfaces
C
C=BLOCK ABSTRACT
C
c it finds the coefficients of a curve corresponding to a given iso-parametric
c line
C
C=END ABSTRACT
C=BLOCK USAGE
c
c call getiso(keypa,aps,n,m,iso,value,isn,xln1d,narc,marc,msel,ito,ierr)
c
C  INPUT      DIMENSION     DESCRIPTION
c
C  aps        3,*           array of coefficents
C  keypa      *             array of pointers to aps
C  n                        n. of knots along u-coordinate
c  m                        n. of knots along v-coordinate
c  iso        character*1   iso-line (u or v)
c  value                    value of the relevant parameter (global coordinate value)
c  nsseg                    number of curves (it is increased by one into the routine)
c  marc                     max. number of arcs
c  msel                     max. dimension in the array of ceofficents xln1d
C  ito                      i/o number associated to a log error file (6 if you
c                           want the output on the standard output)
c
C  OUTPUT     DIMENSION     DESCRIPTION
c
c  isn        marc          array of pointers to xln1d
c  xln1d      3,msel        array of coefficents
c  narc                     n. of archs +1 (number of knots) on curve
c  ierr                     error condition
c
c  * DIAGNOSTICS *
c
c  ierr = 0   no error
c         1    too many arcs on curve : ABORT  action: increase marc
c         2    too many coefficents   : ABORT  action: increase msel
c         3   1+2
c
c  * NOTES *
c
c if ierr /=0 an error message is sent to the unit associated to ito
c=END USAGE
C
C=BLOCK SOURCE
C
      subroutine getiso(keypa,aps,n,m,iso,value,isn,xln1d,
     1                  narc,marc,msel,ito,ierr)
c
      real aps(3,*),x1(3),x2(3),xln1d(3,msel),coef(3)
c
c get isoparametric curve from a surface
c
      integer keypa(*),ndu,ndv,i,j,n,m,kpa,isn(marc)
      character*1 iso
      apatch (id,ii,jj) =aps(id,kpa+(jj-1)*ndu+ii)
C
      ierr =0
      iarc  =1
      ixx =0
c
      if(iso.eq.'v'.or.iso.eq.'V')then
        narc =n
        j    = max(0,min(int(value),m-1))
        v    = value -j
c
        do 10 i=1,n-1
          if(i.le.marc)isn(i) = iarc
          npa = j + (i-1)*(m-1)
          kpa = keypa(npa)
          ndu = aps(1,kpa)
          ndv = aps(2,kpa)
          call evsurg(aps(1,kpa),0.,v,ndu,ndv,x1,xx,xx,xx,xx,xx,1)
          call evsurg(aps(1,kpa),1.,v,ndu,ndv,x2,xx,xx,xx,xx,xx,1)
          xl=sqrt((x2(1)-x1(1))**2+(x2(2)-x1(2))**2+
     1                             (x1(3)-x2(3))**2)
          xln1d(1,iarc)=ndu
          xln1d(2,iarc)=xl
          iarc         = iarc + 1
c
          do 11  ii = 1,ndu
             do 12 id=1,3
               coef(id) =0.
12           continue
             do 13 jj = ndv,1,-1
             do 13 id = 1,3
               coef(id) = coef(id)*v + apatch(id,ii,jj)
13           continue
             if(iarc.gt.msel)then
                ixx =iarc
                iarc=iarc +1
                go to 11
             endif
             xln1d(1,iarc) = coef(1)
             xln1d(2,iarc) = coef(2)
             xln1d(3,iarc) = coef(3)
             iarc          = iarc + 1
11        continue
10      continue
      else
        narc = m
        i    = max(0,min(int(value),n-1))
        u    = value -i
c
        do 20 j=1,m-1
          if(j.le.marc)isn(j)= iarc
          npa = j + (i-1)*(m-1)
          kpa = keypa(npa)
          ndu = aps(1,kpa)
          ndv = aps(2,kpa)
          call evsurg(aps(1,kpa),u,0.,ndu,ndv,x1,xx,xx,xx,xx,xx,1)
          call evsurg(aps(1,kpa),u,1.,ndu,ndv,x2,xx,xx,xx,xx,xx,1)
          xl=sqrt((x2(1)-x1(1))**2+(x2(2)-x1(2))**2+
     1                             (x1(3)-x2(3))**2)
          xln1d(1,iarc)=ndu
          xln1d(2,iarc)=xl
          iarc         = iarc + 1
c
          do 21  jj = 1,ndv
             do 22 id=1,3
               coef(id) =0.
22           continue
             do 23 ii = ndu,1,-1
             do 23 id = 1,3
               coef(id) = coef(id)*v + apatch(id,ii,jj)
23           continue
             if(iarc.gt.msel)then
                ixx =iarc
                iarc=iarc +1
                go to 21
             endif
             xln1d(1,iarc) = coef(1)
             xln1d(2,iarc) = coef(2)
             xln1d(3,iarc) = coef(3)
             iarc          = iarc + 1
21        continue
20      continue
      endif
c
c error handling
c
      if(narc.gt.marc)then
        write(ito,*)' ERROR IN GETISO : 001'
        write(ito,*)' TOO MANY ARCHS ON CURVE '
        write(ito,*)' increase MARC to at least  ',iarc
        ierr =1
      endif

       if(ixx.ne.0)then
         write(ito,*)' ERROR IN GETISO: 002'
         write(ito,*)' TOO MANY COEFFICENTS ON CURVE '
         write(ito,*)' increase MSEL to at least  ',ixx
         ierr = ierr + 2
       endif
c
       return
       end
C
C=END SOURCE
c
