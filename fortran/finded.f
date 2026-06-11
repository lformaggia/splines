C=NAME finded
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/finded.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:36 $
C=TYPE subroutine
C=AUTHOR L.formaggia
C=PURPOSE General
C=KEYWORDS Parametric_polynomial_curves Parametric_polynomial_surfaces
C
C=BLOCK ABSTRACT
C
c it finds the coefficients of a curve corresponding to a surface edge
c
c iedg: edge number according to the following layout
c
c          ......3........
C         !              !
C         !              !
C         4              2
C   (v)^  !              !
C     j!  !......1.......!
C      !
C      !____>
C          i (u)
C
C
C
C=END ABSTRACT
C=BLOCK USAGE
c
c call finded(aps,keypa,n,m,iedg,nsseg,isn,xln1d,narc,marc,msel,ito,ierr)
c
C  INPUT      DIMENSION     DESCRIPTION
c
C  aps        3,*           array of coefficents
C  keypa      *             array of pointers to aps
C  n                        n. of knots along u-coordinate
c  m                        n. of knots along v-coordinate
c  iedg                     edge number (see abstract)
c  nsseg                    number of curves (it is increased by one into the routine)
c  marc                     max. number of arcs
c  msel                     max. dimension in the array of ceofficents xln1d
C  ito                      i/o number associated to a log error file (6 if you
c                           want the output on the standard output)
c
C  OUTPUT     DIMENSION     DESCRIPTION
c  nsseg                    number of curves (it is increased by one into the routine)
c  isn         marc          array of pointers to xln1d
c  xln1d      3,msel        array of coefficents
c  narc                     n. of archs +1 (number of knots) on curve
c  ierr                     error condition
c
c  * DIAGNOSTICS *
c
c  ierr = 0   no error
c         1    too many arcs on curve : ABORT  action: increase marc
c         2    too many coefficents   : ABORT  action: increase msel
c
c  * NOTES *
c
c if ierr /=0 an error message is sent to the unit associated to ito
c=END USAGE
C
C=BLOCK SOURCE
C
      subroutine  finded(aps,keypa,n,m,iedg,nsseg,isn,xln1d,narc,
     1                   marc,msel,ito,ierr)
      implicit none
C
C FINDS THE COEFFICIENTS OF A CURVE CORRESPONDING TO A SURFACE EDGE
C
C IEDG: EDGE NUMBER ACCORDING TO
C
C          ......3........
C         !              !
C         !              !
C         4              2
C         !              !
C     j!  !......1.......!
C      !
C      !____
C          i
C
C
      integer*8 i,iarc,ico,id,iedg,ien,ierr,ist,ito,ixx,j,jen,jst
      integer*8 k,keypa(*),kk,kpa,marc,m,msel,n,na,narc,nd,ndu,ndv
      integer*8 nsseg,isn(marc)
      real*8 a1,a2,a3,aps(3,*),u1,u2,v1,v2,x1(3),x2(3),xl,xln1d(3,msel)
      real*8 xx(3)
C
      ierr =0
      ico  =1
C
      if(iedg.eq.1.or.iedg.eq.3)then
        na=n
        ist=1
        ien=n-1
        u1  =0.d0
        u2  =1.d0
        if (iedg.eq.1)then
            jst=1
            jen=1
            v1 =0.d0
            v2 =0.d0
        else
            jst=m-1
            jen=m-1
            v1 =1.d0
            v2 =1.d0
        endif
      else if(iedg.eq.2.or.iedg.eq.4)then
        na=m
        jst=1
        jen=m-1
        v1 =0.d0
        v2 =1.d0
        if (iedg.eq.2)then
            ist=n-1
            ien=n-1
            u1 = 1.d0
            u2 = 1.d0
        else
            ist=1
            ien=1
            u1 =0.d0
            u2 =0.d0
        endif
      endif
      if(na.gt.marc)then
        write(ito,*)' ERROR IN FINDED : 002'
        write(ito,*)' TOO MANY ARCHS ON CURVE N. ',nsseg
        write(ito,*)' increase MARC to at least  ',na
        ierr =2
        return
      endif
      narc = na
      iarc = 0
      ixx  = 0
      do 10 i=ist,ien
      do 10 j=jst,jen
         iarc = iarc +1
         k = j+(i-1)*(m-1)
         kpa=keypa(k)
         ndu = aps(1,kpa)
         ndv = aps(2,kpa)
         if(iedg.eq.1.or.iedg.eq.3)then
             nd=ndu
         else
             nd=ndv
         endif
c
         call evsurg(aps(1,kpa),u1,v1,ndu,ndv,x1,xx,xx,xx,xx,xx,1)
         call evsurg(aps(1,kpa),u2,v2,ndu,ndv,x2,xx,xx,xx,xx,xx,1)
         xl=sqrt((x2(1)-x1(1))**2+(x2(2)-x1(2))**2+
     1                   (x1(3)-x2(3))**2)
         isn(iarc)=ico
         xln1d(1,ico) = dble(nd)
         xln1d(2,ico) = xl
         ico          = ico +1
         if(iedg.eq.1)then
               do 70 k=1,ndu
                 a1=aps(1,kpa+(1-1)*ndu+k)
                 a2=aps(2,kpa+(1-1)*ndu+k)
                 a3=aps(3,kpa+(1-1)*ndu+k)
                 if (ico.gt.msel)then
                   ixx = ico
                   ico = ico+1
                   go to 70
                 endif
                 xln1d(1,ico) = a1
                 xln1d(2,ico) = a2
                 xln1d(3,ico) = a3
                 ico          = ico +1
70             continue
         else if(iedg.eq.3)then
               do 71 k=1,ndu
                 a1=aps(1,kpa+(1-1)*ndu+k)
                 a2=aps(2,kpa+(1-1)*ndu+k)
                 a3=aps(3,kpa+(1-1)*ndu+k)
                 do 72 kk=2,ndv
                   a1=a1+aps(1,kpa+(kk-1)*ndu+k)
                   a2=a2+aps(2,kpa+(kk-1)*ndu+k)
                   a3=a3+aps(3,kpa+(kk-1)*ndu+k)
72               continue
                 if (ico.gt.msel)then
                   ixx = ico
                   ico = ico+1
                   go to 71
                 endif
                 xln1d(1,ico) = a1
                 xln1d(2,ico) = a2
                 xln1d(3,ico) = a3
                 ico          = ico +1
71              continue
         else if(iedg.eq.4)then
               do 80 k=1,ndv
                 a1=aps(1,kpa+(k-1)*ndu+1)
                 a2=aps(2,kpa+(k-1)*ndu+1)
                 a3=aps(3,kpa+(k-1)*ndu+1)
                 if (ico.gt.msel)then
                   ixx = ico
                   ico = ico+1
                   go to 80
                 endif
                 xln1d(1,ico) = a1
                 xln1d(2,ico) = a2
                 xln1d(3,ico) = a3
                 ico          = ico +1
80             continue
         else if(iedg.eq.2)then
               do 81 k=1,ndv
                 a1=aps(1,kpa+(k-1)*ndu+1)
                 a2=aps(2,kpa+(k-1)*ndu+1)
                 a3=aps(3,kpa+(k-1)*ndu+1)
                 do 82 kk=2,ndu
                   a1=a1+aps(1,kpa+(k-1)*ndu+kk)
                   a2=a2+aps(2,kpa+(k-1)*ndu+kk)
                   a3=a3+aps(3,kpa+(k-1)*ndu+kk)
82               continue
                 if (ico.gt.msel)then
                   ixx = ico
                   ico = ico+1
                   go to 81
                 endif
                 xln1d(1,ico) = a1
                 xln1d(2,ico) = a2
                 xln1d(3,ico) = a3
                 ico          = ico +1
81              continue
         endif
10     continue
c
       if(ixx.ne.0)then
         write(ito,*)' ERROR IN FINDED : 002'
         write(ito,*)' TOO MANY COEFFICENTS ON CURVE N. ',nsseg
         write(ito,*)' increase MSEL to at least  ',ixx
         ierr =3
       endif
c
       return
       end
C
C=END SOURCE
