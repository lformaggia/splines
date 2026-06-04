C=NAME reseca
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/reseca.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:40 $
C=TYPE subroutine
C=PURPOSE general Local_to readca
C=KEYWORDS Parametric_Polynomial_Curves Data_Base_Interface CATIA_neutral_file
C=BLOCK ABSTRACT
C
C   The routine reads a SINGLE curve from a catia neutral file and stores the relevant
c   information into the Database
C
C
C=END ABSTRACT
C=BLOCK USAGE
c
c         call reseca(ioc,n,xln1d,isn,msel,marc,il1,il2,ug1,ug2,icati,ito,ierr,linen)
c
c  INPUT      DIMENSION     DESCRIPTION
c
c msel                      max. dimension of  file of coefficents for curve (xln1d)
c marc                      max. dimension of  file of keys for curve (isn)= max. n. of arcs
c ioc                       i/o unit number associated to catia file
c ito                       i/o unit number associated to log/error file
c linen                     number of last line read from CATIA neutral file (it is incremented by 1 for each
c                           line read)
c
c  OUTPUT     DIMENSION     DESCRIPTION
c
c n                         n. of knots (n. of arcs+1)
c xln1d       3,msel        array of coefficents
c isn         marc          keys to array of  coefficents
c il1,il2                   local  numbering of the surfaces adjacent to the curve
c ug1,ug2                   global coordinates of definition curve end points
c icati                     curve ID number from CATIA (normally taken as global numbering)
c linen                     number of last line read from CATIA neutral file (it is incremented by 1 for each
c                           line read)
c ierr                      error indicator (see DIAGNOSTICS)
c
c *NOTES and DIAGNOSTICS*
c
c ierr =0     no error
c       1     too many arcs.
c              memory overflow in isn: ABORT.
c              action: increase marc
c       2     too many coefficents.
c              memory overflow in xln1d: ABORT.
c              action: increase msel
c       3     1 + 2
c
c  The unit associated to ito is used for the logging of mesages. In standard FORTRAN
c  if you want the messages on the standard output you should put ito=6
c=END USAGE
C
C=BLOCK SOURCE
C
      subroutine reseca(ioc,narc,xln1d,isn,msel,marc,
     1                  il1,il2,ug1,ug2,icati,ito,ierr,linen)
      implicit none
C
C READ CURVE COEFFS FROM CATIA GEOMETRY DEFINITION FILE
C
c
c     ioc = catia file input cahnnel
c     narc =n .of arcs +1
c     xln1d(3,*) = curve coefficents array
c     isn(*)     = pointer to xln1d
c     msel = max. total n. of coefficents
c     marc = total n. of arcs
c     il1,il2= curve - surface links
c     ug1,ug2 = global parametric coordinates of ens points
c     icati = CATIA's curve ID number
c     ito   = standard output channel n.
c     ierr  = error condition n. (0=no error)
c
c
      integer*8 i,icati,icount,ico,icx,ierr,ikk,il(2),il1,il2
      integer*8 io_dummy,ioc,is1,is2,isn(marc),ito,ix,j,linen,marc
      integer*8 msel,narc,ndeg,nli,nlink
      real*8 u1,u2,ug1,ug2,x,x1(3),x2(3),xl,xln1d(3,*),y,z,xx(3)
      character*80 line,cha*1,lin2*80
c
c
      read(ioc,'(A)')line
      linen = linen+1
c
      ierr = 0
      ico =1
c
      if(line(2:5).ne.'NARC')then
         write(ito,*)' ERROR IN RESECA :001'
         write(ito,*)' KEYWORD NOT CORRESPONDING TO A CATIA SEGMENT'
         write(ito,*)' READING ABORTED'
         ierr =1
         return
      endif
c
c read catia curve id number (if it exists)
c
      icati =0
      icount =0
      do 99 i=14,34
       if(line(i:i+3).eq.'*CRV')then
         do 195 j=i+4,34
           cha = line(j:j)
           if(cha.eq.' ')goto 196
           read(cha,'(i1)',err=197)ikk
           icount = icount*10 + ikk
195      continue
       endif
99    continue
196   icati=icount
197   continue
      if(icati.eq.0)then
        write(ito,*) ' WARNING in RESECA'
        write(ito,*) ' CANNOT READ CATIA CURVE ID NUMBER'
      endif
c
c read curve links (if any)
c
      il(1) = 0
      il(2) = 0
      nlink =0
      read(line,'(T39,i2)',err=115)nlink
      go to 125
115   write(ito,*) ' WARNING in RESECA'
      write(ito,*) ' CANNOT READ NUMBER OF LINKS FOR CURVE'
125   continue
c
c read links number
c
      nli =0
      if(nlink.ne.0)then
         do 301 i=42,78
         if(line(i:i+2).eq.'SUR')then
           ikk =0
           do 303 j=i+3,80
            if(line(j:j).ne.' ')then
              cha=line(j:j)
              read(cha,'(i1)',err=432)ix
              ikk = ikk*10 + ix
            else
              nli = nli+1
              il(nli)=ikk
              if(nli.lt.nlink)then
                  go to 301
              else
                  goto 304
              endif
            endif
303        continue
        endif
301     continue
      endif
304   continue
      il1 = il(1)
      il2 = il(2)
c
      if(nli.eq.2)go to 433
432   write(ito,*) ' WARNING in RESECA'
      write(ito,*) ' CANNOT READ 2  LINKS FOR CURVE',icati
c
433   continue
c
      if(nlink.ne.nli)then
        write(ito,*) ' WARNING in RESECA'
        write(ito,*) ' N. OF LINKS DECLARED NOT EQUALS N. FOUND '
        write(ito,*) ' IN CURVE',icati
      endif
c
c
      read(line,'(T7,I6)')narc
      if(narc.gt.marc)then
        write(ito,*)' ERROR IN RESECA :002'
        write(ito,*)' NUMBER OF ARCHS EXCEEDED MAXIMUM'
        write(ito,*)' INCREASE MARC TO AT LEAST',narc
        ierr=2
      endif
      narc = narc+1
      read(ioc,678)is1,u1,is2,u2
      linen=linen+1
      ug1 = is1 + u1
      ug2 = is2 + u2
678   format(t5,i3,f8.4,t21,i3,f8.4)
      icx =ico
      do 10 i=1,narc-1
         isn(i)=ico
         read(ioc,'(T19,I6)')ndeg
         linen = linen+1
         xln1d(1,ico)=dble(ndeg)
         ico = ico +1
         read(ioc,'(A)')line
         linen = linen+1
         read(ioc,'(A)')line
         linen = linen +1
         do 20 j=1,ndeg
           read(ioc,'(9X,3(f16.10,1x))')x,y,z
c           read(ioc,'(9X,3F16.5)')x,y,z
           linen = linen+1
           if(ico.gt.msel)go to 20
           xln1d(1,ico)=x
           xln1d(2,ico)=y
           xln1d(3,ico)=z
           ico = ico + 1
c
20       continue
C
C EVALUATE ARCH CHORD LENGTH
C
         if(ico.gt.msel)go to 10
         call evps1d(xln1d(1,icx),0.d0,x1,xx,xx,xx,x,1)
         call evps1d(xln1d(1,icx),1.d0,x2,xx,xx,xx,x,1)
         xl = sqrt((x1(1)-x2(1))**2+(x1(2)-x2(2))**2+(x1(3)-x2(3))**2)
         xln1d(2,icx)=xl
10    continue
      if(ico.gt.msel)then
        write(ito,*)' ERROR IN RESECA :003'
        write(ito,*)' TOO MANY COEFFICENTS IN CATIA CURVE'
        write(ito,*)' INCRESE MSEL TO AT LEAST',ico
        ierr =3
      endif
      return
      end
C
C=END SOURCE
