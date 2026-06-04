c=NAME surfb
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/surfb.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:40 $
C=TYPE subroutine
C=AUTHOR L.formaggia
C=PURPOSE General
C=KEYWORDS  Parametric_polynomial_surfaces Parametric_bicubic_splines
C
C=BLOCK ABSTRACT
C
c It builds a surface passing trough a lattice of points and stores it directly
c into the genaral array type. The surface is a natural surface ( first
c derivatives zero at the surface border and cross derivative zero at the surface
c corners)
C
C=END ABSTRACT
C=BLOCK USAGE
c
c call surfb(ifa,rest,nssur,n,m,mpatch,msul,ito,aps,keypa,ierr)
c
C  INPUT      DIMENSION     DESCRIPTION
c
c  ifa                      i/o unit connected with the file containing the coordinates
c                           of the knots. The file is read in the format
c                           READ(IFA,*)N,M
C                           READ(IFA,*)(((coord(id,i,j),id=1,3),i=1,n),j=1,m)
c                             the knots must form a lattice n*m, the i-index corresponds
c                           to the u-corrdinate, (consequently j corresponds to the
c                           v-coords)
c  nssur                    number of surfaces (NOT USED)
c  msul                     max. dimension in the array of coefficents aps
c  mpatch                   max. number of patches
C  ito                      i/o number associated to a log error file (6 if you
c                           want the output on the standard output)
c
C  OUTPUT     DIMENSION     DESCRIPTION
c
C  aps         3,msul        array of coefficents
C  keypa      mpatch        array of pointers to aps
C  n                        n. of knots along u-coordinate (read from file)
c  m                        n. of knots along v-coordinate (read from file)
c  ierr                     error condition
c
c  HELP ARRAY
c
c  rest       mx (see note) big vector
c
c  NOTE
c
c  mx must be at leat equal to 14*n*m + 10*max(n,m)
c
c  * DIAGNOSTICS *
c
c  ierr = 0   no error
c         1    too patches            action: increase mpatch
c         2    too many coefficents   action: increase msul
c
c  * NOTES *
c
c if ierr /=0 an error message is sent to the unit associated to ito
c
c=END USAGE
C
C=BLOCK SOURCE
C
      subroutine surfb(ifa,rest,nssur,n,m,mpatch,msul,ito,
     1                 aps,keypa,ierr)
      implicit none
C
C
C  memory mapping for GETTA3
C
C   coor   tanret  choret  q    cs    a    b    c    t
C     3nm    (9nm)     2nm   3n    n    n    n    n    3n  <-- dimensions
C   1      IP1     IP2     IP3  IP4   IP5  IP6  IP7  IP8  <-- REST
C
C total: 14*n*m + 10*max(n,m)
C
      integer*8 i,id,ifa,ierr,ii,ip1,ip2,ip3,ip4,ip5,ip6,ip7,ip8
      integer*8 ipa,ito,j,jj,k,keypa(mpatch),m,mpatch,msul,n,npa
      integer*8 npatch,nssur,nn,ispt3(3,2)
      real*8 apatch(3,4,4),aps(3,msul),ch1,ch2,rest(*),xl
c
      data ispt3/0,0,0,0,0,0/
C
      ierr =0
      read(ifa,*)n,m
      nn=max(n,m)
      npatch = (n-1)*(m-1)
      if(npatch.gt.mpatch)then
        write(ito,*)' ERROR : SURFB_001'
        write(ito,*)' N. OF PATCHES IN SUR.',nssur,' =',npatch
        write(ito,*)' MAX ALLOWED =        ',mpatch
        write(ito,*)' ABORTED              '
        write(ito,*)' INCREASE MPATCH      '
        ierr =1
        return
      endif
      if(17*npatch.gt.msul)then
        write(ito,*)' ERROR : SURFB_002'
        write(ito,*)' TOTAL N. OF COEF.IN SURF.',nssur,'=',17*npatch
        write(ito,*)' MAX ALLOWED =        ',msul
        write(ito,*)' ABORTED              '
        write(ito,*)' INCREASE MSUL       '
        ierr=2
        return
      endif
      ip1 = 1+3*n*m
      ip2 = ip1+9*m*n
      ip3 = ip2 + 2*m*n
      ip4 = ip3 + 3*nn
      ip5 = ip4 + nn
      ip6 = ip5 + nn
      ip7 = ip6 + nn
      ip8 = ip7 + nn
      read(ifa,*)(rest(k),k=1,3*n*m)
      call getta3(3,n,m,nn,rest,rest(ip1),rest(ip2),rest(ip3),
     1            rest(ip4),rest(ip5),rest(ip6),rest(ip7),
     1            rest(ip8),ispt3)
      ipa = 1
      do 10 i=1,n-1
      do 10 j=1,m-1
        npa = j+(i-1)*(m-1)
        keypa(npa)=ipa
        ch1 = rest(12*n*m + 1 + (i-1)*2 + (j-1)*2*n)
        ch2 = rest(12*n*m + 2 + (i-1)*2 + (j-1)*2*n)
c        xl=0.5*(choret(1,i,j)+choret(2,i,j))
        xl=0.5d0*(ch1+ch2)
        call evapa2(n,m,3,rest(ip1),rest(ip2),rest,apatch,i,j)
        aps(1,ipa)=4.d0
        aps(2,ipa)=4.d0
        aps(3,ipa)=xl
        ipa = ipa+1
c
        do 20 jj=1,4
        do 20 ii=1,4
        do 20 id=1,3
           aps(id,ipa)=apatch(id,ii,jj)
           ipa=ipa+1
20      continue
c
10    continue
      return
      end
C
C=END SOURCE
