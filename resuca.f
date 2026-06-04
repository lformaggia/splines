C=NAME resuca
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/resuca.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:40 $
C=TYPE subroutine
C=PURPOSE general Local_to readca
C=KEYWORDS Parametric_Polynomial_Surfaces Data_Base_Interface CATIA_neutral_file
C=BLOCK ABSTRACT
C
C   The routine reads a SINGLE surface from a catia neutral file and stores the relevant
c   information into the Database
C
C
C=END ABSTRACT
C=BLOCK USAGE
c
c         call resuca(ioc,aps,keypa,n,m,mpatch,msul,isoc,ito,ierr,linen)
c
c  INPUT      DIMENSION     DESCRIPTION
c
c msul                      max. dimension of  file of coefficents for surface (aps)
c marc                      max. dimension of  file of keys for curve (isn)= max. n. of arcs
c mpatch                    max. dimension of  file of keys for surface (keypa) = max. n. of patches
c ioc                       i/o unit number associated to catia file
c ito                       i/o unit number associated to log/error file
c linen                     number of last line read from CATIA neutral file (it is incremented by 1 for each
c                           line read)
c
c  OUTPUT     DIMENSION     DESCRIPTION
c
c m,n                       number of knots along u and v for the surface
c aps         3,msul        array of coefficents for surface
c keypa       mpatch        keys to array of  coefficents for surface
c icati                     Catia surface ID number (normally taken as global surface number)
c linen                     number of last line read from CATIA neutral file (it is incremented by 1 for each
c                           line read)
c ierr                      error indicator (see NOTES and DIAGNOSTICS)
c
c
c *NOTES and DIAGNOSTICS*
c
c
c ierr =0     no error
c       1     too many patches.
c              memory overflow in keypa: ABORT.
c              action: increase mpatch
c       2     too many coefficents.
c              memory overflow in aps: ABORT.
c              action: increase msul
c       3     1 + 2
c  The unit associated to ito is used for the logging of mesages. In standard FORTRAN
c  if you want the messages on the standard output you should put ito=6
c
c=END USAGE
C
C=BLOCK SOURCE
C
      subroutine resuca(ioc,aps,keypa,n,m,mpatch,msul,
     1                  isoc,ito,ierr,linen)
c
c Reads Catia surface definition
c
c
c   ioc = CATIA file input channel n.
c   aps(3,msul) = surface coefficent array
c   keypa(mpatch) = pointers to aps
c   n,m      = number of knots in the two parametric directions
c   mpatch   = max n. of patches allowed
c   msul     = max. global n. of coefficents
c   isoc     = CATIA surface ID number
c   ito      = standard output channel
c   ierr     = error condition (0= no error)
c   linen = line number counter
c
      character*6 field,line*80,cha*1
      real aps(3,msul),x1(3),x2(3)
      integer keypa(*)
C
      ipa =1
      ierr =0
      mcoef = msul
      read(ioc,'(A)')line
      linen = linen +1
      write(ito,'(A)')line
      field=line(1:6)
      if(field.ne.' NI NJ')then
        write(ito,*)' THE KEYWORD DOESNT CORRESPOND TO A SURFACE'
        write(ito,*)' RESUCA -001 : WRONG KEYWORD'
        ierr = 1
        return
      endif
      read(line,'(T9,2I6)')ni,nj
      n = ni+1
      m = nj+1
      npatch = ni*nj
      isoc =0
      do 456 i=23,77
        if(line(i:i+3).eq.'*SUR')then
          do 457 j=i+4,80
            if(line(j:j).ne.' ')then
              cha=line(j:j)
              read(cha,'(i1)',err=998)ix
              isoc = isoc*10 + ix
            else
              go to 999
            endif
457       continue
        endif
456     continue
458     continue
998     write(ito,*)' WARNING IN RESUCA '
        write(ito,*)' CANNOT READ CATIA SURFACE ID.'
999     continue
        if(npatch.gt.mpatch)then
        write(ito,*)' ',npatch,' PATCHES ON SURFACE'
        write(ito,*)' MAXIMUM ALLOWED IS: ',mpatch
        ierr = 2
        return
      endif
      read(ioc,'(A)')line
      linen = linen +1
c
c read all patches
c
      do 10 i= 1,npatch
         keypa(i)=ipa
         read(ioc,'(T19,2I6)')nu,nv
         linen = linen+1
         read(ioc,'(A)')line
         linen = linen +1
         read(ioc,'(A)')line
         linen = linen +1
         kpa = ipa
         aps(1,kpa)=float(nu)
         aps(2,kpa)=float(nv)
         ipa = ipa +1
         ixx =0
         do 20 k=1,nu*nv
            if(ipa.gt.mcoef)then
             ixx=1
             go to 21
            endif
            read(ioc,'(9X,3(f16.10,1x))')x,y,z
c           read(ioc,'(9X,3F16.5)')x,y,z
            linen = linen +1
            aps(1,ipa) = x
            aps(2,ipa) = y
            aps(3,ipa) = z
21          ipa = ipa +1
20       continue
c
         if(ixx.eq.1)goto 10
C
C EVALUATE THE PATCH TYPICAL LENGHT XL AS THE DIAGONAL LENGHT
C
         call evsurg(aps(1,kpa),0.,0.,nu,nv,x1,xx,xx,xx,xx,xx,1)
         call evsurg(aps(1,kpa),1.,1.,nu,nv,x2,xx,xx,xx,xx,xx,1)
         xl = sqrt((x1(1)-x2(1))**2+(x1(2)-x2(2))**2
     1        +(x1(3)-x2(3))**2)
         aps(3,kpa)=xl
10    continue
c
c error condition : too many coefficents on surface
c
      if(ixx.eq.1)then
         write(ito,*)' SURFACE CONTAINS TOO MANY COEFFICENTS'
         write(ito,*)' NEEDED AT LEAST MSUL = ',ipa
         write(ito,*)' SURFACE NOT READ!      '
         ierr = 4
      endif
      return
      end
C
C=END SOURCE
