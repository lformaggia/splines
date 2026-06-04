C=NAME readca
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/readca.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:40 $
C=TYPE subroutine
C=PURPOSE general
C=KEYWORDS Parametric_Polynomial_Curves Parametric_Polynomial_Surfaces Data_Base_Interface CATIA_neutral_file
C=BLOCK ABSTRACT
C
C   The routine reads a catia neutral file and stores the relevant information into the Database
C
C
C=END ABSTRACT
C=BLOCK USAGE
C
c       subroutine readca(nssur,nsseg,nbcs,aps,keypa,mssur,msseg,inds,indg,
c                   msel,msul,marc,mpatch,xln1d,isn,filen,xgeod,ibd,inksur,inkseg,iocat,
c                   ierr,ito,ksio,isio,kfio,ifio)
c
c  INPUT      DIMENSION     DESCRIPTION
c
c nssur                     number of support surfaces (it is increased by 1 each time a surface
c                           descrition is read)
c nsseg                     number of support curves (it is increased by 1 each time a curve
c                           descrition is read)
c nbcs                      number of definition curves (it is increased by 1 each time a curve
c                           descrition is read)
c mssur                     max number of support surfaces allowed
c msseg                     max. number of support segment allowed
c msel                      max. dimension of  file of coefficents for curve (xln1d)
c msul                      max. dimension of  file of coefficents for surface (aps)
c marc                      max. dimension of  file of keys for curve (isn)= max. n. of arcs
c mpatch                    max. dimension of  file of keys for surface (keypa) = max. n. of patches
c filen       char*(*)      Catia neutral file name
c xgeod                     default value for control parameter for all curves read (this parameter
c                           is used by the mesh generation algorithm)
c ibd                       default value for surface indicator for all surfaces read
c                           (normally used to indicate the bounday condition associated to the surface)
c inksur                    initial position in the file of coefficient  from which the
c                           information about the surfaces will be stored
c inkseg                    initial position in the file of coefficient from which the
c                           information about the curves will be stored
c iocat                     i/o unit number associated to catia file
c ito                       i/o unit number associated to log/error file
c ksio                      i/o unit number associated to the file of keys
c isio                      i/o unit number associated to the file of coefficents
c kfio                      i/o unit number associated to the file of keys
c ifio                      i/o unit number associated to the file of coefficents
c
c
c  OUTPUT     DIMENSION     DESCRIPTION
c
c nssur                     number of support surfaces (it has been increased by 1 each time a surface
c                           descrition was read)
c nsseg                     number of support curves (it has been increased by 1 each time a curve
c                           descrition was read)
c nbcs                      number of definition curves (it has been increased by 1 each time a curve
c                           descrition was read)
c inds        mssur         local-to-global index for surfaces : inds(i) contains the global number associated to surface i
c indg        msseg         local-to-global index for curves : indg(i) contains the global number associated to curve i
c inkcur                    new available position in the file of coefficents(see NOTES)
c inksur                    new available position in the file of coefficents(see NOTES)
c
c HELP ARRAYS
c
c xln1d       3,msel        array of coefficents for curve
c isn         marc          keys to array of  coefficents for curve
c aps         3,msul        array of coefficents for surface
c keypa       mpatch        keys to array of  coefficents for surface
c
c NOTES
c
C (1) The routine checks if the i/o channel associated to iocat in the  following way:
c
c   if  file=filen does not exist
c   then
c          ABORT
c   endif
c
c   if file is already opened
c   then
c        if unit associated to file is iocat
c        then
c               DO NOTHING
c        else
c               close file=filen
c               if unit=iocat is opened
c               then
c                    close iocat
c               endif
c               open unit=iocat file=filen
c        endif
c   else
c        if unit=iocat is open
c        then
c             close iocat
c        endif
c        open unit=iocat file=filen
c   endif
c
c Thus, if the file was already opened the routiner starts reading from the next record.
c
c (2) The files associated to the data base ARE NOT OPENED inside the routine. It must be done outside it
c with the appropriate commands, i.e.
c
c       open(unit=ksio,file=keyseg,access='DIRECT',form='UNFORMATTED',
c     1      recl=8)
c       open(unit=isio,file=infseg,access='DIRECT',form='UNFORMATTED',
c     1      recl=3)
c       open(unit=kfio,file=keysur,access='DIRECT',form='UNFORMATTED',
c     1      recl=5)
c       open(unit=ifio,file=infsur,access='DIRECT',form='UNFORMATTED',
c     1      recl=3)
c
c (keyseg,infseg,... are the names of the files of keys and coefficents...)
c
c (3) The value of inkcur and inksur returned by the routine correspond to the next available position
c in the file of coefficients, after storing the information relative to the curves and surfaces read from catia file.
c After having read the CATIA file, the routine DOES NOT write anything in the first field of the next available
c record in the files of keys. This is normally required to signal the end of information on
c those files and to give the next available position in the files of coefficients. It is up to
c the user to do it, if necessary. Normally the local number of the last curcve and surface in the database after calling
c readca is nssur and nsseg respectively. To correctly store the data base then the following commands are required
c (After Calling readca!)
c
c      WRITE(KSIO,REC=nsseg+1)-inkcur
c      WRITE(KFIO,REC=nssur+1)-inksur
c
c (4) The unit associated to ito is used for the logging of mesages. In standard FORTRAN
c  if you want the messages on the standard output you should put ito=6
c
C=END USAGE
C
C=BLOCK SOURCE
C
       subroutine readca(nssur,nsseg,nbcs,aps,keypa,
     1                  mssur,msseg,inds,indg,
     1                  msel,msul,marc,mpatch,xln1d,isn,
     1                  filen,xgeod,ibd,inksur,inkseg,iocat,
     1                  ierr,ito,ksio,isio,kfio,ifio)
C
C ROUTINE DI LETTURA INPUT FILE E DI ORGANIZZAZIONE DATA BASE.
C
C VARIABILE                CONTENUTO                DIMENSIONE
C
C NSSUR       N. OF SUPPORT SURFACES
C NSSEG       N. OF SUPPORT SEGMENTS
C NBCS        N. OF DEFINITION SEGMENTS
C APS         ARRAY DEI COEFFICENTI PER SUPERFICE         (3,*)
C KEYPA       ARRAY DI POINTERS SU APS                    (*)
C ICHS        INFORMAZIONI SU I SEGMENTI DI DEFINIZIONE   (8,*)
C XB1D        PARAMETRIC COORDINATE END POINTS OF DEF.SEG.(2,*)
C IBOUN       BOUNDARY CONDITION MARKER
C MSSUR       MAX.N. SUPERFICI DI SUPPORTO
C MSSEG       MAX.N SEGMENTI DI SUPPORTO
C MBCSG       MAX. N. DEFINITION SEGMENTS
C INDS        LOCAL SURFACE NUMBER->GLOBAL SURFACE NUMBER
C INDG        LOCAL SEGMENT NUMBER->GLOBAL SEGMENT NUMBER
C ISASS       N. OF SURFACES ADJACENT TO A GIVEN CURVE     (2,*)
C MSEL        MAX N. OF   COEFFICENTS ON A SURFACE
C MSUL        MAX N. OF   COEFFICENTS ON A CURVE SEGMENT
C MARC        MAX N. OF   ARCS ON A CURVE
C MPATCH      MAX N. OF   PATCHES ON A SURFACE
C XB1D3       COORDINATES OF END POINTS FOR DEFINITION SEGMENT (3,2,NBCS)
C
C.............................................................................
C
       real aps(*),xln1d(3,*)
       integer inds(*),indg(*),keypa(*)
       integer isn(*)
C
       character*(*) filen
       logical isopen,exist,isof
       character*4 line*80
c
       ns0 = nssur
       ng0 = nsseg
       ierr =0
       linen=0
C
C OPEN FILES
C
       inquire(unit=iocat,opened=isopen)
       inquire(file=filen,exist=exist,number=iof,opened=isof)
       if(.not.exist)then
         write(ito,*)' READCA : ERROR 001 '
         write(ito,*)' CATIA Neutral File ' , filen
         write(ito,*)' DOES NOT EXIST '
         write(ito,*)' reding aborted **'
         ierr = 1
         return
       endif
       if(isof)then
         if(iocat.ne.iof)then
             close(iof)
             if(isopen)close(iocat)
             open(unit=iocat,file=filen,status='old',
     1        form='formatted')
         endif
       else
         if(isopen)close(iocat)
         open(unit=iocat,file=filen,status='old',form='formatted')
       endif
c
150    read(iocat,'(a)',end=1000,err=1000)line
       if(line(2:5).eq.'NARC')then
         backspace(iocat)
         call reseca(iocat,narc,xln1d,isn,msel,marc,il1,il2,ug1,
     1               ug2,icati,ito,ierr,linen)
         if(ierr.ne.0)then
            write(ito,*)' ERROR WHILE READING CURVE IN CATIA'
            write(ito,*)' NEUTRAL FILE ' ,filen
            write(ito,*)' READING ABORTED at line',linen
            write(ito,*)' error code ',ierr+10
            ierr = ierr + 10
         endif
         nbcs = nbcs +1
         nsseg =  nsseg + 1
         call buicvf(narc,nsseg,icati,xln1d,isn,ksio,isio,inkseg,
     1               ug1,ug2,il1,il2,xgeod)
         indg(nsseg) = icati
c
c
      else if (line(1:6) .eq. ' NI NJ')then
c
c
         backspace(iocat)
         call resuca(iocat,aps,keypa,n,m,mpatch,msul,
     1               isoc ,ito,ierr,linen)
         if(ierr.ne.0)then
            write(ito,*)' ERROR WHILE READING SURFACE IN CATIA'
            write(ito,*)' NEUTRAL FILE ' ,filen
            write(ito,*)' READING ABORTED at line',linen
            write(ito,*)' error code ',ierr + 20
            ierr = ierr + 20
         endif
         nssur = nssur +1
         inds(nssur)  = isoc
         call buisgf(n,m,nssur,ibn,isoc,keypa,aps,kfio,ifio,inksur)
       else
            write(ito,*)' ERROR WHILE READING  CATIA'
            write(ito,*)' NEUTRAL FILE ' ,filen
            write(ito,*)' READING ABORTED at line',linen
            ierr = 2
            write(ito,*)' error code ',2
       endif
c
                                                          GO TO 150
c
1000   continue
       write(ito,*)' CATIA neutral file ',filen,' read'
       write(ito,*)' n. of surfaces read =',nssur - ns0
       write(ito,*)' n. of curves   read =',nsseg - ng0
       return
       end
C
C=END SOURCE
