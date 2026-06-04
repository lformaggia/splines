C=NAME buisgf
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/buisgf.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:35 $
C=TYPE subroutine
C=PURPOSE General Local_to readca
C=KEYWORDS Parametric_polinomial_surfaces Data_base_interface
C=BLOCK ABSTRACT
c
c The routine writes a surface definition into the data base. It starts writing
c the coefficent information from position inksur in the file of coefficent
c
c
c=END ABSTRACT
c
c=BLOCK USAGE
c
c call buisgf(n,m,nlsur,ib,ngsur,keypa,aps,kfio,ifio,inksur)
c
c  INPUT      DIMENSION     DESCRIPTION
c
c n                         n. of knots along u
c m                         n. of knots along v
c ib                        surface indicator (normally used to indicte the bounday condition
c                           associated to the surface)
c nlsur                     local number to be attributed to surface
c ngsur                     global numbering to be attributed to the surface
c aps         3,*           array of coefficents
c keypa       (n-1)(m-1)    keys to array of  coefficents
c kfio                      i/o unit number associated to the file of keys
c ifio                      i/o unit number associated to the file of coefficents
c inksur                    position in the file of coefficient from which the
c                           information about the surface will be stored
c
c  OUTPUT     DIMENSION     DESCRIPTION
c
c  inksur                   new available position in the file of coefficents(see NOTES)
c
c NOTES
c--------
c
c The file associated to the file of keys and the file of coefficents must have already
c been opened. NO CHECKING is done in buisgf about it.
c
c The value of inksur returned by the routine correspond to the next available position
c in the file of coefficients, after storing the information relative to the surface. The
c routine DOES NOT write the value (-inksur) in the first field of the next available
c record in the file of keys. This  is required to signal the end of information on
c that file and the next available position in the file of coefficients. It is up to
c the user to do it, if necessary. If the next position in the files of keys is nlsur +1,
c after calling buisgf the user should write
c
c      WRITE(KFIO,REC=nlsur+1)-inksur
c
c=END USAGE
C
C=BLOCK SOURCE
C
      subroutine buisgf(n,m,nlsur,ib,ngsur,keypa,aps,
     1                  kfio,ifio,inksur)
c
c It writes surface description into data base
c
c     n,m = n. of knots in the two parametric directions
c     nlsur = local number for surface
c     nhsur = global surface ID number
c     ib    = boundary marker
c     aps(3,*) = surface coefficents array
c     keypa(*) = pointers to aps
c     kfio    = "key" file i/o channel number
c     ifio    = "data" file i/o channel number
c     inksur  = first available location in data file
c
      dimension aps(3,*)
      dimension keypa(*)
c
      ipa =0
c      write(kfio,rec=nlsur)inksur,n,m
      write(kfio,rec=nlsur)inksur,n,m,ngsur,ib
c
c start loop over patches
c
      do 10 i=1,(n-1)*(m-1)
        ipa = i
        kpa = keypa(ipa)
        xl  = aps(3,kpa)
        nu  = int(aps(1,kpa))
        nv  = int(aps(2,kpa))
        kpa = kpa + 1
        write(ifio,rec=inksur)nu,nv,xl
        inksur=inksur+1
c
c loop over corfficents
c
        do 20 jj=1,nv
        do 20 ii=1,nu
           write(ifio,rec=inksur)(aps(id,kpa),id=1,3)
           inksur = inksur+ 1
           kpa    = kpa   + 1
20      continue
10    continue
      return
      end
C
C=END SOURCE
