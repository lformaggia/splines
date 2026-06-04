C=NAME builtg
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/builtg.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:35 $
C=TYPE subroutine
C=AUTHOR L.Formaggia
C=PURPOSE General
C=KEYWORDS Parametric_polynimial_surfaces Parametric_polynomial_curves Data_Base_Interface
C
C=BLOCK ABSTRACT
C
c It scans the files of keys for the curves and surfaces in the Data base and it builds
c the array which give the link between local and global numbering
C
C=END ABSTRACT
C=BLOCK USAGE
C
C    call builtg(kfio,ksio,indg,inds,nsseg,nssur,mseg,msur,ierr,iendc)
c
C  INPUT      DIMENSION     DESCRIPTION
C
c kfio                      i/o unit number associated to the file of keys for the
c                           surfaces (0 if surface file is not considered)
c ksio                      i/o unit number associated to the file of keys for the
c                           curves  (0 if curvee file is not considered)
c mseg                      max number of curves allowed
c msur                      max number of surfaces allowed
c
C OUTPUT      DIMENSION     DESCRIPTION
c
c indg        mseg          locat-to-global link for curves: indg(i) is the global
C                           numbering of curve i
c indg        msur          locat-to-global link for surfaces: inds(i) is the global
C                           numbering of surface i
c nssur                     number of surfaces found
c nsseg                     number of curves found
c ierr                       error indicator (see DIAGNOSTICS)
c iendc                     end condition indicator (see DIAGNOSTICS)
C
C   * DIAGNOSTICS *
c
c ierr = 0    no error
c ierr = 1    error while reading curve file
c ierr = 2    max number of arcs exceeded
c ierr = 4     error reading surface file
c ierr = 8    max number of surfaces exceeded
c
c iendc = 0    OK
c iendc = 1    file of keys for curve does not end with a negative integer number
c             indicating the next available position in xln1d
c iendc = 2    file of keys for surfaces does not end with a negative integer number
c             indicating the next available position in aps
c
c ierr and iend error numbers sum up to indicate multiple error conditions
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine builtg(kfio,ksio,indg,inds,nsseg,nssur,mseg,msur,ierr,
     1                  iendc)
      integer indg(mseg),inds(msur)
      integer ios
c
c it sweeps over the file of keys in the DB and reports the
c number of records and the local to global array
c kfio = surface i/o channel
c ksio = curve segment i/o channel
c inds = locat to globag for surfaces
c indg = global to local for segments
c
      nssur =0
      nsseg =0
      ierr  =0
      iendc =0
c
      if(kfio.ne.0)then
10      continue
        read(kfio,rec=nsseg+1,iostat=ios)key
        if(ios.ne.0)then
           iendc = iendc + 1
           goto 130
        endif
        if(key.lt.0)then
           go to 130
        else
           nsseg = nsseg +1
           if(nsseg.gt.mseg)then
               ierr = ierr + 2
               goto 130
           endif
           read(kfio,rec=nsseg,iostat=ios)key,n,ng
           if(ios.ne.0)then
              ierr = ierr + 1
              goto 130
           endif
           indg(nsseg)=ng
           go to 10
        endif
130     continue
      endif
c
c
      if(ksio.ne.0)then
20      continue
        read(ksio,rec=nssur+1,iostat=ios)key
        if(ios.ne.0)then
           iendc = iendc + 2
           goto 230
        endif
        if(key.lt.0)then
           go to 230
        else
           nssur = nssur +1
           if(nssur.gt.msur)then
                 ierr = ierr + 8
                 goto 230
           endif
           read(ksio,rec=nssur,iostat=ios)key,n,m,ng
           if(ios.ne.0)then
              ierr = ierr + 4
              goto 230
           endif
           inds(nssur)=ng
           go to 20
        endif
230     continue
      endif
c
      return
      end
C
C=END SOURCE
c
