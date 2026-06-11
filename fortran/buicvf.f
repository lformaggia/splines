C=NAME buicvf
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/buicvf.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:35 $
C=TYPE subroutine
C=PURPOSE General Local_to readca
C=KEYWORDS Parametric_polinomial_curves Data_base_interface
C=BLOCK ABSTRACT
c
c The routine writes a curve definition into the data base. It starts writing
c the coefficent information from position inkcur in the file of coefficent
c
c
c=END ABSTRACT
c
c=BLOCK USAGE
c
c  call buicvf(n,nlcur,ngcur,xln1d,isn,ksio,isio,inkcur,ug1,ug2,il1,il2,xg)
c
c  INPUT      DIMENSION     DESCRIPTION
c
c n                         n. of knots (n. of arcs+1)
c nlcur                     local number to be attributed to curve
c ngcur                     global numbering to be attributed to the curve
c xln1d       3,*           array of coefficents
c isn         n-1           keys to array of  coefficents
c ksio                      i/o unit number associated to the file of keys
c isio                      i/o unit number associated to the file of coefficents
c inkcur                    position in the file of coefficient from which the
c                           information about the curve will be stored
c ug1,ug2                   global coordinates of end points for definition segment
c il1,il2                   local  numbering of the surfaces adjacent to the curve
c xg                        control parameter (optional)
c
c  OUTPUT     DIMENSION     DESCRIPTION
c
c  inkcur                   new available position in the file of coefficents(see NOTES)
c
c
c NOTES
c--------
c
c The file associated to the file of keys and the file of coefficents must have already
c been opened. NO CHECKING is done in buicvf about it.
c
c The value of inkcur returned by the routine correspond to the next available position
c in the file of coefficients, after storing the information relative to the curve. The
c routine DOES NOT write the value (-inkcur) in the first field of the next available
c record in the file of keys. This is required to signal the end of information on
c that file and the next available position in the file of coefficients. It is up to
c the user to do it, if necessary. If the next position in the files of keys is nlcur +1,
c after calling buicvf the user should write
c
c      WRITE(KSIO,REC=nlcur+1)-inkcur
c
c
c=END USAGE
C
C=BLOCK SOURCE
C
      subroutine buicvf(narc,nlcur,ngcur,xln1d,isn,ksio,isio,inkcur,
     1                  ug1,ug2,il1,il2,xg)
      implicit none
c
c writes a curve definition into the data base
c
c  narc = n. of arcs +1
c  nlcur= local n. for curve
c  ngcur= global ID n.
c  xln1d(3,*)=curve coefficents

c  isn(narc) = pointers to xln1d
c  ksio      = "key" file i/o channel
c  isio      = "data" file i/o channel
c  inkcur    = pointer to first available location in "data" file
c  ug1,ug2   = parametric coordinates of curve limit points
c  il1,il2   = curve-surface links
c  xg        = curvature control parameter
c
      integer*8 i,id,ii,il1,il2,inkcur,ipa,isio,isn(*),kpa,ksio
      integer*8 narc,ngcur,nlcur,nu
      real*8 ug1,ug2,xg,xl,xln1d(3,*)
c
      ipa =0
c
      write(ksio,rec=nlcur)inkcur,narc,ngcur,ug1,ug2,il1,il2,xg
c
c start loop over arcs
c
      do 10 i=1,narc-1
        ipa = ipa+1
        kpa = isn(ipa)
        xl  = xln1d(2,kpa)
        nu  = int(xln1d(1,kpa))
        kpa = kpa + 1
        write(isio,rec=inkcur)nu,xl,0.d0
        inkcur=inkcur+1
c
c loop over coefficents
c
        do 20 ii=1,nu
           write(isio,rec=inkcur)(xln1d(id,kpa),id=1,3)
           inkcur = inkcur+ 1
           kpa    = kpa   + 1
20      continue
10    continue
      return
      end
C
C=END SOURCE
