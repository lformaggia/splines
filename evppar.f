C=NAME evppar
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/evppar.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:36 $
C=TYPE subroutine
C=AUTHOR $Author: forma $
C=PURPOSE Local
C=KEYWORDS Parametric_bicubic_spline Parametric_polynomial_surfaces Data_structure_interface
C
C=BLOCK ABSTRACT
C
c It converts the rapresentation of a patch ofa  bicubic spline surface
c into the general rapresentation for a polynomial patch.
c
c E' una vecchia versione ormai obsoleta.. tenuta solo per compatibilita'
c Usare evapa2 o evapa3 a cui si rimanda per la descrizione delle variabili
C
C=END ABSTRACT
C=BLOCK USAGE
C
C    call evppar(n,m,tanret,choret,coors,apatch,i,j,xl)
C
C=END USAGE
C
C=BLOCK SOURCE
C
      subroutine evppar(n,m,tanret,choret,coors,apatch,i,j,xl)
      dimension tanret(3,3,n,m),choret(2,n,m),coors(3,*)
      dimension apatch(3,4,4)
      xl=0.5*(choret(1,i,j)+choret(2,i,j))
      call evapa2(n,m,3,tanret,choret,coors,apatch,i,j)
      return
      end
