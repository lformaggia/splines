C=NAME ilocli
C=RCSfile $Source: /u/flosys3d/software/cvsroot/splines/ilocli.f,v $
C=VERSION $Revision: 1.1.1.1 $ $Date: 1996/03/05 17:16:38 $
C=TYPE integer function
C=AUTHOR L.formaggia
C=PURPOSE General
C=KEYWORDS Parametric_polynomial_curves Parametric_polynomial_surfaces
C
C=BLOCK ABSTRACT
C
c It locates the local number of a surface or curve, given the local-to-global pointer
c array
C
C=END ABSTRACT
C=BLOCK USAGE
c
c  ilocal = ilocli(item,list,nlist)
c
C  INPUT      DIMENSION     DESCRIPTION
c
c item                      item to be fount in list (global number)
c list                      list to be sarched
c nlist                     n. of elements in list
c
c
C  OUTPUT     DIMENSION     DESCRIPTION
c
c ilocli                    location in list corrisponding to item, i.e.
c                           item = list(ilocli)
c                           if ilocli=0 -> item not found
c
c=END USAGE
C
C=BLOCK SOURCE
C
      integer*8 function ilocli(item ,list,nlist)
      implicit none
      integer*8 i,item,list(*),nlist
      ilocli=0
      if(item.eq.0)return
      do 10 i=1,nlist
        if(list(i).eq.item)then
          ilocli =i
          return
        endif
10    continue
      return
      end
C
C=END SOURCE
