        subroutine TESTRE(ndimn,i1,i2,timp,tanret,n,m)
        real timp(ndimn,2),tanret(ndimn,3,n,m)
c
        do 530 id=1,ndimn
           timp(id      ,1) = tanret(id,3,i,i1)
           timp(id      ,2) = tanret(id,3,i,i2)
530     continue
        return
        end
