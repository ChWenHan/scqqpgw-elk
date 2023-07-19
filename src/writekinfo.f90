subroutine writekinfo
use modmain
use modgw
implicit none
! arguments
integer ik,i1,i2,i3,jk
integer(16), allocatable :: ivklist(:)
! allocatable arrays
call init0
call init1
allocate(ivklist(nkptnr))
open(50,file='REDUCEDK.OUT',form='FORMATTED')
do ik=1,nkpt
  write(50,'(I8,3G18.10)') ik,vkl(:,ik)
end do
close(50)
!
open(51,file='KMAP.OUT',form='FORMATTED')
do jk=1,nkptnr
  i1=ivk(1,jk);i2=ivk(2,jk);i3=ivk(3,jk)
  ivklist(jk)=ivkik(i1,i2,i3)
  write(51,'(I8,3G18.10,I8)') jk,vkl(:,jk),ivklist(jk)
end do
close(51)
deallocate(ivklist)
return
end subroutine

