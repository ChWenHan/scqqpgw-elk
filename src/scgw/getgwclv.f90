
! Copyright (C) 2017 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine getgwclv(ik,v)
use modmain
use modgw
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(out) :: v(nstsv,nstsv)
! local variables
integer recl,nstsv_
real(8) vkl_(3),t1
! find the record length
inquire(iolength=recl) vkl_,nstsv_,v
!$OMP CRITICAL(u280)
open(280,file='GWCLV.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
read(280,rec=ik) vkl_,nstsv_,v
close(280)
!$OMP END CRITICAL(u280)
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getgwclv): differing vectors for k-point ",I8)') ik
  write(*,'(" current    : ",3G18.10)') vkl(:,ik)
  write(*,'(" GWSEFM.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(getgwsefm): differing nstsv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nstsv
  write(*,'(" GWSEFM.OUT : ",I8)') nstsv_
  write(*,*)
  stop
end if
return
end subroutine

