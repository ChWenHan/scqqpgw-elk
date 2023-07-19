
! Copyright (C) 2017 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine putgwclv(ik,v)
use modmain
use modgw
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: v(nstsv,nstsv)
! local variables
integer recl
! find the record length
inquire(iolength=recl) vkl(:,ik),nstsv,v
!$OMP CRITICAL(u280)
open(280,file='GWCLV.OUT',form='UNFORMATTED',access='DIRECT',recl=recl)
write(280,rec=ik) vkl(:,ik),nstsv,v
close(280)
!$OMP END CRITICAL(u280)
return
end subroutine
  
  