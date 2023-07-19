
! Copyright (C) 2017 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine writegwclv(ik,v,fpre)
use modmain
use modomp
use modgw
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: v(nstsv,nstsv)
character(128),intent(in) :: fpre
! local variables
integer ist1,ist2
character(512) fname1
logical twrite
call checkwrite(twrite)
!write(*,*)'',fpre
! find the record length
!$OMP CRITICAL(u98_)
if ((twrite).or.(nwrite.ge.1).and.(mod(iscl,nwrite).eq.0).or.(iscl.eq.1)) then
  write(fname1,'("_K",I6.6,"_ISCL",I4.4".OUT")') ik,iscl
  !write(*,*)'',trim(fpre)//trim(fname1)
  !write(*,*) fname1,ik
  open(98,file=trim(fpre)//trim(fname1),form='FORMATTED',STATUS='REPLACE')
  do ist1=1,nstsv
    do ist2=1,nstsv
      write(98,'(I6.3,I6.3,2G18.10)') ist1,ist2,v(ist1,ist2)
    end do
  end do
  close(98)
end if
!$OMP END CRITICAL(u98_)
return
end subroutine
  
  