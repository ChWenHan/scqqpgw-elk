
! Copyright (C) 2017 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine acgwse_whch(ist,jst,nr,sem,wr,ser)
use modmain
use modgw
implicit none
! arguments
integer,intent (in) :: ist,jst
integer , intent (in) :: nr
complex(8), intent(in) :: sem(nstsv,nstsv,0:nwfm)
real(8), intent(in) :: wr(nr)
complex(8), intent(out) :: ser(nstsv,nstsv,nr)
! allocatable arrays
complex(8), allocatable :: zm(:),zwr(:),zr(:)
allocate(zm(0:nwfm),zwr(nr),zr(nr))
!whch ---- ATTENTION. It's different with the original version.
!nr is the number of real w points. usually = nwplotgw + nstsv
!nr is a new input var.
zm(:)=sem(ist,jst,:)
zwr(:)=wr(:)
!@write(*,*)'in acgwse_whch'
select case(actype)
case(1)
! fit a multipole model
  call acpole_whch(nr,zm,zwr,zr)
case(10)
! stabilised Pade approximant
  call pades(nspade,swidth,nwfm+1,wfm,zm,nr,zwr,zr)
case default
  write(*,*)
  write(*,'("Error(acgwse): actype not defined : ",I8)') actype
  write(*,*)
  stop
end select
ser(ist,jst,:)=zr(:)
deallocate(zm,zwr,zr)
return
end subroutine

