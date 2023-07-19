subroutine eveqngwdftk2(ik,evalsvp,evecsv,vse)
use modmain
use modgw
use moddftu
use modomp
implicit none
! arguments
integer,intent(in) :: ik
complex(8), intent(in) :: vse(nstsv,nstsv)
real(8), intent(out) :: evalsvp(nstsv)
complex(8), intent(out) :: evecsv(nstsv,nstsv)
! local variables
integer ld,ist,jst,ispn,i,j
integer jspn
real(8) eigtemp(nstsv),eig0(nstsv,nstsv)
real(8) ts0,ts1
complex(8) z1
character(256) fpre
character(256) fname
! allocatable arrays
! external functions
complex(4), external :: cdotc
complex(8), external :: zdotc
! no calculation of second-variational eigenvectors

call getevalsv(filext,ik,vkl(:,ik),eigtemp(1:nstsv))
!$OMP CRITICAL(u99_)
write(fname,'("EIG0PREV_K",I6.6,"_ISCL",I4.4".OUT")') ik,iscl
open(99,file=trim(fname),form='FORMATTED')
do i=1,nstsv
  write(99,'(I6.4,G18.10)')i,eigtemp(i)
end do
close(99)
!$OMP END CRITICAL(u99_)
eig0=zzero
do ist=1,nstsv
  eig0(ist,ist)=eigtemp(ist)
end do
!

fpre=trim('VSE_in_eveqn_for_check')
call writegwclv(ik,vse,fpre)
! previous eigven value plus ΔSe, H = ε_ii + SERx+SERc-Vxc
! stored as evecsv.
evecsv=zzero
evecsv=eig0+vse
fpre=trim('EIG0_as_totpot')
call writegwclv(ik,evecsv,fpre)
! !
if (spcpl.or.(.not.spinpol)) then
!spins are coupled; or spin-unpolarised: full diagonalisation
  call eveqnzh(nstsv,nstsv,evecsv,evalsvp)
else
! spins not coupled: block diagonalise H
  call eveqnzh(nstfv,nstsv,evecsv,evalsvp)
  i=nstfv+1
  call eveqnzh(nstfv,nstsv,evecsv(i,i),evalsvp(i))
  do j=1,nstfv
    do i=1,nstfv
      evecsv(i,j+nstfv)=0.d0
      evecsv(i+nstfv,j)=0.d0
    end do
  end do
end if
call timesec(ts1)
!$OMP ATOMIC
timesv=timesv+ts1-ts0
!$OMP CRITICAL(u99_)
write(fname,'("EIGAFTERDIAG_K",I6.6,"_ISCL",I4.4".OUT")') ik,iscl
open(99,file=trim(fname),form='FORMATTED')
do i=1,nstsv
  write(99,'(I6.4,G18.10)')i,evalsvp(i)
end do
close(99)
!$OMP END CRITICAL(u99_)
end subroutine
