
! Copyright (C) 2006 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnhfgw(vmt,vir,bmt,bir)
use modmain
use modgw
use modomp
use modmpi
implicit none
!
real(8), intent(in) :: vmt(npcmtmax,natmtot),vir(ngtot)
real(8), intent(in) :: bmt(npcmtmax,natmtot,ndmag),bir(ngtot,ndmag)
! arguments
! local variables
integer ik,nthd,lp
! automatic arrays
! allocatable arrays
complex(8), allocatable :: evecsv(:,:)
! external functions
complex(8) zfinp
external zfinp
!$OMP CRITICAL(eveqnhf_)
!write(*,'("Info(eveqnhf): ",I6," of ",I6," k-points")') ikp,nkpt
!$OMP END CRITICAL(eveqnhf_)
call mpi_barrier(mpicom,ierror)
!--------------------------------------------
call holdthd(nkpt/np_mpi,nthd)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evecsv) &
!$OMP NUM_THREADS(nthd)
allocate(evecsv(nstsv,nstsv))
!$OMP DO
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  call getevecsv(filext,ik,vkl(:,ik),evecsv)
! whch---- solve the Hartree-Fock eigenvalue equation with Σ_c
  ! the electron hamiltonian now is: Ve-n + Vh + Vx + Σ_c
  call eveqnhfgwk(ik,vmt,vir,bmt,bir,evecsv)
! write the eigenvalues/vectors to file
  call putevalsv(filext,ik,evalsv(:,ik))
  call putevecsv(filext,ik,evecsv)
end do
!$OMP END DO
deallocate(evecsv)
!$OMP END PARALLEL
call freethd(nthd)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! broadcast eigenvalue array to every process
do ik=1,nkpt
  lp=mod(ik-1,np_mpi)
  call mpi_bcast(evalsv(:,ik),nstsv,mpi_double_precision,lp,mpicom,ierror)
end do
return
end subroutine

