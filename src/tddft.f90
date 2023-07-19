
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tddft
use modmain
use modtddft
use moddftu
use modmpi
use modomp
use modstore
use modtest
implicit none
if (tshift) then
  write(*,*)
  write(*,'("Error(tddft): use tshift = .false. for the ground-state run")')
  write(*,*)
  stop
end if
if ((task.eq.460).or.(task.eq.462)) then
! calculation starts at t=0
  tdt0=.true.
else
! calculation restarts
  tdt0=.false.
end if
! determine the static density and charge and write to file
if (tdt0) call rhostatic
! average force can be non-zero (allow for translation of atomic basis)
tfav0_=tfav0
tfav0=.false.
tforce_=tforce
tjr_=tjr
! initialise TDDFT variables
call tdinit
! set the stop signal to .false.
tstop=.false.
!---------------------------------!
!    main loop over time steps    !
!---------------------------------!
if (mp_mpi) write(*,*)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
do itimes=itimes0,ntimes-1
  if (mp_mpi) then
    write(*,'("Info(tddft): time step ",I8," of ",I8,",   t = ",G18.10)') &
     itimes,ntimes-1,times(itimes)
  end if
! reset the OpenMP thread variables
  call omp_reset
! evolve the wavefunctions across a single time step
  call timestep
! generate the density and magnetisation at current time step
  call rhomag
! compute the gauge-invariant current j(r) if required
  if (tjr) call genjr
! time step the induced A-field
  if (tafindt) call afindtstep
! calculate the electric field
  call genefieldt
! compute the time-dependent Kohn-Sham potentials and magnetic fields
  call potkst
! add the fixed spin moment effective field if required
  call addbfsm
! DFT+U
  if (dftu.ne.0) then
    call gendmatmt
    call genvmatmt
    call vmatmtsc
  end if
! compute the total energy
  call energytd
! write muffin-tin L, S and J if required
  if (tdlsj) call writetdlsj
! calculate the atomic forces if required
  if (tforce) then
    if ((itimes.eq.itimes0).or.(mod(itimes-1,ntsforce).eq.0)) then
      call force
    end if
  end if
! time step the atomic positions for Ehrenfest dynamics using forces calculated
! during the previous TDDFT run
  if (tdatpos) call atptstep
! write TDDFT output
  if (mp_mpi) call writetddft
! check for STOP file
  call checkstop
  if (tstop) exit
end do
filext='.OUT'
tfav0=tfav0_
tforce=tforce_
tjr=tjr_
tdatpos=.false.
! write the total current of the last step to test file
call writetest(460,'total current of last time step',nv=3,tol=5.d-4,rva=jtot)
end subroutine
