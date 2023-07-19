subroutine eveqngwdft(vmt,vir,bmt,bir)
use modmain
use modgw
use modomp
use modmpi
implicit none
!
real(8), intent(in) :: vmt(npcmtmax,natmtot),vir(ngtot)
real(8), intent(in) :: bmt(npcmtmax,natmtot,ndmag),bir(ngtot,ndmag)

!
integer ik,nthd
!
real(8), allocatable :: evalfv(:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
complex(8), allocatable :: apwalm(:,:,:,:,:),vse(:,:)

!
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evalfv,evecfv,evecsv,vse) &
!$OMP NUM_THREADS(nthd)
allocate(evalfv(nstfv,nspnfv),vse(nstsv,nstsv))
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
!$OMP DO
do ik=1,nkpt
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  call gwsefmkdft(ik,vmt,vir,bmt,bir,vse)

  call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,1,ik),evecfv)

  call getevalfv(filext,ik,vkl(:,ik),evalfv)

  ! call eveqngwdftk(ik,ngk(1,ik),igkig(:,1,ik),vgkc(:,:,1,ik),evalfv,evecfv, &
  ! evalsv(:,ik),evecsv,vse)
  call eveqngwdftk2(ik,evalsv(:,ik),evecsv,vse)
  call putevalsv(filext,ik,evalsv(:,ik))
  call putevecsv(filext,ik,evecsv)
end do
!$OMP END DO
deallocate(evalfv,vse,evecfv,evecsv)
!$OMP END PARALLEL
call freethd(nthd)
end subroutine
