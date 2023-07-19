
! Copyright (C) 2006 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnhfgwk(ikp,vmt,vir,bmt,bir,evecsvp)
use modmain
use modomp
implicit none
! arguments
integer, intent(in) :: ikp
real(8), intent(in) :: vmt(npcmtmax,natmtot),vir(ngtot)
real(8), intent(in) :: bmt(npcmtmax,natmtot,ndmag),bir(ngtot,ndmag)
complex(8), intent(inout) :: evecsvp(nstsv,nstsv)
! local variables
integer ist1
character(512) fpre
! automatic arrays
integer idx(nstsv)
! allocatable arrays
real(8), allocatable :: vgqc(:,:),gqc(:),gclgq(:),jlgqrmt(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: ylmgq(:,:),sfacgq(:,:)
complex(8), allocatable :: h(:,:),v(:,:),kmat(:,:)
complex(8), allocatable :: wfmt1(:,:,:,:),wfir1(:,:,:)
complex(8), allocatable :: wfmt2(:,:,:,:),wfir2(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:),zrhoir(:,:)
! external functions
complex(8), external :: zfinp
!$OMP CRITICAL(eveqnhf_)
write(*,'("Info(eveqnhfgwk): ",I6," of ",I6," k-points")') ikp,nkpt
!$OMP END CRITICAL(eveqnhf_)
! allocate local arrays
allocate(vgqc(3,ngvc),gqc(ngvc),gclgq(ngvc))
allocate(jlgqrmt(0:lnpsd,ngvc,nspecies))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
allocate(ylmgq(lmmaxo,ngvc),sfacgq(ngvc,natmtot))
allocate(h(nstsv,nstsv),v(nstsv,nstsv))
allocate(wfmt1(npcmtmax,natmtot,nspinor,nstsv),wfir1(ngtc,nspinor,nstsv))
allocate(wfmt2(npcmtmax,natmtot,nspinor,nstsv),wfir2(ngtc,nspinor,nstsv))
allocate(zrhomt(npcmtmax,natmtot,nstsv),zrhoir(ngtc,nstsv))
! get the first-variational eigenvectors from file for input reduced k-point
call getevecfv(filext,ikp,vkl(:,ikp),vgkl(:,:,1,ikp),evecfv)
! find the matching coefficients
call match(ngk(1,ikp),vgkc(:,:,1,ikp),gkc(:,1,ikp),sfacgk(:,:,1,ikp),apwalm)
! index to all states
do ist1=1,nstsv
  idx(ist1)=ist1
end do
! calculate the wavefunctions for all states of the input k-point
call genwfsv(.false.,.true.,nstsv,idx,ngdgc,igfc,ngk(1,ikp),igkig(:,1,ikp), &
  apwalm,evecfv,evecsvp,wfmt1,ngtc,wfir1)
!-----------------------------------------!
!     local potential matrix elements     !
!-----------------------------------------!
if (hybrid.and.spinpol) then
! magnetic field matrix elements in hybrid case
  call genvbmatk(vmt,vir,bmt,bir,ngk(1,ikp),igkig(:,1,ikp),wfmt1,ngtc,wfir1,h)
else
  !col matrix
  call genvmatk(vmt,vir,ngk(1,ikp),igkig(:,1,ikp),wfmt1,ngtc,wfir1,h)
end if
! Fourier transform wavefunctions to real-space
call zftwfir(ngk(1,ikp),igkig(:,1,ikp),wfir1)
!---------------------------------!
!     kinetic matrix elements     !
!---------------------------------!
allocate(kmat(nstsv,nstsv))
call getkmat(ikp,kmat)
call zgemm('N','N',nstsv,nstsv,nstsv,zone,kmat,nstsv,evecsvp,nstsv,zzero,v, &
  nstsv)
call zgemm('C','N',nstsv,nstsv,nstsv,zone,evecsvp,nstsv,v,nstsv,zone,h,nstsv)
deallocate(kmat)
fpre='HF_GWPOTH0'
call writegwclv(ikp,h,fpre)
!------------------------------!
!     Fock matrix elements     !
!------------------------------!
v(:,:)=0.d0
!generate Se_xc & Vh
!write(*,*)'gen se for ik: ',ikp
call gwsefmkhf(ikp,v)
fpre='HF_GWPOTVXC'
call writegwclv(ikp,v,fpre)
! add the Coulomb matrix elements to Hamiltonian
h(:,:)=h(:,:)+v(:,:)
fpre='HF_GWPOT'
call writegwclv(ikp,h,fpre)
!call writekscl2dfile(ikp,iscl,'POTH_XC',nstsv,nstsv,h)
!----------------------------------------------!
!     diagonalise Hartree-Fock Hamiltonian     !
!----------------------------------------------!
call eveqnzh(nstsv,nstsv,h,evalsv(:,ikp))
! apply unitary transformation to the third-variational states so that they
! refer to the first-variational basis
evecsv(:,:)=evecsvp(:,:)
call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,h,nstsv,zzero,evecsvp, &
  nstsv)
deallocate(evecsv,h,v)
end subroutine

