
! Copyright (C) 2016 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine se2vsk2(thr,twr,ik,sec,vse)
use modmain
use modgw
use modomp
implicit none
! arguments
logical, intent(in) :: thr,twr !make it hermitian or take the real part.
integer, intent(in) :: ik
complex(8), intent(in) :: sec(nstsv,nstsv,0:nwfm)
complex(8), intent(out) :: vse(nstsv,nstsv)
! local variables
integer ist,jst,iw
integer nthd,nac
real(8) ts0,ts1
character(512) fname1
! allocatable arrays
complex(8), allocatable :: ser(:,:,:)
real(8), allocatable :: wr(:)
!
!try k-dependent se
if (twr) then
!$OMP CRITICAL(ser_)
  if (((nwrite.ge.1).and.(mod(iscl,nwrite).eq.0)).or.(iscl.eq.1)) then
    write(fname1,'("SEIW_K",I6.6,"_ISCL",I4.4,".OUT")') ik,iscl
    open(101,file=trim(fname1),form='FORMATTED',STATUS='REPLACE')
    do iw=0,nwfm
      do ist=1,nstsv
        do jst=1,nstsv
          write(101,'(I6.6,2I6.3,2G18.10)') iw,ist,jst,sec(ist,jst,iw)
        end do
      end do
    end do
    close(101)
  end if
!$OMP END CRITICAL(ser_)
end if
nac=nstsv
allocate(wr(nac))
do iw = 1 , nac
  !wr(iw)=evalsv(iw,ik)-efermi
  wr(iw)=evalsv(iw,ik)
end do
!$OMP CRITICAL(ser_)
if (((nwrite.ge.1).and.(mod(iscl,nwrite).eq.0)).or.(iscl.eq.1)) then
  write(fname1,'("WR_K",I6.6,"_ISCL",I4.4,".OUT")') ik,iscl
  open(101,file=trim(fname1),form='FORMATTED',STATUS='REPLACE')
  do iw=1,nac
        write(101,'(I6.6,G18.10)') iw,wr(iw)
  end do
  close(101)
end if
!$OMP END CRITICAL(ser_)


! sec now  imaginary w Σ.
allocate(ser(nstsv,nstsv,nac))
ser(:,:,:)=0.d0
!write(*,*)'in se2vsk2'
call timesec(ts0)
call holdthd(nstsv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ist) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do jst=1,nstsv
  do ist=1,nstsv
    if (tsediag.and.(ist.ne.jst)) cycle
! perform analytic continuation from the imaginary to the real axis
    call acgwse_whch(ist,jst,nac,sec,wr,ser)
  end do
end do
!$OMP END DO
!$OMP END PARALLEL
call freethd(nthd)
call timesec(ts1)
!$OMP ATOMIC
timeac=timeac+ts1-ts0
if (twr) then
!$OMP CRITICAL(ser_)
  if (((nwrite.ge.1).and.(mod(iscl,nwrite).eq.0)).or.(iscl.eq.1)) then
    write(fname1,'("SER_K",I6.6,"_ISCL",I4.4,".OUT")') ik,iscl
    open(101,file=trim(fname1),form='FORMATTED',STATUS='REPLACE')
    do iw=1,nac
      do ist=1,nstsv
        do jst=1,nstsv
          write(101,'(I6.6,2I6.3,2G18.10)') iw,ist,jst,ser(ist,jst,iw)
        end do
      end do
    end do
    close(101)
  end if
!$OMP END CRITICAL(ser_)
end if
vse(:,:)=0.d0
! w-independent se.
if (thr) then
  do ist=1,nstsv
    do jst=1,nstsv
      ! make SE w-independent. method from DOI: 10.1103/PhysRevLett.93.126406 | eq 2
      !Se(ist,jst)=(Se_ij(e_i) + Se_ij(e_j))/2
      vse(ist,jst)=0.5d0*(ser(ist,jst,ist)+ser(ist,jst,jst))
    end do 
  end do
else
  do ist=1,nstsv
    vse(ist,ist)=real(ser(ist,ist,ist))
  end do
end if

! ! make it Hermitian~
! do ist=1,nstsv
!   do jst=1,nstsv
!     vse(ist,jst)=0.5d0*(vse(ist,jst)+conjg(vse(jst,ist)))
!   end do
! end do
! mixing. it's now linear mixing
! allocate(vse0(nstsv,nstsv))
! if (iscl.gt.1) then
!   dv=0.d0
!   z10=0.d0
!   call getgwclv(ik,vse0)
!   do ist=1,nstsv
!     do jst=1,nstsv
!       z1=vse(ist,jst)-vse0(ist,jst)
!       if ((real(z1)**2+aimag(z1)**2).gt.(real(z10)**2+aimag(z10)**2)) z10=z1
!       dv=dv+sqrt((real(z1)**2+aimag(z1)**2))
!     end do
!   end do
! else
!   vse0=0.d0
! end if
! dvse(ik,iscl)=dv
! vse=alphagw*vse+(1-alphagw)*vse0
! ! try adapt linear
! !if (dv.lt.dvse(ik,iscl)) betavse=betavse+0.05.d0
! call putgwclv(ik,vse)
! call writegwclv(ik,vse)
deallocate(ser,wr)
return
end subroutine

