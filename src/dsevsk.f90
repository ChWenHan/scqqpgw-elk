
! Copyright (C) 2016 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine dsevsk(tdiag,ik,sec,dse)
use modmain
use modgw
use modomp
implicit none
! arguments
logical ,intent(in) :: tdiag !if just calc the diag of dse
integer, intent(in) :: ik
complex(8), intent(in) :: sec(nstsv,nstsv,0:nwfm)
complex(8), intent(out) :: dse(nstsv,nstsv)
! local variables
integer ist,jst,iw
integer nthd,nac
!real(8) ts0,ts1
character(512) fname1
real(8) eta
! allocatable arrays
complex(8), allocatable :: ser(:,:,:)
real(8), allocatable :: wr(:),dsei(:),dsej(:)
integer, allocatable :: idx1(:),idx2(:)
!
!try k-dependent se
nac=nstsv * 2
allocate(wr(nac),idx1(nstsv),idx2(nstsv),dsei(nstsv))
eta=swidth*1.d-4
do iw = 1,nstsv
  idx1(iw)=(iw-1)*2+1
  idx2(iw)=(iw-1)*2+2
  wr((iw-1)*2+1)=evalsv(iw,ik)-efermi-eta/2
  wr((iw-1)*2+2)=evalsv(iw,ik)-efermi+eta/2
end do
! sec now  imaginary w Σ.
allocate(ser(nstsv,nstsv,nac))
ser(:,:,:)=0.d0
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
!$OMP CRITICAL(ser_)
if (((nwrite.ge.1).and.(mod(iscl,nwrite).eq.0)).or.(iscl.eq.1)) then
  write(fname1,'("DSER_K",I6.6,"_SCL",I4.4,".OUT")') ik,iscl
  open(102,file=trim(fname1),form='FORMATTED',STATUS='REPLACE')
  do iw=1,nac
    do ist=1,nstsv
      do jst=1,nstsv
        write(102,'(I6.6,2I6.3,2G18.10)') iw,ist,jst,ser(ist,jst,iw)
      end do
    end do
  end do
  close(102)
end if
!$OMP END CRITICAL(ser_)
!write(*,*)'',swidth
dse(:,:)=1.d0
if (tdiag) then
  do ist=1,nstsv
    dsei(ist)=real(ser(ist,ist,idx2(ist))-ser(ist,ist,idx1(ist)))
    dse(ist,ist)=1/(1-(dsei(ist)/eta))
  end do
else
  allocate(dsej(nstsv))
  do ist=1,nstsv
    do jst=1,nstsv
      dsei(ist)=real(ser(ist,jst,(ist-1)*2+2)-ser(ist,jst,(ist-1)*2+1))/eta
      dsej(jst)=real(ser(ist,jst,(jst-1)*2+2)-ser(ist,jst,(jst-1)*2+1))/eta
      dse(ist,jst)=0.5d0*((1/(1-dsei(ist)))+(1/(1-dsej(jst))))
    end do
  end do
end if

!$OMP CRITICAL(dsei_)
if (((nwrite.ge.1).and.(mod(iscl,nwrite).eq.0)).or.(iscl.eq.1)) then
  write(fname1,'("DSEI_K",I6.6,"_SCL",I4.4,".OUT")') ik,iscl
  open(103,file=trim(fname1),form='FORMATTED',STATUS='REPLACE')
  do jst=1,nstsv
    write(103,'(I6.3,G18.10)') jst,dsei(jst)
  end do
  close(103)
end if
!$OMP END CRITICAL(dsei_)
deallocate(ser,wr,idx1,idx2)
return
end subroutine

