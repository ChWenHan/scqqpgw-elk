
! Copyright (C) 2021 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine wfmtsv_sp(tsh,lrstp,is,ias,nst,idx,ngp,apwalm,evecfv,evecsv,ld,wfmt)
use modmain
use modomp
implicit none
! arguments
logical, intent(in) :: tsh
integer, intent(in) :: lrstp,is,ias,nst,idx(*),ngp(nspnfv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv)
integer, intent(in) :: ld
complex(4), intent(out) :: wfmt(ld,nspinor,nst)
! local variables
logical tasv
integer io,ilo,ispn,jspn
integer nr,nri,nro,iro
integer l,lm,np,npi
integer n,p,i,j,k,nthd
complex(8) zq(2),z1
! automatic arrays
complex(8) x(nstfv,nspnfv),y(nlmwf(is),nspinor,nst)
complex(4) cfmt(npmtmax)
! external functions
complex(8), external :: zdotu
iro=nrmti(is)+lrstp
if (lrstp.eq.1) then
  nr=nrmt(is)
  nri=nrmti(is)
  np=npmt(is)
  npi=npmti(is)
else
  nr=nrcmt(is)
  nri=nrcmti(is)
  np=npcmt(is)
  npi=npcmti(is)
end if
nro=nr-nri
! de-phasing factor for spin-spirals
if (ssdph) then
  zq(1)=zqss(ias)
  zq(2)=conjg(zq(1))
end if
! check if all the second-variational wavefunctions should be calculated
if (idx(1).eq.0) then
  tasv=.true.
else
  tasv=.false.
end if
!-----------------------!
!     APW functions     !
!-----------------------!
p=0
do l=0,lmaxo
  do lm=l**2+1,(l+1)**2
    do io=1,apword(l,is)
      p=p+1
      if (tevecsv) then
        do jspn=1,nspnfv
          n=ngp(jspn)
          do j=1,nstfv
            x(j,jspn)=zdotu(n,evecfv(:,j,jspn),1,apwalm(:,io,lm,ias,jspn),1)
          end do
        end do
! loop only over required states
        do j=1,nst
! index to state in evecsv
          if (tasv) then; k=j; else; k=idx(j); end if
          y(p,1,j)=zdotu(nstfv,evecsv(1,k),1,x,1)
          if (spinpol) then
            jspn=jspnfv(2)
            y(p,2,j)=zdotu(nstfv,evecsv(nstfv+1,k),1,x(1,jspn),1)
          end if
        end do
      else
        do j=1,nst
          if (tasv) then; k=j; else; k=idx(j); end if
          y(p,1,j)=zdotu(ngp(1),evecfv(:,k,1),1,apwalm(:,io,lm,ias,1),1)
        end do
      end if
    end do
  end do
end do
call holdthd(nst,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ispn,p,l,lm,i,io,z1) &
!$OMP NUM_THREADS(nthd)
do j=1,nst
  wfmt(1:np,:,j)=0.d0
  do ispn=1,nspinor
    p=0
    do l=0,lmaxo
      do lm=l**2+1,(l+1)**2
        i=npi+lm
        do io=1,apword(l,is)
          p=p+1
          z1=y(p,ispn,j)
          if (ssdph) z1=z1*zq(ispn)
          if (l.le.lmaxi) then
            call cfzrf(nri,z1,lrstp,apwfr(1,1,io,l,ias),lmmaxi,wfmt(lm,ispn,j))
          end if
          call cfzrf(nro,z1,lrstp,apwfr(iro,1,io,l,ias),lmmaxo,wfmt(i,ispn,j))
        end do
      end do
    end do
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
!---------------------------------!
!     local-orbital functions     !
!---------------------------------!
p=0
do ilo=1,nlorb(is)
  l=lorbl(ilo,is)
  do lm=l**2+1,(l+1)**2
    p=p+1
    i=idxlo(lm,ilo,ias)
    if (tevecsv) then
      do jspn=1,nspnfv
        n=ngp(jspn)
        do j=1,nstfv
          x(j,jspn)=evecfv(n+i,j,jspn)
        end do
      end do
      do j=1,nst
        if (tasv) then; k=j; else; k=idx(j); end if
        y(p,1,j)=zdotu(nstfv,evecsv(1,k),1,x,1)
        if (spinpol) then
          jspn=jspnfv(2)
          y(p,2,j)=zdotu(nstfv,evecsv(nstfv+1,k),1,x(1,jspn),1)
        end if
      end do
    else
      do j=1,nst
        if (tasv) then; k=j; else; k=idx(j); end if
        y(p,1,j)=evecfv(ngp(1)+i,k,1)
      end do
    end if
  end do
end do
call holdthd(nst,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(ispn,p,ilo,l,lm,i,z1) &
!$OMP NUM_THREADS(nthd)
do j=1,nst
  do ispn=1,nspinor
    p=0
    do ilo=1,nlorb(is)
      l=lorbl(ilo,is)
      do lm=l**2+1,(l+1)**2
        p=p+1
        i=npi+lm
        z1=y(p,ispn,j)
        if (ssdph) z1=z1*zq(ispn)
        if (l.le.lmaxi) then
          call cfzrf(nri,z1,lrstp,lofr(1,1,ilo,ias),lmmaxi,wfmt(lm,ispn,j))
        end if
        call cfzrf(nro,z1,lrstp,lofr(iro,1,ilo,ias),lmmaxo,wfmt(i,ispn,j))
      end do
    end do
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
if (tsh) return
! convert to spherical coordinates if required
call holdthd(nst,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(cfmt,ispn) &
!$OMP NUM_THREADS(nthd)
do j=1,nst
  do ispn=1,nspinor
    cfmt(1:np)=wfmt(1:np,ispn,j)
    call cbsht(nr,nri,cfmt,wfmt(:,ispn,j))
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
return

contains

pure subroutine cfzrf(n,z,ld1,rf,ld2,cf)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: z
integer, intent(in) :: ld1
real(8), intent(in) :: rf(ld1,n)
integer, intent(in) :: ld2
complex(4), intent(inout) :: cf(ld2,n)
cf(1,:)=cf(1,:)+z*rf(1,:)
end subroutine

end subroutine

