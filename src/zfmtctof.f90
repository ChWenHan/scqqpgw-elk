
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine zfmtctof(zfmt)
use modmain
use modomp
implicit none
! arguments
complex(8), intent(inout) :: zfmt(npmtmax,natmtot)
! local variables
integer is,ias,lm
integer nr,nri,nro
integer iro,ir,npi
integer nrc,nrci,nrco
integer irco,irc,npci
integer i,nthd
! automatic arrays
real(8) fi1(nrcmtmax),fi2(nrcmtmax)
real(8) fo1(nrmtmax),fo2(nrmtmax)
complex(8) zfmt1(npcmtmax)
if (lradstp.eq.1) return
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(zfmt1,fi1,fi2,fo1,fo2) &
!$OMP PRIVATE(is,nr,nri,nro,iro,npi) &
!$OMP PRIVATE(nrc,nrci,nrco,irco,npci) &
!$OMP PRIVATE(lm,i,irc,ir) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  nro=nr-nri
  iro=nri+1
  npi=npmti(is)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  nrco=nrc-nrci
  irco=nrci+1
  npci=npcmti(is)
! copy the input function
  zfmt1(1:npcmt(is))=zfmt(1:npcmt(is),ias)
! interpolate up to lmaxi over entire muffin-tin
  do lm=1,lmmaxi
    i=lm
    do irc=1,nrci
      fi1(irc)=dble(zfmt1(i)); fi2(irc)=aimag(zfmt1(i))
      i=i+lmmaxi
    end do
    do irc=irco,nrc
      fi1(irc)=dble(zfmt1(i)); fi2(irc)=aimag(zfmt1(i))
      i=i+lmmaxo
    end do
    call rfinterp(nrc,rcmt(:,is),wcrcmt(:,:,is),fi1,nr,rlmt(:,1,is),fo1)
    call rfinterp(nrc,rcmt(:,is),wcrcmt(:,:,is),fi2,nr,rlmt(:,1,is),fo2)
    i=lm
    do ir=1,nri
      zfmt(i,ias)=cmplx(fo1(ir),fo2(ir),8)
      i=i+lmmaxi
    end do
    do ir=iro,nr
      zfmt(i,ias)=cmplx(fo1(ir),fo2(ir),8)
      i=i+lmmaxo
    end do
  end do
! interpolate up to lmaxo on outer part of muffin-tin
  do lm=lmmaxi+1,lmmaxo
    i=npci+lm
    do irc=irco,nrc
      fi1(irc)=dble(zfmt1(i)); fi2(irc)=aimag(zfmt1(i))
      i=i+lmmaxo
    end do
    call rfinterp(nrco,rcmt(irco,is),wcrcmt(:,irco,is),fi1(irco),nro, &
     rsp(iro,is),fo1(iro))
    call rfinterp(nrco,rcmt(irco,is),wcrcmt(:,irco,is),fi2(irco),nro, &
     rsp(iro,is),fo2(iro))
    i=npi+lm
    do ir=iro,nr
      zfmt(i,ias)=cmplx(fo1(ir),fo2(ir),8)
      i=i+lmmaxo
    end do
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end subroutine

