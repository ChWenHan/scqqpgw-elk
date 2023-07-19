
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rfmtctof
! !INTERFACE:
subroutine rfmtctof(rfmt)
! !USES:
use modmain
use modomp
! !INPUT/OUTPUT PARAMETERS:
!   rfmt : real muffin-tin function (in,real(npmtmax,natmtot))
! !DESCRIPTION:
!   Converts a real muffin-tin function from a coarse to a fine radial mesh by
!   using cubic spline interpolation. See {\tt rfinterp} and {\tt spline}.
!
! !REVISION HISTORY:
!   Created October 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(inout) :: rfmt(npmtmax,natmtot)
! local variables
integer is,ias,lm,nthd
integer nr,nri,nro
integer iro,ir,npi,i
integer nrc,nrci,nrco
integer irco,irc,npci
! automatic arrays
real(8) rfmt1(npcmtmax),fr(nrmtmax),frc(nrcmtmax)
if (lradstp.eq.1) return
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt1,fr,frc,is,nr,nri,nro) &
!$OMP PRIVATE(iro,npi,nrc,nrci,nrco) &
!$OMP PRIVATE(irco,npci,lm,i,irc,ir) &
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
  rfmt1(1:npcmt(is))=rfmt(1:npcmt(is),ias)
! interpolate up to lmaxi over entire muffin-tin
  do lm=1,lmmaxi
    i=lm
    do irc=1,nrci
      frc(irc)=rfmt1(i)
      i=i+lmmaxi
    end do
    do irc=irco,nrc
      frc(irc)=rfmt1(i)
      i=i+lmmaxo
    end do
    call rfinterp(nrc,rcmt(:,is),wcrcmt(:,:,is),frc,nr,rlmt(:,1,is),fr)
    i=lm
    do ir=1,nri
      rfmt(i,ias)=fr(ir)
      i=i+lmmaxi
    end do
    do ir=iro,nr
      rfmt(i,ias)=fr(ir)
      i=i+lmmaxo
    end do
  end do
! interpolate up to lmaxo on outer part of muffin-tin
  do lm=lmmaxi+1,lmmaxo
    i=npci+lm
    do irc=irco,nrc
      frc(irc)=rfmt1(i)
      i=i+lmmaxo
    end do
    call rfinterp(nrco,rcmt(irco,is),wcrcmt(:,irco,is),frc(irco),nro, &
     rsp(iro,is),fr(iro))
    i=npi+lm
    do ir=iro,nr
      rfmt(i,ias)=fr(ir)
      i=i+lmmaxo
    end do
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end subroutine
!EOC

