
! Copyright (C) 2002-2016 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: hmlrad
! !INTERFACE:
subroutine hmlrad
! !USES:
use modmain
use modomp
! !DESCRIPTION:
!   Calculates the radial Hamiltonian integrals of the APW and local-orbital
!   basis functions. In other words, for atom $\alpha$, it computes integrals of
!   the form
!   $$ h^{\alpha}_{qq';ll'l''m''}=\begin{cases}
!    \int_0^{R_i}u^{\alpha}_{q;l}(r)H u^{\alpha}_{q';l'}(r)r^2dr & l''=0 \\
!    \int_0^{R_i}u^{\alpha}_{q;l}(r)V^{\alpha}_{l''m''}(r)
!    u^{\alpha}_{q';l'}(r)r^2dr & l''>0 \end{cases}, $$
!   where $u^{\alpha}_{q;l}$ is the $q$th APW radial function for angular
!   momentum $l$; $H$ is the Hamiltonian of the radial Schr\"{o}dinger equation;
!   and $V^{\alpha}_{l''m''}$ is the muffin-tin Kohn-Sham potential. Similar
!   integrals are calculated for APW-local-orbital and
!   local-orbital-local-orbital contributions.
!
! !REVISION HISTORY:
!   Created December 2003 (JKD)
!   Updated for compressed muffin-tin functions, March 2016 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ias,nthd
integer nr,nri,iro
integer ir,npi,i
integer l1,l2,l3,m2,lm2
integer io,jo,ilo,jlo
real(8) sm,t1
! begin loops over atoms and species
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(is,nr,nri,iro,npi) &
!$OMP PRIVATE(l1,l2,l3,io,jo,sm) &
!$OMP PRIVATE(lm2,m2,i,ir,t1,ilo,jlo) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  iro=nri+1
  npi=npmti(is)
!---------------------------!
!     APW-APW integrals     !
!---------------------------!
  do l1=0,lmaxapw
    do io=1,apword(l1,is)
      do l3=0,lmaxapw
        do jo=1,apword(l3,is)
          if (l1.eq.l3) then
            sm=sum(apwfr(1:nr,1,io,l1,ias)*apwfr(1:nr,2,jo,l3,ias) &
             *wrmt(1:nr,is))
            haa(1,jo,l3,io,l1,ias)=sm/y00
          else
            haa(1,jo,l3,io,l1,ias)=0.d0
          end if
          if (l1.ge.l3) then
            lm2=1
            do l2=1,lmaxi
              do m2=-l2,l2
                lm2=lm2+1
                sm=0.d0
                i=lm2
                do ir=1,nri
                  t1=apwfr(ir,1,io,l1,ias)*apwfr(ir,1,jo,l3,ias)*wrmt(ir,is)
                  sm=sm+t1*vsmt(i,ias)
                  i=i+lmmaxi
                end do
                do ir=iro,nr
                  t1=apwfr(ir,1,io,l1,ias)*apwfr(ir,1,jo,l3,ias)*wrmt(ir,is)
                  sm=sm+t1*vsmt(i,ias)
                  i=i+lmmaxo
                end do
                haa(lm2,jo,l3,io,l1,ias)=sm
                haa(lm2,io,l1,jo,l3,ias)=sm
              end do
            end do
            do l2=lmaxi+1,lmaxo
              do m2=-l2,l2
                lm2=lm2+1
                sm=0.d0
                i=npi+lm2
                do ir=iro,nr
                  t1=apwfr(ir,1,io,l1,ias)*apwfr(ir,1,jo,l3,ias)*wrmt(ir,is)
                  sm=sm+t1*vsmt(i,ias)
                  i=i+lmmaxo
                end do
                haa(lm2,jo,l3,io,l1,ias)=sm
                haa(lm2,io,l1,jo,l3,ias)=sm
              end do
            end do
          end if
        end do
      end do
    end do
  end do
!-------------------------------------!
!     local-orbital-APW integrals     !
!-------------------------------------!
  do ilo=1,nlorb(is)
    l1=lorbl(ilo,is)
    do l3=0,lmaxapw
      do io=1,apword(l3,is)
        if (l1.eq.l3) then
          sm=sum(lofr(1:nr,1,ilo,ias)*apwfr(1:nr,2,io,l3,ias)*wrmt(1:nr,is))
          hloa(1,io,l3,ilo,ias)=sm/y00
        else
          hloa(1,io,l3,ilo,ias)=0.d0
        end if
        lm2=1
        do l2=1,lmaxi
          do m2=-l2,l2
            lm2=lm2+1
            sm=0.d0
            i=lm2
            do ir=1,nri
              t1=lofr(ir,1,ilo,ias)*apwfr(ir,1,io,l3,ias)*wrmt(ir,is)
              sm=sm+t1*vsmt(i,ias)
              i=i+lmmaxi
            end do
            do ir=nri+1,nr
              t1=lofr(ir,1,ilo,ias)*apwfr(ir,1,io,l3,ias)*wrmt(ir,is)
              sm=sm+t1*vsmt(i,ias)
              i=i+lmmaxo
            end do
            hloa(lm2,io,l3,ilo,ias)=sm
          end do
        end do
        do l2=lmaxi+1,lmaxo
          do m2=-l2,l2
            lm2=lm2+1
            sm=0.d0
            i=npi+lm2
            do ir=iro,nr
              t1=lofr(ir,1,ilo,ias)*apwfr(ir,1,io,l3,ias)*wrmt(ir,is)
              sm=sm+t1*vsmt(i,ias)
              i=i+lmmaxo
            end do
            hloa(lm2,io,l3,ilo,ias)=sm
          end do
        end do
      end do
    end do
  end do
!-----------------------------------------------!
!     local-orbital-local-orbital integrals     !
!-----------------------------------------------!
  do ilo=1,nlorb(is)
    l1=lorbl(ilo,is)
    do jlo=1,nlorb(is)
      l3=lorbl(jlo,is)
      if (l1.eq.l3) then
        sm=sum(lofr(1:nr,1,ilo,ias)*lofr(1:nr,2,jlo,ias)*wrmt(1:nr,is))
        hlolo(1,jlo,ilo,ias)=sm/y00
      else
        hlolo(1,jlo,ilo,ias)=0.d0
      end if
      lm2=1
      do l2=1,lmaxi
        do m2=-l2,l2
          lm2=lm2+1
          sm=0.d0
          i=lm2
          do ir=1,nri
            t1=lofr(ir,1,ilo,ias)*lofr(ir,1,jlo,ias)*wrmt(ir,is)
            sm=sm+t1*vsmt(i,ias)
            i=i+lmmaxi
          end do
          do ir=iro,nr
            t1=lofr(ir,1,ilo,ias)*lofr(ir,1,jlo,ias)*wrmt(ir,is)
            sm=sm+t1*vsmt(i,ias)
            i=i+lmmaxo
          end do
          hlolo(lm2,jlo,ilo,ias)=sm
        end do
      end do
      do l2=lmaxi+1,lmaxo
        do m2=-l2,l2
          lm2=lm2+1
          sm=0.d0
          i=npi+lm2
          do ir=iro,nr
            t1=lofr(ir,1,ilo,ias)*lofr(ir,1,jlo,ias)*wrmt(ir,is)
            sm=sm+t1*vsmt(i,ias)
            i=i+lmmaxo
          end do
          hlolo(lm2,jlo,ilo,ias)=sm
        end do
      end do
    end do
  end do
! end loops over atoms and species
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end subroutine
!EOC

