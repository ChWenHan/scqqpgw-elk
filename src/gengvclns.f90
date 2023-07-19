
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gengvclns
use modmain
use modtddft
implicit none
! local variables
integer is,ias
integer nr,nri,ir
integer np,i
! automatic arrays
real(8) rfmt(npmtmax)
complex(8) zrhomt(npmtmax),zvclmt(npmtmax)
! allocate global array
if (allocated(gvclns)) deallocate(gvclns)
allocate(gvclns(npmtmax,3,natmtot))
! loop over atoms
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  np=npmt(is)
! take the average of the static density over the three directions
  rfmt(1:np)=(1.d0/3.d0) &
   *(rhosmt(1:np,ias,1)+rhosmt(1:np,ias,2)+rhosmt(1:np,ias,3))
! convert to complex spherical harmonics expansion
  call rtozfmt(nr,nri,rfmt,zrhomt)
! solve Poisson's equation in the muffin-tin
  call zpotclmt(nr,nri,nrmtmax,rlmt(:,:,is),wprmt(:,:,is),zrhomt,zvclmt)
! add the nuclear Coulomb potential
  i=1
  do ir=1,nri
    zvclmt(i)=zvclmt(i)+vcln(ir,is)
    i=i+lmmaxi
  end do
  do ir=nri+1,nr
    zvclmt(i)=zvclmt(i)+vcln(ir,is)
    i=i+lmmaxo
  end do
! compute the gradient of the potential and store in global array
  call gradzfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),zvclmt,npmtmax, &
   gvclns(:,:,ias))
end do
end subroutine

