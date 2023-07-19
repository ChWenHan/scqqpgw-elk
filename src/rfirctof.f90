
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfirctof(rfirc,rfir)
use modmain
implicit none
! arguments
real(8), intent(in) :: rfirc(ngtc)
real(8), intent(out) :: rfir(ngtot)
! local variables
integer ig
! allocatable arrays
complex(8), allocatable :: zfftc(:),zfft(:)
allocate(zfftc(ngtc),zfft(ngtot))
! Fourier transform function on coarse grid to G-space
zfftc(:)=rfirc(:)
call zfftifc(3,ngdgc,-1,zfftc)
! Fourier transform to fine real-space grid
zfft(:)=0.d0
do ig=1,ngvc
  zfft(igfft(ig))=zfftc(igfc(ig))
end do
call zfftifc(3,ngridg,1,zfft)
! output real function
rfir(:)=dble(zfft(:))
deallocate(zfftc,zfft)
end subroutine

