
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure complex(8) function zfmtint(nr,nri,wr,zfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: wr(nr)
complex(8), intent(in) :: zfmt(*)
! local variables
integer ir,i
zfmtint=0.d0
i=1
do ir=1,nri
  zfmtint=zfmtint+wr(ir)*zfmt(i)
  i=i+lmmaxi
end do
do ir=nri+1,nr
  zfmtint=zfmtint+wr(ir)*zfmt(i)
  i=i+lmmaxo
end do
zfmtint=fourpi*y00*zfmtint
end function

