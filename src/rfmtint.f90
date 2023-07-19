
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure real(8) function rfmtint(nr,nri,wr,rfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: wr(nr),rfmt(*)
! local variables
integer ir,i
rfmtint=0.d0
i=1
do ir=1,nri
  rfmtint=rfmtint+wr(ir)*rfmt(i)
  i=i+lmmaxi
end do
do ir=nri+1,nr
  rfmtint=rfmtint+wr(ir)*rfmt(i)
  i=i+lmmaxo
end do
rfmtint=fourpi*y00*rfmtint
end function

