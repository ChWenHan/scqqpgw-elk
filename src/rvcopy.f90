
! Copyright (C) 2022 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

pure subroutine rvcopy(n,x,y)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: x(n)
real(8), intent(out) :: y(n)
y(:)=x(:)
end subroutine
