
! Copyright (C) 2016 A. Davydov, A. Sanna, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

module modgw

! maximum Matsubara frequency for the GW calculation
real(8) wmaxgw
! maximum number of Matsubara frequencies
integer nwgw
! integer grid intervals for Matsubara frequencies
integer intwgw(2)
! map from frequency index to FFT array
integer, allocatable :: iwfft(:)
!temp evalsv for constructing new G
real(8), allocatable :: evalsv0(:,:)
! maximum fermionic Matsubara frequency index to be used for the GW calculation
integer nwfm
! maximum bosonic frequency index
integer nwbs
! number of real points W of analytical continuation
integer nwplotgw
! energy distance of AC
real(8) dwgw
! dvse: residu of V_correlation_se
real(8), allocatable :: dvse(:,:),dvsediag(:,:)
! alphagw. mixing ration of SE
real(8) alphagw
!
integer(4) fscgw
!
logical mixse
!
logical calcgw
!if renomrlaize the qp by factor of (1-(d(se(ei))/dw))^-1
logical qprenorm
!
logical readepsinv
!if perform hf or dft gw
logical hfgw
!if write seiw
logical writeseiw
!if the occupancy from dyson equation for Se_iw
logical dmatfromseiw
!real(8) maxdvse(3)
! dw_sp=minval(abs(evalsv(:,ik)-efermi))
!real(8) dw_sp
! time spent for AC
real(8) timeac
!
real(8), allocatable :: vsemt(:,:)
real(8), allocatable :: vseir(:)
real(8), allocatable :: bsemt(:,:)
real(8), allocatable :: bseir(:)
!
logical tvse
!

! energy range of ac
real(8) wplotgw(2)
! imaginary frequencies used for the GW calculation
real(8), allocatable :: wgw(:)
! complex fermionic frequencies
complex(8), allocatable :: wfm(:)
! index mapping w= evalsv_i
integer, allocatable :: seidx(:,:)
! twdiag is .true. if the screened interaction W is taken to be diagonal
logical twdiag
! tsediag is .true. if the GW self-energy is taken to be diagonal
logical tsediag
! whehter diagonalise interacting green's function
logical gdiag
! whether half scrho
integer sclevel
! whether perform dyson in gwdmatk2
logical dysongwdmatk2
! type of analytic continuation to be used for determining the self-energy on
! the real axis
integer actype
! number of poles used for fitting the self-energy matrix elements
integer npole
! number of complex shifts used in averaging the Pade approximant for the
! analytic continuation of the self-energy to the real axis
integer nspade
!
integer maxsclgw
!
real(8) tempkgw


end module

