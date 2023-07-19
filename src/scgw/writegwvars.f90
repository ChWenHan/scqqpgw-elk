subroutine writegwvars
use modmain
use modgw
use modmpi
implicit none
! arguments
! local variables
integer it1,it2
!real(8) dw,w,e
character(256) fname
if (mp_mpi) then
  write(fname,'("GW_VAR1.OUT")') 
  open(50,file=trim(fname),form='FORMATTED')
  write(50,'(A10,5X,G18.10)') 'wmaxgw', wmaxgw
  write(50,'(A10,5X,G18.10)') 'nwgw', nwgw
  write(50,'(A10,5X,G18.10)') 'intwgw(1)', intwgw(1)
  write(50,'(A10,5X,G18.10)') 'intwgw(2)', intwgw(2)
  write(50,'(A10,5X,G18.10)') 'nwfm', nwfm
  write(50,'(A10,5X,G18.10)') 'nwbs', nwbs
  write(50,'(A10,5X,G18.10)') 'ngrf', ngrf
  write(50,'(A10,5X,G18.10)') 'nwrf', nwrf
  close(50)
  !whch---- modified
  write(fname,'("GW_iwfft.OUT")') 
  open(51,file=trim(fname),form='FORMATTED')
  do it1=intwgw(1),intwgw(2)
    write(51,'(2I10.8)') it1,iwfft(it1)
  end do
  close(51)

  write(fname,'("GW_wgw.OUT")') 
  open(52,file=trim(fname),form='FORMATTED')
  do it1=intwgw(1),intwgw(2)
    write(52,'(I10.8,G18.10)') it1, wgw(it1)
  end do
  close(52)

  write(fname,'("GW_wfm.OUT")') 
  open(53,file=trim(fname),form='FORMATTED')
  do it1=0,nwfm
    write(53,'(2G18.10)') wfm(it1)
  end do
  close(53)
end if
return
end subroutine