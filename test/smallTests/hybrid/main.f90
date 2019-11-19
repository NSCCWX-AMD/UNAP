program main
	implicit none
	integer (kind=8) :: ptr
	integer (kind=4) :: i = 10
	real(kind = 8) ::v = 0.9
	logical :: symm = .false.
	call cppFunc()
	call show1(ptr, i, v, symm)
	call show2(ptr)
end program main
