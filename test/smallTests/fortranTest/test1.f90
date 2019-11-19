! test for sin(pi)
! results shows that sin(pi) can not close to 0
! different from the one in c++/c

program main
	implicit none

	real (kind=8), parameter:: pi = 3.1415926536

	write(*,*) sin(pi)
	stop
end program main