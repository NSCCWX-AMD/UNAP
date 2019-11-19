! test for Fortran interface
! solver can use CG or MG

subroutine au(auf, u, rhx2, rhy2, i, j, N)
implicit none
	integer :: N
	integer, parameter :: real_size = 8
	real (kind=real_size) :: auf
	real (kind=real_size) :: u(0:(N+1)*(N+1)-1)
	real (kind=real_size) :: rhx2
	real (kind=real_size) :: rhy2
	integer :: i
	integer :: j


	auf = ( &
			rhx2*(u((i-1)*(N+1) + j  ) - 2*u(i*(N+1) + j) + u((i+1)*(N+1) + j  )) + &
			rhy2*(u(    i*(N+1) + j-1) - 2*u(i*(N+1) + j) + u(    i*(N+1) + j+1)) &
		  )
end subroutine au

! program main
program ex6f
implicit none

	integer, parameter :: real_size = 8

	integer (kind=8) :: APtr

	integer :: precond = 4
	integer :: aggl = 1
	integer :: smoother = 1
	integer :: solver

	integer (kind=4), parameter :: N = 1000
	real (kind=real_size), parameter :: pi = 3.1415926536

	integer, parameter :: nCells = (N-1)*(N-1)
	real (kind=real_size) :: x(0:nCells-1) = (/ (0.0, i=0, nCells-1) /)

	! domain solved
	real (kind=real_size), parameter :: x0 = 0.0;
	real (kind=real_size), parameter :: x1 = 1.0;
	real (kind=real_size), parameter :: y0 = 0.0;
	real (kind=real_size), parameter :: y1 = 1.0;

	real (kind=real_size), parameter :: dx = (x1 - x0) / N;
	real (kind=real_size), parameter :: dy = (y1 - y0) / N;

	real (kind=real_size) :: bAll(0:(N+1)*(N+1)-1)
	real (kind=real_size) :: xAll(0:(N+1)*(N+1)-1)
	real (kind=real_size) :: b(0:nCells-1)

	real (kind=real_size), parameter :: rhx2 = 1.0/dx/dx;
	real (kind=real_size), parameter :: rhy2 = 1.0/dy/dy;

	integer :: nZeros = 0
	integer :: nZerosCells(0:nCells)

	integer :: i, j
	integer :: n0, n1
	integer :: nStart = 0

	integer, allocatable :: lowerAddr(:)
	integer, allocatable :: upperAddr(:)

	real (kind=real_size), allocatable :: upper(:)
	! real, pointer :: lower(:)
	real (kind=real_size) :: diag(0:nCells-1) = (/ (-4*rhx2, i=0, nCells-1) /)

	real (kind=real_size) :: ax

	integer :: num_iterations = 0, maxiter = 100, miniter = 1
	real(kind=real_size) :: final_res_norm = 0.0, tol = 1D-10, relTol = 1D-4

	do i=0, N
		do j=0, N
			bAll(i*(N+1)+j) = -8.0*pi*pi*sin(2*pi*i*dx)*sin(2*pi*j*dy)
		end do
	end do

	! write(*,*) 2.0*pi
	! write(*,*) sin(2.0*pi)

	do i=0, N
		do j=0, N
			if ((i==0) .or. (i==N) .or. (j==0) .or. (j==N)) then
				xAll(i*(N+1)+j) = sin(2*pi*i*dx)*sin(2*pi*j*dy)
			else
				xAll(i*(N+1)+j) = 0.0
			end if
		end do
	end do

	! do i=0, ((N+1)*(N+1)-1)
	! 	print '("xAll: at i = " i4, ", xAll = " (e12.5))', i, xAll(i)
	! end do

	do i=1, N-1
		do j=1, N-1
			call au(ax, xAll, rhx2, rhy2, i, j, N)
			b((i-1)*(N-1) + j -1) = bAll(i*(N+1)+j) - ax
		end do
	end do

	do i=1, N-1
		do j=1, N-1
			n0 = nZeros
			if((i+1)<N) then
				nZeros = nZeros+1
			end if
			if((j+1)<N) then
				nZeros = nZeros+1
			end if
			n1 = nZeros
			nZerosCells((i-1)*(N-1) + j - 1 + 1) = n1 - n0 + nZerosCells((i-1)*(N-1) + j - 1)
		end do
	end do

	print '("nZeros in upper is ", i)', nZeros

	allocate(lowerAddr(0:nZeros-1))
	allocate(upperAddr(0:nZeros-1))

	do i=1, N-1
		do j=1, N-1
			if((j+1)<N) then
				lowerAddr(nStart) = (i-1)*(N-1) + j - 1
				upperAddr(nStart) = (i-1)*(N-1) + j - 1 + 1
				nStart = nStart + 1
			end if

			if((i+1)<N) then
				lowerAddr(nStart) = (i-1)*(N-1) + j - 1
				upperAddr(nStart) = (i-1)*(N-1) + j - 1 + N - 1
				nStart = nStart + 1
			end if
		end do
	end do

	allocate(upper(0:nZeros-1))
	! lower => upper

	do i=0, nZeros-1
		upper(i) = rhx2
	end do

	! do i=0, nZeros-1
	! 	print '("At i = " i4, ", lAddr = " i4 ", uAddr = ", i4)', i, lowerAddr(i), upperAddr(i)
	! end do

	! do i=0, nZeros-1
	! 	print '("At i = " i4 ", u = " (e12.5))', i, upper(i)
	! end do

	! do i=0, nCells-1
	! 	print '("At i = " i4 ", diag = " (e12.5) ", b = " (e12.5))', i, diag(i), b(i)
	! end do

	! call lduMatrixCreat(Aptr, nCells, nZeros, lowerAddr, upperAddr, upper, diag, upper)
	! call PCGSolverSolve(x, Aptr, b, nCells, precond, tol, relTol, maxiter, miniter, num_iterations, final_res_norm)

	! call lduMatrixCreat(Aptr, nCells, nZeros, lowerAddr, upperAddr, upper, diag, upper)
	! call PBiCGStabSolverSolve(x, Aptr, b, nCells, precond, tol, relTol, maxiter, miniter, num_iterations, final_res_norm)

	call lduMatrixCreat(Aptr, nCells, nZeros, lowerAddr, upperAddr, upper, diag, upper)
	call MGSolverSolve(x, Aptr, b, nCells, aggl, smoother, tol, relTol, maxiter, miniter, num_iterations, final_res_norm)



stop
end program ex6f