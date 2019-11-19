program ex9f
implicit none
	INCLUDE 'mpif.h'
	integer  ierr

	integer, parameter :: real_size = 8

	integer (kind = 8) :: APtr

	integer :: precond = 3
	integer :: aggl = 1
	integer :: smoother = 1

	integer, parameter :: NNZ = 25338, Ncel = 5150

	real (kind = real_size), allocatable :: x(:), b(:)
	integer, allocatable :: Arow(:), Acol(:)
	real (kind = real_size), allocatable :: Acoo(:), faceArea(:)

	logical :: symm = .false.

	integer :: num_iterations = 0, maxiter = 1000, miniter = 1
	real(kind = real_size) :: final_res_norm = 0.0, tol = 0, relTol = 1D-8

	integer :: i, j, f

	call mpi_init(ierr)

	allocate(x(1:Ncel))
	allocate(b(1:Ncel))
	allocate(Arow(1:NNZ))
	allocate(Acol(1:NNZ))
	allocate(Acoo(1:NNZ))
	allocate(faceArea(1:NNZ))


	write(*, *) "Start reading A"
	f = 1
	open (unit = f, file = "../test/exData/702NaVi2X/matrix.m", status='old')
	! read(f, *)
	do i=1, NNZ
		read(f, *) Arow(i), Acol(i), Acoo(i)
		! write(*, *) Arow(i), Acol(i), Acoo(i)
	end do
	close(f)
	write(*, *) "Finish reading A"


	write(*, *) "Start reading b"
	f = 2
	open (unit = f, file = "../test/exData/702NaVi2X/vector.m", status='old')
	! read(f, *)
	do i=1, Ncel
		! read(f, *) j, b(i)
		read(f, *) b(i)
		! write(*, *) j, b(i)
	end do
	close(f)
	write(*, *) "Finish reading b"

	! do i=1,


	do i=1, Ncel
		x(i) = 0
	end do

	call coo2ldumatrixcreat(APtr, Acoo, Arow, Acol, Ncel, NNZ, symm)
	write(*, *) "Finish creating lduMatrix"

	! call pcgsolversolve(x, APtr, b, Ncel, precond, tol, relTol, maxiter, miniter, num_iterations, final_res_norm)
	! call pbicgstabsolversolve(x, APtr, b, Ncel, precond, tol, relTol, maxiter, miniter, num_iterations, final_res_norm)
	call mgsolversolve(x, APtr, b, Ncel, aggl, smoother, tol, relTol, maxiter, miniter, num_iterations, final_res_norm, b)

	write(*,*) "Num of iterations is ", num_iterations
	write(*,*) "Final residual is ", final_res_norm

stop
end program ex9f