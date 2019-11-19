program ex12f

implicit none
INCLUDE 'mpif.h'

integer :: ierr

integer, parameter :: real_size = 8
integer, parameter :: int_size = 8

integer (kind = 8) :: APtr

integer :: precond = 3
integer :: aggl = 2
integer :: smoother = 2

integer :: NNZ, Ncel, upperSize

real (kind = real_size), allocatable :: x(:), b(:)
integer, allocatable :: Arow(:), Acol(:)
real (kind = real_size), allocatable :: Acoo(:), faceArea(:)

logical :: symm

integer :: num_iterations = 0, maxiter = 200, miniter = 1
real(kind = real_size) :: final_res_norm = 0.0, tol = 0, relTol = 1D-2

integer :: i, j, f, mype, ranks, iter
character(len = 200) :: filename
character(len = 5) :: str1
logical :: prun
character(len = 38):: fileDir =  "../test/exData/compass/cavityp4/step3/"

integer :: Nneigh, offdiag_size, locSize, locStart
integer, allocatable :: neighidx(:), interfaceidx(:), offdiag_row(:)
real (kind = real_size), allocatable :: offdiag_coeffs(:)

call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD, mype, ierr)
call mpi_comm_size(MPI_COMM_WORLD, ranks, ierr)

if(ranks > 1) then
	prun = .true.
else
	prun = .false.
endif

write(str1,"(i5.5)") mype

#ifdef SW_SLAVE
	call swacc_init();
#endif

f = 1
filename = fileDir//"A_p"//"_"//str1//".txt"
open(unit = f, file = filename, status='old')
read(f, *) Ncel, NNZ, symm
close(f)

upperSize = (NNZ-Ncel)/2

allocate(x(1:Ncel))
allocate(b(1:Ncel))
allocate(Arow(1:NNZ))
allocate(Acol(1:NNZ))
allocate(Acoo(1:NNZ))
allocate(faceArea(1:upperSize))

if(mype == 0) then
	write(*, *) "Start reading A"
end if
f = 1
filename = fileDir//"A_p"//"_"//str1//".txt"
open(unit = f, file = filename, status='old')
read(f, *)
do i=1, NNZ
	read(f, *) Arow(i), Acol(i), Acoo(i)
	! write(*, *) Arow(i), Acol(i), Acoo(i)
end do
close(f)
if(mype == 0) then
	write(*, *) "Finish reading A"
end if

if(mype == 0) then
	write(*, *) "Start reading b"
end if
filename = fileDir//"b_p"//"_"//str1//".txt"
open(unit = f, file = filename, status='old')
read(f, *)
do i=1, Ncel
	read(f, *) b(i)
	! write(*, *) j, b(i)
end do
close(f)
if(mype == 0) then
	write(*, *) "Finish reading b"
end if

if(mype == 0) then
	write(*, *) "Start reading faceArea"
end if
filename = fileDir//"faceArea_p"//"_"//str1//".txt"
open(unit = f, file = filename, status='old')
read(f, *)
do i=1, upperSize
	read(f, *) faceArea(i)
	! if(mype  == 1) then
	! 	write(*, *) i, faceArea(i)
	! end if
end do
close(f)
if(mype == 0) then
	write(*, *) "Finish reading faceArea"
end if

if(prun) then
	if(mype == 0) then
		write(*, *) "Start reading interfaces"
	end if
	filename = fileDir//"interfaces_p"//"_"//str1//".txt"
	open(unit = f, file = filename, status='old')
	read(f, *) Nneigh, offdiag_size
	close(f)

	allocate(neighidx(1:Nneigh))
	allocate(interfaceidx(1:Nneigh+1))
	allocate(offdiag_row(1:offdiag_size))
	allocate(offdiag_coeffs(1:offdiag_size))

	interfaceidx(1) = 1

	open(unit = f, file = filename, status='old')
	read(f, *)
	do i=1, Nneigh
		read(f, *) neighidx(i), locSize
		do j=1, locSize
			locStart = interfaceidx(i) + j - 1
			read(f,*) offdiag_row(locStart), offdiag_coeffs(locStart)
		end do
		interfaceidx(i+1) = interfaceidx(i) + locSize
	end do
	close(f)

	if(mype == 0) then
		write(*, *) "Finish reading interfaces"
	end if
end if

do iter=1, 5

	do i=1, Ncel
		x(i) = 0
	end do

	call swtimerstart("mg solve")

	call coo2ldumatrixcreat(APtr, Acoo, Arow, Acol, Ncel, NNZ, symm)

	if(prun) then
		call matrixinterfacescreat(APtr, Nneigh, neighidx, interfaceidx, offdiag_row, offdiag_coeffs)
	end if

	! call pcgsolversolve(x, APtr, b, Ncel, precond, tol, relTol, maxiter, miniter, num_iterations, final_res_norm)
	! call pbicgstabsolversolve(x, APtr, b, Ncel, precond, tol, relTol, maxiter, miniter, num_iterations, final_res_norm)
	call mgsolversolve(x, APtr, b, Ncel, aggl, smoother, tol, relTol, maxiter, miniter, num_iterations, final_res_norm, faceArea)

	call swtimerend("mg solve")

	if(mype == 0) then
		write(*,*) "Finish step ", iter
		write(*,*) "Num of iterations is ", num_iterations
		write(*,*) "Final residual is ", final_res_norm
	end if
end do

deallocate(x)
deallocate(b)
deallocate(Arow)
deallocate(Acol)
deallocate(Acoo)
deallocate(faceArea)

if(prun) then
	deallocate(neighidx)
	deallocate(interfaceidx)
 	deallocate(offdiag_row)
 	deallocate(offdiag_coeffs)
endif

call swtimerprintall()

# ifdef SW_SLAVE
	call swacc_end();
# endif

stop
end program ex12f
