! test for reading data from a text file

program main
	implicit none
	integer :: f
	integer :: NNZ, Ncel
	integer :: i, j
	integer :: row, col
	real    :: val

	NNZ = 2712000
	Ncel = 400000

	f = 1
	open (unit = f, file = "A.txt", status='old')
	read(f, *)
	do i=1, NNZ
		read(1, *) row, col, val
		write(*,*) row, col, val
	end do
	close(f)
	stop
end program main