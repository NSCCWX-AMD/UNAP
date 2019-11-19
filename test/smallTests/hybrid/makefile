# CPPC = icpc
# CC   = icc
# FC	 = ifort
# F_LINKER = ${FC}

CPPC = swg++453
FC   = mpif90
F_LINKER = swld453-fort


cppTest1.o: cppTest1.cpp
	${CPPC} -c -o $@ $^

#cppTest
cppMain: cppMain.o cppTest1.o
	${CPPC} $^ -o $@
cppMain.o: main.cpp
	${CPPC} -c -o $@ $^

#cTest
cMain: cMain.o cppTest1.o
	${CC} $^ -o $@
cMain.o: main.c
	${CC} -c -o $@ $^

#f90Test
f90Main: f90Main.o cppTest1.o
	${F_LINKER} $^ -o $@ -lstdc++
f90Main.o: main.f90
	${FC} -c -o $@ $^


clean:
	rm -f *.o

