ROOT = ${PWD}

default:
	cd ./src && make all -j 12 PLA=sw
	# cp ./build/lib/libswunap.a -f ../../../../swlib/libunap.a

clean:
	cd ${ROOT}/src  && make clean
	cd ${ROOT}/test && make clean

link:
	cd ${ROOT}/src  && make link

cleanlink:
	rm -rf build/lnInclude

cleanall:
	rm -rf build
	rm -f test/*.o
