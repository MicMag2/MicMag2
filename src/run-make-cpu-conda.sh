#!/bin/bash
thisDir=`pwd`
if [[ `basename $thisDir` != "src" ]]; then
	echo "current working directory is NOT a src directory. exiting for safety reasons."
else
	mkdir build 
	cd build
	echo "clearing build directory..."
	make clean 
	rm CMakeCache.txt 
	rm config.h
	
	echo "starting compiler..."
	cmake .. -DPYTHON_LIBRARY=$(which python3) -DPYTHON_INCLUDE_DIR=${CONDA_PREFIX}/include/python3.8 -DCMAKE_CXX_FLAGS=-isystem\ ${CONDA_PREFIX}/include 
	make $@ 
	#cmake .. && make $@ 
	make $@  && echo "CPU compiling successful."
	cp magneto_cpu.py ../magnum/magneto_cpu.py 2>/dev/null
	cp _magneto_cpu.so ../magnum/_magneto_cpu.so 2>/dev/null
	cp magneto_cuda.py ../magnum/magneto_cuda.py 2>/dev/null
	cp _magneto_cuda.so ../magnum/_magneto_cuda.so 2>/dev/null
        cd ..
	export PYTHONPATH=$(pwd)
	echo "python path exported temporarilly. Consider adding it in you .bashrc"
fi
