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
	if nvcc --version 2>/dev/null; then
		echo ""
		echo "CUDA found. trying to compile for CPU and GPU..."
		# double compiling necessary on some systems due to issue in nvcc compiler. cannot hurt, so just doing it by default...
		cmake .. -DPYTHON_LIBRARY=$(which python3)  
		cmake .. && make $@
		make $@ echo "CPU compiling successful."
		cp magneto_cpu.py ../magnum/magneto_cpu.py 2>/dev/null
		cp _magneto_cpu.so ../magnum/_magneto_cpu.so 2>/dev/null
		echo "trying GPU compiling..." && rm CMakeCache.txt && cmake .. -DENABLE_CUDA_64=on
		make $@ 
		make $@ && echo "CPU and GPU compiling successful."
		cp magneto_cpu.py ../magnum/magneto_cpu.py 2>/dev/null
		cp _magneto_cpu.so ../magnum/_magneto_cpu.so 2>/dev/null
		cp magneto_cuda.py ../magnum/magneto_cuda.py 2>/dev/null
		cp _magneto_cuda.so ../magnum/_magneto_cuda.so 2>/dev/null
	else
		echo "No CUDA found. compiling for CPU only..."
		cmake .. -DPYTHON_LIBRARY=$(which python3) 
		cmake .. && make $@ 
		make $@  && echo "CPU compiling successful."
		cp magneto_cpu.py ../magnum/magneto_cpu.py 2>/dev/null
		cp _magneto_cpu.so ../magnum/_magneto_cpu.so 2>/dev/null
		cp magneto_cuda.py ../magnum/magneto_cuda.py 2>/dev/null
		cp _magneto_cuda.so ../magnum/_magneto_cuda.so 2>/dev/null
	fi
fi
