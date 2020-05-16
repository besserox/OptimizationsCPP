all: main main-omp main-acc main-acc-cpu main-acc-cpu-omp main-cuda

main: main.cpp
	g++ $^ -o $@ -O3

main-omp: main-omp.cpp
	g++ $^ -o $@ -O3 -fopenmp

main-acc: main-acc.cpp
	pgc++ $^ -o $@ -O3 -acc -Minfo=accel

main-acc-cpu: main-acc.cpp
	pgc++ $^ -o $@ -O3

main-acc-cpu-omp: main-acc.cpp
	pgc++ $^ -o $@ -O3 -mp

main-cuda: main-cuda.cpp
	nvcc -x cu $^ -o $@ -O3
