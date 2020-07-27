CXX=mpicxx
BLAS=CBLAS
LAPACK=LAPACKE
INCLUDES=-I$(HOME)/critter/include -I/opt/intel/compilers_and_libraries_2018.2.199/linux/mkl/include
CXXFLAGS=-g -O2 -DMKL -D$(BLAS) -D$(LAPACK) -std=c++0x -mkl -fPIC
LDFLAGS= 
