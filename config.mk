CXX=mpicxx
DEFS=-DMKL -DCBLAS -DLAPACKE
INCLUDES=-I/opt/intel/compilers_and_libraries_2018.2.199/linux/mkl/include
CXXFLAGS=-g -O2 $(DEFS) -std=c++0x -mkl -fPIC $(INCLUDES)
LDFLAGS= 
