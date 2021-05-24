cc=gfortran
flags=-O3 -g -fopenmp
include_dir=./include
source_dir=./src

all: csrmat utils blas
	$(cc) $(flags) -I$(include_dir) $(source_dir)/*.f90 $(source_dir)/*.f

blas: $(source_dir)/blas.f90 csrmat
	$(cc) $(flags) -c $< -o $(source_dir)/class_blas.mod

csrmat: $(source_dir)/CSRMAT.f90 precision
	$(cc) $(flags) -c $< -o $(source_dir)/class_csrmat.mod

precision: $(source_dir)/precision.f90
	$(cc) $(flags) -c $< -o $(source_dir)/class_precision.mod

utils: $(source_dir)/utils.f90 precision
	$(cc) $(flags) -c $< -o $(source_dir)/utils.mod

clean:
	rm *.mod *.o
