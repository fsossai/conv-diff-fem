cc=gfortran
flags=-O0 -g -fopenmp
include_dir=./include
source_dir=./src

all: csrmat utils
	$(cc) $(flags) -I$(include_dir) $(source_dir)/*.f90 $(source_dir)/*.f

csrmat: precision blas $(source_dir)/CSRMAT.f90
	$(cc) $(flags) -c $(source_dir)/CSRMAT.f90

precision: $(source_dir)/precision.f90
	$(cc) $(flags) -c $^

utils: $(source_dir)/utils.f90
	$(cc) $(flags) -c $^

blas: $(source_dir)/blas.f90
	$(cc) $(flags) -c $^

clean:
	rm *.mod *.o
