cc=gfortran
flags=-Wall -O3 -g -fopenmp
include_dir=./include
source_dir=./src

all: csrmat fem utils blas
	$(cc) $(flags) -I$(include_dir) $(source_dir)/*.f90 $(source_dir)/*.f -o solver.out

fem: $(source_dir)/fem.f90 csrmat blas utils
	$(cc) $(flags) -I$(include_dir) -c $<

blas: $(source_dir)/blas.f90 csrmat
	$(cc) $(flags) -c $<

csrmat: $(source_dir)/CSRMAT.f90 precision
	$(cc) $(flags) -c $<

utils: $(source_dir)/utils.f90 precision
	$(cc) $(flags) -c $<

precision: $(source_dir)/precision.f90
	$(cc) $(flags) -c $<

clean:
	rm *.mod *.o
