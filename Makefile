cc=gfortran
flags=-O0 -g
include_dir=./include
source_dir=./src
all:
	$(cc) $(flags) -fopenmp -I$(include_dir) $(source_dir)/*.f90 $(source_dir)/*.f