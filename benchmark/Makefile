cc=gfortran
flags=-O3 -g -fopenmp
include_dir=../include
source_dir=../src

amxpby: csrmat blas
	$(cc) $(flags) -I$(include_dir)	$(source_dir)/blas.f90 \
		$(source_dir)/readmat.f90 \
		$(source_dir)/CSRMAT.f90 \
		$(source_dir)/irow2iat.f90 \
		$(source_dir)/*.f \
		amxpby.f90 \
		-o amxpby.out

blas: $(source_dir)/blas.f90 csrmat
	$(cc) $(flags) -c $<

csrmat: $(source_dir)/CSRMAT.f90 precision
	$(cc) $(flags) -c $<

precision: $(source_dir)/precision.f90
	$(cc) $(flags) -c $<

utils: $(source_dir)/utils.f90 precision
	$(cc) $(flags) -c $<

clean:
	rm *.mod *.o
