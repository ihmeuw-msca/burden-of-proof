# make file for the fortran module

.PHONY: clean

special_mat: special_mat.f90
	f2py -c -m special_mat special_mat.f90 -L/usr/local/lib -llapack -lblas

clean:
	find . -name '__pycache__' | xargs rm -rf
	find . -name '*.so' | xargs rm -rf
