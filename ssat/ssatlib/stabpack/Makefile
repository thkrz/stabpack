SUFFIX = $(shell python3-config --extension-suffix)
SRC = fmin.f90
OBJ = ${SRC:.f90=${SUFFIX}}

%${SUFFIX}: %.f90
	@echo FORT $<
	@f2py3 --quiet -c -DNPY_NO_DEPRECATED_API=13 -m $* $<

all: ${OBJ}

clean:
	@rm -f *.mod ${OBJ}

.PHONY: clean
