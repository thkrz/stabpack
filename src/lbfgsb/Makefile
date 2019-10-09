AR  = ar
FC  = gfortran

MKLROOT = /opt/intel/mkl
MKLINC  = -I${MKLROOT}/include
MKLLIB  = -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl

FFLAGS = -std=legacy -O3 \
				 -ffixed-form -fmax-errors=1 \
				 -pedantic -Wall \
				 ${MKLINC}
LDFLAGS = ${MKLLIB}

all: lbfgsb

lbfgsb: lbfgsb.f
	@echo FC $<
	@${FC} ${FFLAGS} -c $<
