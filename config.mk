AR  = ar
CAF = caf
FC  = gfortran

# DISLINROOT = /usr/local/dislin
# DISLININC  = -I${DISLINROOT}/gf/real64
# DISLINLIB  = -L${DISLINROOT} -ldislin_d

# MKLROOT = /opt/intel/mkl
# MKLINC  = -I${MKLROOT}/include
# MKLLIB  = -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl

NETCDFINC = -I/usr/include
NETCDFLIB = -lnetcdff

FFLAGS = -std=f2008 -O3 \
				 -fdefault-real-8 -fdefault-double-8 \
				 -ffree-form -fmax-errors=1 \
				 -pedantic -Wall \
				 -m64 -J./src \
				 ${NETCDFINC}
LDFLAGS = -s
