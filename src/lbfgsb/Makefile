AR  = ar
FC  = gfortran

FFLAGS = -std=legacy -O3 \
				 -ffixed-form -fmax-errors=1 \
				 -pedantic -Wall
LDFLAGS = -lblas -llapack

%.o: %.f
	@echo FC $<
	@${FC} -o $@ -c ${FFLAGS} $<

all: lbfgsb

lbfgsb: lbfgsb.o
	@echo AR $@.a
	@${AR} rcs $@.a $<

clean:
	rm -f *.o *.a

.PHONY: clean
