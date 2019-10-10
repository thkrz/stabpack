include config.mk

OBJ = src/aquapack.o src/numpack.o src/stabpack.o \
      src/ssat_env.o src/scx.o src/main.o

%.o: %.f90
	@echo FC $<
	@${FC} -o $@ -c ${FFLAGS} $<

all: ssat

lbfgsb:
	@${MAKE} -C src/$@

ssat: ${OBJ}
	@echo LD $@
	@${FC} -o $@ ${OBJ} ${LDFLAGS}

clean:
	@echo cleaning...
	@find . \( -name '*.mod' -o -name '*.o' \) -exec rm {} \;
	@rm -f ssat

.PHONY: all clean
