include config.mk

OBJ = src/ssat_env.o src/bezier.o src/fmin.o \
			src/interp1d.o src/rtfind.o src/slope.o \
			src/main.o

%.o: %.f90
	@echo FC $<
	@${FC} -o $@ -c ${FFLAGS} $<

all: ssat

ssat: ${OBJ}
	@echo LD $@
	@${FC} -o $@ ${OBJ} ${LDFLAGS}

clean:
	@echo cleaning...
	@find . \( -name '*.mod' -o -name '*.o' \) -exec rm {} \;
	@rm -f ssat

.PHONY: all clean
