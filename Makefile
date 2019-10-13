include config.mk

OBJ = src/mospack.o src/numpack.o src/wapack.o \
      src/ssat_env.o src/scx.o src/wasim.o \
      src/main.o

%.o: %.f90
	@echo FC $<
	@${FC} -o $@ -c ${FFLAGS} $<

all: ssat

lbfgsb:
	@${MAKE} -C src/$@

ssat: ${OBJ}
	@echo LD $@
	@${FC} -o $@ $^ ${LDFLAGS}

test: hyp2f1

hyp2f1: src/numpack.o test/hyp2f1.o
	@echo LD $@
	@${FC} -o test/$@ $^ ${LDFLAGS}

clean:
	@echo cleaning...
	@find . \( -name '*.mod' -o -name '*.o' \) -exec rm {} \;
	@rm -f ssat

.PHONY: all clean
