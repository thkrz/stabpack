include config.mk

PACK = src/bempack.o src/mospack.o src/numpack.o \
	src/statpack.o src/wapack.o

%.o: %.f90
	@echo FC $<
	@${FC} -o $@ -c ${FFLAGS} $<

all: seep stab

lbfgsb:
	@${MAKE} -C src/$@

seep: ${PACK} src/ssat_env.o src/scx.o src/seep.o
	@echo LD $@
	@${FC} -o $@ $^ ${LDFLAGS}

stab: ${PACK} src/ssat_env.o src/scx.o src/stab.o
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
