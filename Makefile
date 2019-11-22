include config.mk

CMD = seep stab
STABPACK = src/bempack.o src/mospack.o src/numpack.o \
	src/statpack.o src/wapack.o
OBJ = src/ssat_env.o src/scx.o

%.o: %.f90
	@echo FC $<
	@${FC} -o $@ -c ${FFLAGS} $<

all: ${CMD}

lbfgsb:
	@${MAKE} -C src/$@

seep: ${STABPACK} ${OBJ} src/seep.o
	@echo LD $@
	@${FC} -o $@ $^ ${LDFLAGS}

stab: ${STABPACK} ${OBJ} src/stab.o
	@echo LD $@
	@${FC} -o $@ $^ ${LDFLAGS}

test: hyp2f1

hyp2f1: src/numpack.o test/hyp2f1.o
	@echo LD $@
	@${FC} -o test/$@ $^ ${LDFLAGS}

clean:
	@echo cleaning...
	@find . \( -name '*.mod' -o -name '*.o' \) -exec rm {} \;
	@rm -f ${CMD}

.PHONY: all clean
