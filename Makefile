include config.mk

CMD = debkin seep stab
STABPACK = src/stabpack/bem.o src/stabpack/mos.o src/stabpack/num.o \
	src/stabpack/stat.o src/stabpack/wa.o
OBJ = src/ssat_env.o src/scx.o

%.o: %.f90
	@echo FC $<
	@${FC} -o $@ -c ${FFLAGS} $<

all: ${CMD}

lbfgsb:
	@${MAKE} -C src/$@

debkin: ${STABPACK} ${OBJ} src/debkin.o
	@echo LD $@
	@${FC} -o $@ $^ ${LDFLAGS}

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
