include config.mk

CMD = debkin seep stab
STABPACK = src/stabpack/bem.o src/stabpack/mos.o src/stabpack/num.o \
	src/stabpack/stat.o src/stabpack/wa.o
OBJ = src/ssat_env.o src/scx.o
LDFLAGS += -L./ -lstabpack ${NETCDFLIB}

%.o: %.f90
	@echo FC $<
	@${FC} -o $@ -c ${FFLAGS} $<

all: stabpack ${CMD}

lbfgsb:
	@${MAKE} -C src/$@

stabpack: ${STABPACK}
	@echo AR lib${@}.a
	@${AR} rcs lib${@}.a $^

debkin: ${OBJ} src/debkin.o
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
	@rm -f libstabpack.a ${CMD}

.PHONY: all clean
