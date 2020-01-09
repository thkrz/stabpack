include config.mk

CMD = #debkin seep stab
STABPACK = $(wildcard src/stabpack/*.f90)
LIBSSAT = $(wildcard src/libssat/*.f90)
LDFLAGS += -L./ -lssat -lstabpack ${NETCDFLIB}

%.o: %.f90
	@echo FC $<
	@${FC} -o $@ -c ${FFLAGS} $<

all: stabpack libssat ${CMD}

lbfgsb:
	@${MAKE} -C src/$@

stabpack: ${STABPACK:.f90=.o}
	@echo AR lib${@}.a
	@${AR} rcs lib${@}.a $^

libssat: ${LIBSSAT:.f90=.o}
	@echo AR ${@}.a
	@${AR} rcs ${@}.a $^

.SECONDEXPANSION:
${CMD}: src/cmd/$$@.o
	@echo LD $@
	@${FC} -o $@ $< ${LDFLAGS}

test: hyp2f1

hyp2f1: src/stabpack/num.o test/hyp2f1.o
	@echo LD $@
	@${FC} -o test/$@ $^ ${LDFLAGS}

clean:
	@echo cleaning...
	@find . \( -name '*.mod' -o -name '*.o' \) -exec rm {} \;
	@rm -f libstabpack.a libssat.a ${CMD}

install:
	@echo installing

.PHONY: all clean install
