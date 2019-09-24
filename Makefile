include config.mk

OBJ = src/slope.o src/main.o

%.o: %.f90
	@echo FC $<
	@${FC} -o $@ -c ${FFLAGS} $<

all: ssat

ssat: ${OBJ}
	@echo LD $@
	@${FC} -o $@ ${OBJ} ${LDFLAGS}

clean:
	find . \( -name '*.mod' -o -name '*.o' \) -exec rm {} \;
	rm -f ssat

.PHONY: all clean
