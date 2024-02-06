.POSIX:
.SUFFIXES:
.SUFFIXES: .o .f90

VERSION = 1.0
SONUM := $(shell echo $(VERSION) | cut -d '.' -f 1)

PREFIX = /usr/local
INCDIR = $(PREFIX)/include
LIBDIR = $(PREFIX)/lib
MANDIR = $(PREFIX)/share/man

AR = ar
FC = gfortran

# DISLINROOT = /usr/local/dislin
# DISLININC  = -I$(DISLINROOT)/gf/real64
# DISLINLIB  = -L$(DISLINROOT) -ldislin_d

NETCDFINC = -I/usr/include
NETCDFLIB = -lnetcdff

FFLAGS = -std=f2008 -O3 \
				 -fdefault-real-8 -fdefault-double-8 \
				 -ffree-form -fmax-errors=1 \
				 -pedantic -Wall \
				 -m64 -J./include \
				 $(NETCDFINC)
LDFLAGS = -s

SRC = $(wildcard src/*.f90)
OBJ = $(SRC:.f90=.o)

%.o: %.f90
	@echo FC $<
	@$(FC) -o $@ -c $(FFLAGS) $<

all: libstabpack

libstabpack: $(OBJ)
	@echo AR $(@).a
	@$(AR) rcs $(@).a $^
	@echo LD $(@).so
	@$(FC) -fPIC -shared -o $(@).so.$(VERSION) $^
	@ln -s $(@).so.$(VERSION) $(@).so.$(SONUM)
	@ln -s $(@).so.$(SONUM) $(@).so

help: doc/bez3.scd
	scdoc < $< > bez.3

clean:
	find . \( -name '*.mod' -o -name '*.o' \) -exec rm {} \;
	rm -f libstabpack.*

install:
	@echo installing

.PHONY: all clean install
