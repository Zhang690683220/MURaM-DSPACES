.SILENT:

INCDIR = ../include
LIBDIR = ../lib
CFDLIB = $(LIBDIR)/libcfd.a
LIBS   = -L$(LIBDIR) -lcfd -lm

OPT = -O2

CCC = mpic++
CCFLAGS = $(OPT) -I$(INCDIR)

SHELL = /bin/sh

TEMPLATEDSOURCES = \
	euler_lf.C \
	euler_roe.C \
	euler_hlle.C \
	euler_hllel.C \
	euler_godunov.C \
	mhd_lf.C \
	mhd_hlle.C \
	mhd_hllel.C \
	mhd_boltzmann.C \
	mhd_roe.C

COMMONOBJECTS = \
	euler_common.o \
	mhd_common.o

library: $(CFDLIB)
	echo "$(CFDLIB) is up to date."

$(CFDLIB): objects $(COMMONOBJECTS)
	 ar cr $(CFDLIB) *.o
	-ranlib $(CFDLIB)
	 echo "Created $(CFDLIB)"

objects: $(TEMPLATEDSOURCES)
	 touch objects
	 for file in $?; do \
	   echo "Compiling $${file} using float"; \
	   $(CCC) $(CCFLAGS) -DPRECISION=float  -o $${file%.C}_f.o -c $${file}; \
	   echo "Compiling $${file} using double"; \
	   $(CCC) $(CCFLAGS) -DPRECISION=double -o $${file%.C}_d.o -c $${file}; \
	 done

install: $(CFDLIB)
	-mv -f $(CFDLIB) $(LIBDIR)
	 echo "Moved $(CFDLIB) to $(LIBDIR)"

all: library install

clean:
	-rm -rf *~ *.a *.o *.i *.ii *.out ii_files core objects

distclean: clean
	-rm -f $(LIBDIR)/$(CFDLIB)

######### Suffix rules ########################################

.SUFFIXES : .o .cpp .c++ .C
.cpp.o:
	echo "Compiling $<"
	$(CCC) $(CCFLAGS) -c $<
.c++.o:
	echo "Compiling $<"
	$(CCC) $(CCFLAGS) -c $<
.C.o:
	echo "Compiling $<"
	$(CCC) $(CCFLAGS) -c $<
