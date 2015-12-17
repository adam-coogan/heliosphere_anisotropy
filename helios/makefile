ifeq ($(shell hostname),hyades.ucsc.edu)
	CC=icpc -std=c++11 -I/pfs/sw/serial/gcc/boost-1.57.0/include # Broken?
else
	CC=clang++ -std=c++11 -I /usr/local/opt/boost/include
endif

CFLAGS=-c -O2
LFLAGS=-lboost_program_options
OBJS=$(BUILDDIR)/Parameters.o $(BUILDDIR)/Basic3DParams.o
BINDIR=bin
BUILDDIR=build
SRCDIR=src

all: $(BINDIR)/Helios # $(BINDIR)/HeliosTest $(BINDIR)/StraussUnitsModified

# Messy program using Strauss' units

#$(BINDIR)/StraussUnitsModified: $(BUILDDIR)/StraussUnitsModified.o
#	$(CC) $(LFLAGS) $(BUILDDIR)/StraussUnitsModified.o -o $(BINDIR)/StraussUnitsModified
#
#$(BUILDDIR)/StraussUnitsModified.o: $(SRCDIR)/StraussUnitsModified.cpp
#	$(CC) $(CFLAGS) $< -o $@
#
## Test program
#
#$(BINDIR)/HeliosTest: $(OBJS) $(BUILDDIR)/HeliosTest.o $(SRCDIR)/Basic3D.h
#	$(CC) $(LFLAGS) $(OBJS) $(BUILDDIR)/HeliosTest.o -o $(BINDIR)/HeliosTest
#
#$(BUILDDIR)/HeliosTest.o: $(SRCDIR)/HeliosTest.cpp $(SRCDIR)/Basic3D.h
#	$(CC) $(CFLAGS) $< -o $@

# Main program

$(BINDIR)/Helios: $(OBJS) $(BUILDDIR)/Helios.o $(SRCDIR)/TrajectoryBase.h $(SRCDIR)/Simulation.h
	$(CC) $(LFLAGS) $(OBJS) $(BUILDDIR)/Helios.o -o $(BINDIR)/Helios

$(BUILDDIR)/Helios.o: $(SRCDIR)/Helios.cpp $(SRCDIR)/Basic3D.h $(SRCDIR)/Basic3DParams.h
	$(CC) $(CFLAGS) $< -o $@

# Objects

$(BUILDDIR)/Basic3DParams.o: $(SRCDIR)/Parameters.h $(SRCDIR)/Basic3DParams.h $(SRCDIR)/Basic3DParams.cpp
	$(CC) $(CFLAGS) $(SRCDIR)/Basic3DParams.cpp -o $@

$(BUILDDIR)/Parameters.o: $(SRCDIR)/Parameters.h $(SRCDIR)/Parameters.cpp
	$(CC) $(CFLAGS) $(SRCDIR)/Parameters.cpp -o $@

clean:
	\rm $(BUILDDIR)/*o


