ifeq ($(shell hostname),hyades.ucsc.edu)
	CC=icpc -std=c++11 -I/pfs/sw/serial/gcc/boost-1.57.0/include
	LFLAGS=/pfs/sw/serial/gcc/boost-1.57.0/lib/libboost_program_options.so.1.57.0
else
	CC=clang++ -std=c++11 -I /usr/local/opt/boost/include
	LFLAGS=-lboost_program_options
endif

CFLAGS=-c -O2
OBJS=$(BUILDDIR)/Parameters.o $(BUILDDIR)/Basic3DParams.o $(BUILDDIR)/Point.o $(BUILDDIR)/HeliosCLParser.o \
	 $(BUILDDIR)/Wavy3DParams.o
BINDIR=bin
BUILDDIR=build
SRCDIR=src

all: Helios

# Messy program using Strauss' units

StraussUnitsModified: $(BUILDDIR)/StraussUnitsModified.o
	$(CC) $(LFLAGS) $(BUILDDIR)/StraussUnitsModified.o -o $(BINDIR)/StraussUnitsModified

$(BUILDDIR)/StraussUnitsModified.o: $(SRCDIR)/StraussUnitsModified.cpp
	$(CC) $(CFLAGS) $< -o $@

# Test program

HeliosTest: $(OBJS) $(BUILDDIR)/HeliosTest.o $(SRCDIR)/Basic3D.h $(SRCDIR)/Wavy3D.h
	$(CC) $(LFLAGS) $(OBJS) $(BUILDDIR)/HeliosTest.o -o $(BINDIR)/HeliosTest

$(BUILDDIR)/HeliosTest.o: $(SRCDIR)/HeliosTest.cpp $(SRCDIR)/Basic3D.h
	$(CC) $(CFLAGS) $< -o $@

# Main program

Helios: $(OBJS) $(BUILDDIR)/Helios.o
	$(CC) $(LFLAGS) $(OBJS) $(BUILDDIR)/Helios.o -o $(BINDIR)/Helios

$(BUILDDIR)/Helios.o: $(SRCDIR)/Helios.cpp $(SRCDIR)/Basic3D.h $(SRCDIR)/Simulation.h $(SRCDIR)/Wavy3D.h \
	$(SRCDIR)/Wavy3D.h
	$(CC) $(CFLAGS) $< -o $@

# Objects

$(BUILDDIR)/HeliosCLParser.o: $(SRCDIR)/HeliosCLParser.h $(SRCDIR)/HeliosCLParser.cpp
	$(CC) $(CFLAGS) $(SRCDIR)/HeliosCLParser.cpp -o $@

$(BUILDDIR)/Wavy3DParams.o: $(SRCDIR)/Basic3DParams.h $(SRCDIR)/Wavy3DParams.h $(SRCDIR)/Wavy3DParams.cpp
	$(CC) $(CFLAGS) $(SRCDIR)/Wavy3DParams.cpp -o $@

$(BUILDDIR)/Basic3DParams.o: $(SRCDIR)/Parameters.h $(SRCDIR)/Basic3DParams.h $(SRCDIR)/Basic3DParams.cpp
	$(CC) $(CFLAGS) $(SRCDIR)/Basic3DParams.cpp -o $@

$(BUILDDIR)/Parameters.o: $(SRCDIR)/Parameters.h $(SRCDIR)/Parameters.cpp
	$(CC) $(CFLAGS) $(SRCDIR)/Parameters.cpp -o $@

$(BUILDDIR)/Point.o: $(SRCDIR)/Point.h $(SRCDIR)/Point.cpp
	$(CC) $(CFLAGS) $(SRCDIR)/Point.cpp -o $@

clean:
	\rm $(BUILDDIR)/*o


