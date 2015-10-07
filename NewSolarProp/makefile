ifeq ($(shell hostname),hyades.ucsc.edu)
	CC=icpc -std=c++11 -I/pfs/sw/serial/gcc/boost-1.57.0/include # Broken
else
	CC=clang++ -std=c++11
endif
CFLAGS=-c -O2
LFLAGS=
OBJS=$(BUILDDIR)/StraussUnits.o
BINDIR=bin
BUILDDIR=build
SRCDIR=src

$(BINDIR)/StraussUnits: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o $(BINDIR)/StraussUnits

$(BUILDDIR)/StraussUnits.o: $(SRCDIR)/StraussUnits.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	\rm $(BUILDDIR)/*o


