CC=clang++ -std=c++11 -O2
CFLAGS=-c
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


