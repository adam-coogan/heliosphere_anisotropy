CC=clang++ -std=c++11 -O2
CFLAGS=-c
LFLAGS=
OBJS=$(BUILDDIR)/NewSolarProp.o
BINDIR=bin
BUILDDIR=build
SRCDIR=src

$(BINDIR)/NewSolarProp: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o $(BINDIR)/NewSolarProp

$(BUILDDIR)/NewSolarProp.o: $(SRCDIR)/NewSolarProp.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	\rm $(BUILDDIR)/*o


