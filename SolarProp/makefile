CC=clang++ -std=c++11
CFLAGS=-c
LFLAGS=
OBJS=$(BUILDDIR)/PPTrajectory.o $(BUILDDIR)/SolarProp.o $(BUILDDIR)/ParamExceptions.o $(BUILDDIR)/SolarPropMC.o
BINDIR=bin
BUILDDIR=build
SRCDIR=src

$(BINDIR)/SolarProp: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o $(BINDIR)/SolarProp

$(BUILDDIR)/SolarProp.o: $(SRCDIR)/SolarProp.cpp $(SRCDIR)/PPTrajectory.h $(SRCDIR)/RunConfig.h
	$(CC) $(CFLAGS) $(SRCDIR)/SolarProp.cpp -o $@

$(BUILDDIR)/SolarPropMC.o: $(SRCDIR)/SolarProp.cpp $(SRCDIR)/PPTrajectory.h $(SRCDIR)/RunConfig.h
	$(CC) $(CFLAGS) $(SRCDIR)/SolarPropMC.cpp -o $@

$(BUILDDIR)/PPTrajectory.o: $(SRCDIR)/PPTrajectory.h $(SRCDIR)/PPTrajectory.cpp $(SRCDIR)/HelioParams.h
	$(CC) $(CFLAGS) $(SRCDIR)/PPTrajectory.cpp -o $@

$(BUILDDIR)/ParamExceptions.o: $(SRCDIR)/ParamExceptions.h $(SRCDIR)/ParamExceptions.cpp
	$(CC) $(CFLAGS) $(SRCDIR)/ParamExceptions.cpp -o $@

$(SRCDIR)/PPTrajectory.h: $(SRCDIR)/PPPoint.h $(SRCDIR)/SphericalTensor.h $(SRCDIR)/SphericalVector.h $(SRCDIR)/HelioParams.h

$(SRCDIR)/HelioParams.h: $(SRCDIR)/Params.h

$(SRCDIR)/RunConfig.h: $(SRCDIR)/Params.h

clean:
	\rm $(BUILDDIR)/*o


