CC=clang++ -std=c++11
CFLAGS=-c
LFLAGS=
OBJS=PPTrajectory.o SolarProp.o ParamExceptions.o

SolarProp: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o SolarProp

SolarProp.o: SolarProp.cpp PPTrajectory.h RunConfig.h
	$(CC) $(CFLAGS) SolarProp.cpp

PPTrajectory.o: PPTrajectory.h PPTrajectory.cpp HelioParams.h
	$(CC) $(CFLAGS) PPTrajectory.cpp

ParamExceptions.o: ParamExceptions.h ParamExceptions.cpp
	$(CC) $(CFLAGS) ParamExceptions.cpp

PPTrajectory.h: PPPoint.h SphericalTensor.h SphericalVector.h HelioParams.h

HelioParams.h: Params.h

RunConfig.h: Params.h

clean:
	\rm *o


