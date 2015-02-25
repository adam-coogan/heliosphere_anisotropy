CC=clang++ -std=c++11
CFLAGS=-c
LFLAGS=
OBJS=PPTrajectory.o SolarProp.o PPTrajectoryParams.o

SolarProp: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o SolarProp

SolarProp.o: SolarProp.cpp PPTrajectory.h
	$(CC) $(CFLAGS) SolarProp.cpp

PPTrajectory.o: PPTrajectory.h PPTrajectory.cpp PPTrajectoryParams.h
	$(CC) $(CFLAGS) PPTrajectory.cpp

PPTrajectoryParams.o: PPTrajectoryParams.h PPTrajectoryParams.cpp
	$(CC) $(CFLAGS) PPTrajectoryParams.cpp

PPTrajectory.h: PPPoint.h SphericalTensor.h SphericalVector.h PPTrajectoryParams.h

clean:
	\rm *o


