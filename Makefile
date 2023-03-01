CXX=g++
INCLUDES=
CFLAGS=

all: addNH mmannot

mmannot: mmannot.cpp
	$(CXX) mmannot.cpp $(INCLUDES) $(CFLAGS) -DMMSTANDALONE -Wpedantic -std=c++11 -pthread -lz -O3 -o mmannot

addNH: addNH.cpp
	$(CXX) addNH.cpp $(INCLUDES) $(CFLAGS) -DMMSTANDALONE -Wpedantic -std=c++11 -pthread -lz -O3 -o addNH

clean:
	\rm -f mmannot addNH
