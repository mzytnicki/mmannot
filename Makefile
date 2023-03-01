CXX = g++
CFLAGS = -Wpedantic -std=c++11 -pthread -lz
LDFLAGS = 

ifdef STATIC
	LDFLAGS += -static
endif

ifdef DEBUG
	CFLAGS += -O0 -g -ggdb
else ifdef PROF
		CFLAGS += -O3 -g -pg
else
		CFLAGS += -O3
endif

all: addNH mmannot

mmannot: mmannot.cpp
	$(CXX) mmannot.cpp $(LDFLAGS) $(CFLAGS) -o mmannot

addNH: addNH.cpp
	$(CXX) addNH.cpp $(LDFLAGS) $(CFLAGS) -o addNH

clean:
	\rm -f mmannot addNH
