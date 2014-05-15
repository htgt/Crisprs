CC=g++
CXXFLAGS=-c -std=c++0x -O3 -W -Wall

all: find_off_targets

find_off_targets: utils.o crisprutil.o find_off_targets.o
	$(CC) utils.o crisprutil.o find_off_targets.o -o find_off_targets

utils.o: find_off_targets.cpp
	$(CC) $(CXXFLAGS) utils.cpp

crisprutil.o: crisprutil.cpp
	$(CC) $(CXXFLAGS) crisprutil.cpp

find_off_targets.o: find_off_targets.cpp
	$(CC) $(CXXFLAGS) find_off_targets.cpp

clean:
	rm -f *.o find_off_targets
