# Makefile for lake problem
CC = g++
CFLAGS = -O3 -Wall -Wno-unused-local-typedefs -ggdb
INCL = -I boost_1_56_0 

SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
EXE = portfolio.exe

all: $(SOURCES) $(EXE)
	rm $(OBJECTS)

.cpp.o:
	$(CC) -c $(CFLAGS) $^ -o $@ $(INCL)
	
$(EXE): $(OBJECTS)
	$(CC) $(OBJECTS) $(CFLAGS) -o $@ $(INCL)

clean:
	rm -f $(OBJECTS) $(EXE)
