CC=g++
CFLAGS=-std=c++11 -I /home/k-class/rdcunha2/eigen-eigen-67e894c6cd8f
SOURCES=2.cpp
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=2

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $< -o $@
