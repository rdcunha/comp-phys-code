CC=g++
CFLAGS=-std=c++11 -I /home/k-class/rdcunha2/eigen-eigen-67e894c6cd8f
SOURCES=1.cpp
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=1

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(CFLAGS) $< -o $@
