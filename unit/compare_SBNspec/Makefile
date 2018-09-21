CC=g++

CXXFLAGS=-std=c++11

INCFLAGS=`root-config --cflags`
LDFLAGS=`root-config --libs`

SOURCES=compare_cv.cxx
OBJECTS=compare_cv.o
EXECUTABLE=compare_cv

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $^ -o $@ $(LDFLAGS)

$(OBJECTS): $(SOURCES)
	$(CC) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@

clean:
	rm *.o compare_cv
