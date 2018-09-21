CC=g++

CXXFLAGS=-std=c++11

INCFLAGS=`root-config --cflags`
LDFLAGS=`root-config --libs`

SOURCES=compare_covariance.cxx
OBJECTS=compare_covariance.o
EXECUTABLE=compare_covariance

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $^ -o $@ $(LDFLAGS)

$(OBJECTS): $(SOURCES)
	$(CC) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@

clean:
	rm *.o compare_covariance
