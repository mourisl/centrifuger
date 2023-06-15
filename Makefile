CXX = g++
CXXFLAGS= -Wall -g -O3 -msse4.2 #-pg -g #-Wall #-O3
LINKPATH=
LINKFLAGS = -lpthread -lz 
DEBUG=
OBJECTS = main.o 

#asan=1
ifneq ($(asan),)
	CXXFLAGS+=-fsanitize=address -g
	LDFLAGS+=-fsanitize=address -ldl -g
endif

all: centrifuger-class centrifuger-build centrifuger-inspect

centrifuger-build: CentrifugerBuild.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)

centrifuger-class: CentrifugerClass.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)

centrifuger-inspect: CentrifugerInspect.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)

CentrifugerBuild.o: CentrifugerBuild.cpp Builder.hpp ReadFiles.hpp Taxonomy.hpp defs.h compactds/*.hpp 
CentrifugerClass.o: CentrifugerClass.cpp Classifier.hpp ReadFiles.hpp Taxonomy.hpp defs.h ResultWriter.hpp compactds/*.hpp 
CentrifugerInspect: CentrifugerInspect.cpp Taxonomy.hpp 
clean:
	rm -f *.o centrifuger-build centrifuger-class centrifuger-inspect
