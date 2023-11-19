CXX = g++
CXXFLAGS= -Wall -g -msse4.2 -O3 #-pg -g #-Wall #-O3
LINKPATH= 
LINKFLAGS = -lpthread -lz 
DEBUG=
OBJECTS =  #BaseReads.o Alignment.o 
HEADERS = *.hpp

#asan=1
ifneq ($(asan),)
	CXXFLAGS+=-fsanitize=address -g
	LDFLAGS+=-fsanitize=address -ldl -g
endif

#all: bitvector-benchmark #test #bitvector-benchmark
all: test #rbbwt #bitvector-benchmark

test: test.o $(OBJECTS)
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $< $(OBJECTS) $(LINKFLAGS)

bitvector-benchmark: bitvector_benchmark.cpp $(HEADERS)
	$(CXX) -o $@ $(LINKPATH) $< $(LINKFLAGS) -std=c++11 -O3 -ffast-math -funroll-loops -msse4.2 -march=native -DHAVE_CXA_DEMANGLE
	#$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $< $(OBJECTS) -std=c++11 $(LINKFLAGS)

rbbwt: rbbwt.cpp $(HEADERS)
	$(CXX) -o $@ $(LINKPATH) $< $(LINKFLAGS) -g -std=c++11 -Ofast -march=native -fstrict-aliasing 

test.o: test.cpp $(HEADERS) 

clean:
	rm -f *.o *.gch test
