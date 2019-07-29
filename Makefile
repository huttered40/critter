include config.mk

all: lib/libcritter.a

test: bin/test_bcast

lib/libcritter.a: obj/critter.o
	ar -crs lib/libcritter.a obj/critter.o

obj/critter.o: src/critter.h src/critter.cxx
	$(CXX) src/critter.cxx -c -o obj/critter.o $(CXXFLAGS)

bin/test_bcast: lib/libcritter.a test/test_bcast.cxx
	$(CXX) test/test_bcast.cxx -o bin/test_bcast $(CXXFLAGS) -L./lib -lcritter $(LDFLAGS)

clean:
	rm -f obj/critter.o lib/libcritter.a bin/test_bcast
