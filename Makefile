include config.mk

all: lib/libcritter.a

lib/libcritter.a: obj/critter.o
	ar -crs lib/libcritter.a obj/critter.o

obj/critter.o: critter.h critter.cxx
	$(CXX) critter.cxx -c -o obj/critter.o $(CXXFLAGS)

clean:
	rm obj/critter.o lib/libcritter.a
