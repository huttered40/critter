include config.mk

all: lib/libcritter.a

test: lib/libcritter.a

lib/libcritter.a: obj/critter.o
	ar -crs lib/libcritter.a obj/critter.o

lib/libcritter.so: obj/critter.o
	gcc -shared -o lib/libcritter.so obj/critter.o

obj/critter.o: src/critter.cxx
	$(CXX) src/critter.cxx -c -o obj/critter.o $(CXXFLAGS)

clean:
	rm -f obj/critter.o lib/libcritter.a lib/libcritter.so
