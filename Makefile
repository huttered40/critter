include config.mk

all: lib/libcritter.a

test: bin/test_bcast bin/test_bcast bin/test_reduce bin/test_allreduce bin/test_allgather bin/test_sendrecv_replace

lib/libcritter.a: obj/critter.o
	ar -crs lib/libcritter.a obj/critter.o

obj/critter.o: src/critter.h src/critter.cxx
	$(CXX) src/critter.cxx -c -o obj/critter.o $(CXXFLAGS)

bin/test_bcast: lib/libcritter.a bmark/test_bcast.cxx
	$(CXX) bmark/test_bcast.cxx -o bin/test_bcast $(CXXFLAGS) -L./lib -lcritter $(LDFLAGS)
bin/test_reduce: lib/libcritter.a bmark/test_reduce.cxx
	$(CXX) bmark/test_reduce.cxx -o bin/test_reduce $(CXXFLAGS) -L./lib -lcritter $(LDFLAGS)
bin/test_allreduce: lib/libcritter.a bmark/test_allreduce.cxx
	$(CXX) bmark/test_allreduce.cxx -o bin/test_allreduce $(CXXFLAGS) -L./lib -lcritter $(LDFLAGS)
bin/test_allgather: lib/libcritter.a bmark/test_allgather.cxx
	$(CXX) bmark/test_allgather.cxx -o bin/test_allgather $(CXXFLAGS) -L./lib -lcritter $(LDFLAGS)
bin/test_sendrecv_replace: lib/libcritter.a bmark/test_sendrecv_replace.cxx
	$(CXX) bmark/test_sendrecv_replace.cxx -o bin/test_sendrecv_replace $(CXXFLAGS) -L./lib -lcritter $(LDFLAGS)

clean:
	rm -f obj/critter.o lib/libcritter.a bin/test_bcast bin/test_reduce bin/test_allreduce bin/test_allgather bin/test_sendrecv_replace
