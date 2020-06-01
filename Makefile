include config.mk

all: lib/libcritter.a

test: lib/libcritter.a bin/test_bcast bin/test_gather bin/test_reduce bin/test_reducescatter bin/test_scatter bin/test_allreduce bin/test_allgather bin/test_alltoall bin/test_sendrecv_replace bin/test_send_recv bin/test_sendrecv bin/test_waitany

lib/libcritter.a: obj/critter.o
	ar -crs lib/libcritter.a obj/critter.o

lib/libcritter.so: obj/critter.o
	gcc -shared -o lib/libcritter.so obj/critter.o

obj/critter.o: src/critter.cxx
	$(CXX) src/critter.cxx -c -o obj/critter.o $(CXXFLAGS)

bin/test_bcast: lib/libcritter.a bmark/test_bcast.cxx
	$(CXX) bmark/test_bcast.cxx -o bin/test_bcast $(CXXFLAGS) -L$(HOME)/critter/lib -lcritter $(LDFLAGS)
bin/test_gather: lib/libcritter.a bmark/test_gather.cxx
	$(CXX) bmark/test_gather.cxx -o bin/test_gather $(CXXFLAGS) -L$(HOME)/critter/lib -lcritter $(LDFLAGS)
bin/test_reduce: lib/libcritter.a bmark/test_reduce.cxx
	$(CXX) bmark/test_reduce.cxx -o bin/test_reduce $(CXXFLAGS) -L$(HOME)/critter/lib -lcritter $(LDFLAGS)
bin/test_reducescatter: lib/libcritter.a bmark/test_reducescatter.cxx
	$(CXX) bmark/test_reducescatter.cxx -o bin/test_reducescatter $(CXXFLAGS) -L$(HOME)/critter/lib -lcritter $(LDFLAGS)
bin/test_scatter: lib/libcritter.a bmark/test_scatter.cxx
	$(CXX) bmark/test_scatter.cxx -o bin/test_scatter $(CXXFLAGS) -L$(HOME)/critter/lib -lcritter $(LDFLAGS)
bin/test_allgather: lib/libcritter.a bmark/test_allgather.cxx
	$(CXX) bmark/test_allgather.cxx -o bin/test_allgather $(CXXFLAGS) -L$(HOME)/critter/lib -lcritter $(LDFLAGS)
bin/test_allreduce: lib/libcritter.a bmark/test_allreduce.cxx
	$(CXX) bmark/test_allreduce.cxx -o bin/test_allreduce $(CXXFLAGS) -L$(HOME)/critter/lib -lcritter $(LDFLAGS)
bin/test_alltoall: lib/libcritter.a bmark/test_alltoall.cxx
	$(CXX) bmark/test_alltoall.cxx -o bin/test_alltoall $(CXXFLAGS) -L$(HOME)/critter/lib -lcritter $(LDFLAGS)
bin/test_sendrecv_replace: lib/libcritter.a bmark/test_sendrecv_replace.cxx
	$(CXX) bmark/test_sendrecv_replace.cxx -o bin/test_sendrecv_replace $(CXXFLAGS) -L$(HOME)/critter/lib -lcritter $(LDFLAGS)
bin/test_send_recv: lib/libcritter.a bmark/test_send_recv.cxx
	$(CXX) bmark/test_send_recv.cxx -o bin/test_send_recv $(CXXFLAGS) -L$(HOME)/critter/lib -lcritter $(LDFLAGS)
bin/test_sendrecv: lib/libcritter.a bmark/test_sendrecv.cxx
	$(CXX) bmark/test_sendrecv.cxx -o bin/test_sendrecv $(CXXFLAGS) -L$(HOME)/critter/lib -lcritter $(LDFLAGS)
bin/test_waitany: lib/libcritter.a bmark/test_waitany.cxx
	$(CXX) bmark/test_waitany.cxx -o bin/test_waitany $(CXXFLAGS) $(HOME)/critter/lib/libcritter.a $(LDFLAGS)

clean:
	rm -f obj/critter.o lib/libcritter.a lib/libcritter.so bin/test_bcast bin/test_gather bin/test_reduce bin/test_reducescatter bin/test_scatter bin/test_allreduce bin/test_allgather bin/test_alltoall bin/test_sendrecv_replace bin/test_send_recv bin/test_sendrecv bin/test_waitany
