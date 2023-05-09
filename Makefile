include config.mk

all: lib/libcritter.a lib/libcritter.so

static: lib/libcritter.a

shared: lib/libcritter.so

test: lib/libcritter.a

lib/libcritter.a:\
		obj/interface.o\
		obj/util.o\
		obj/intercept_comp.o\
		obj/intercept_comm.o\
		obj/record.o\
		obj/container_comm_tracker.o\
		obj/container_kernel_tracker.o\
		obj/profiler.o
	ar -crs lib/libcritter.a\
					obj/interface.o\
					obj/util.o\
					obj/intercept_comp.o\
					obj/intercept_comm.o\
					obj/record.o\
					obj/container_comm_tracker.o\
					obj/container_kernel_tracker.o\
					obj/profiler.o

lib/libcritter.so:\
		obj/interface.o\
		obj/util.o\
		obj/intercept_comp.o\
		obj/intercept_comm.o\
		obj/record.o\
		obj/container_comm_tracker.o\
		obj/container_kernel_tracker.o\
		obj/profiler.o
	$(CXX) -shared -o lib/libcritter.so\
					obj/interface.o\
					obj/util.o\
					obj/intercept_comp.o\
					obj/intercept_comm.o\
					obj/record.o\
					obj/container_comm_tracker.o\
					obj/container_kernel_tracker.o\
					obj/profiler.o

obj/interface.o: src/interface.cxx
	$(CXX) src/interface.cxx -c -o obj/interface.o $(CXXFLAGS)

obj/util.o: src/util.cxx
	$(CXX) src/util.cxx -c -o obj/util.o $(CXXFLAGS)

obj/intercept_comp.o: src/intercept/comp.cxx
	$(CXX) src/intercept/comp.cxx -c -o obj/intercept_comp.o $(CXXFLAGS)

obj/intercept_comm.o: src/intercept/comm.cxx
	$(CXX) src/intercept/comm.cxx -c -o obj/intercept_comm.o $(CXXFLAGS)

obj/record.o: src/record.cxx
	$(CXX) src/record.cxx -c -o obj/record.o $(CXXFLAGS)

obj/container_comm_tracker.o: src/container/comm_tracker.cxx
	$(CXX) src/container/comm_tracker.cxx -c -o obj/container_comm_tracker.o $(CXXFLAGS)

obj/container_kernel_tracker.o: src/container/kernel_tracker.cxx
	$(CXX) src/container/kernel_tracker.cxx -c -o obj/container_kernel_tracker.o $(CXXFLAGS)

obj/profiler.o: src/profiler.cxx
	$(CXX) src/profiler.cxx -c -o obj/profiler.o $(CXXFLAGS)

clean:
	rm -f obj/*.o lib/libcritter.a lib/libcritter.so
