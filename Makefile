include config.mk

all: lib/libcritter.a

test: lib/libcritter.a

lib/libcritter.a:\
		obj/util.o\
		obj/intercept_comm.o\
		obj/intercept_symbol.o\
		obj/record_record.o\
		obj/intercept_symbol.o\
		obj/container_comm_tracker.o\
		obj/container_symbol_tracker.o\
		obj/mechanism_volumetric_volumetric.o\
		obj/mechanism_per-process_per-process.o\
		obj/mechanism_path_dispatch.o\
		obj/mechanism_path_forward_pass.o\
		obj/intercept_symbol.o
	ar -crs lib/libcritter.a obj/util.o obj/intercept_comm.o obj/intercept_symbol.o obj/intercept_symbol.o obj/record_record.o\
					obj/container_comm_tracker.o obj/container_symbol_tracker.o\
					obj/mechanism_volumetric_volumetric.o obj/mechanism_per-process_per-process.o\
					obj/mechanism_path_dispatch.o obj/mechanism_path_forward_pass.o obj/intercept_symbol.o

lib/libcritter.so: obj/critter.o
	gcc -shared -o lib/libcritter.so obj util.o obj/critter.o

obj/util.o: src/util.cxx
	$(CXX) src/util.cxx -c -o obj/util.o $(CXXFLAGS)

obj/intercept_comm.o: src/intercept/comm.cxx
	$(CXX) src/intercept/comm.cxx -c -o obj/intercept_comm.o $(CXXFLAGS)

obj/intercept_symbol.o: src/intercept/symbol.cxx
	$(CXX) src/intercept/symbol.cxx -c -o obj/intercept_symbol.o $(CXXFLAGS)

obj/record_record.o: src/record/record.cxx
	$(CXX) src/record/record.cxx -c -o obj/record_record.o $(CXXFLAGS)

obj/container_comm_tracker.o: src/container/comm_tracker.cxx
	$(CXX) src/container/comm_tracker.cxx -c -o obj/container_comm_tracker.o $(CXXFLAGS)

obj/container_symbol_tracker.o: src/container/symbol_tracker.cxx
	$(CXX) src/container/symbol_tracker.cxx -c -o obj/container_symbol_tracker.o $(CXXFLAGS)

obj/mechanism_volumetric_volumetric.o: src/mechanism/volumetric/volumetric.cxx
	$(CXX) src/mechanism/volumetric/volumetric.cxx -c -o obj/mechanism_volumetric_volumetric.o $(CXXFLAGS)

obj/mechanism_per-process_per-process.o: src/mechanism/per-process/per-process.cxx
	$(CXX) src/mechanism/per-process/per-process.cxx -c -o obj/mechanism_per-process_per-process.o $(CXXFLAGS)

obj/mechanism_path_dispatch.o: src/mechanism/path/dispatch.cxx
	$(CXX) src/mechanism/path/dispatch.cxx -c -o obj/mechanism_path_dispatch.o $(CXXFLAGS)

obj/mechanism_path_forward_pass.o: src/mechanism/path/forward_pass.cxx
	$(CXX) src/mechanism/path/forward_pass.cxx -c -o obj/mechanism_path_forward_pass.o $(CXXFLAGS)

clean:
	rm -f obj/*.o lib/libcritter.a lib/libcritter.so
