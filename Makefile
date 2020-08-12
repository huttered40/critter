include config.mk

all: lib/libcritter.a

test: lib/libcritter.a

lib/libcritter.a:\
		obj/util_util.o\
		obj/intercept_comp.o\
		obj/intercept_comm.o\
		obj/intercept_symbol.o\
		obj/decomposition_util_util.o\
		obj/decomposition_record_record.o\
		obj/decomposition_container_comm_tracker.o\
		obj/decomposition_container_symbol_tracker.o\
		obj/decomposition_volumetric_volumetric.o\
		obj/decomposition_path_path.o\
		obj/discretization_parameterization_parameterization.o\
		obj/discretization_util_util.o\
		obj/discretization_record_record.o\
		obj/discretization_container_comm_tracker.o\
		obj/discretization_container_symbol_tracker.o\
		obj/discretization_volumetric_volumetric.o\
		obj/discretization_path_path.o\
		obj/dispatch_dispatch.o
	ar -crs lib/libcritter.a obj/util_util.o obj/intercept_comp.o obj/intercept_comm.o obj/intercept_symbol.o\
					obj/decomposition_util_util.o obj/decomposition_record_record.o\
					obj/decomposition_container_comm_tracker.o obj/decomposition_container_symbol_tracker.o\
					obj/decomposition_volumetric_volumetric.o  obj/decomposition_path_path.o\
					obj/discretization_util_util.o obj/discretization_parameterization_parameterization.o\
					obj/discretization_record_record.o obj/discretization_container_comm_tracker.o\
					obj/discretization_container_symbol_tracker.o obj/discretization_volumetric_volumetric.o\
					obj/discretization_path_path.o obj/dispatch_dispatch.o

lib/libcritter.so: obj/critter.o
	gcc -shared -o lib/libcritter.so obj util.o obj/critter.o

obj/util_util.o: src/util/util.cxx
	$(CXX) src/util/util.cxx -c -o obj/util_util.o $(CXXFLAGS)

obj/intercept_comp.o: src/intercept/comp.cxx
	$(CXX) src/intercept/comp.cxx -c -o obj/intercept_comp.o $(CXXFLAGS)

obj/intercept_comm.o: src/intercept/comm.cxx
	$(CXX) src/intercept/comm.cxx -c -o obj/intercept_comm.o $(CXXFLAGS)

obj/intercept_symbol.o: src/intercept/symbol.cxx
	$(CXX) src/intercept/symbol.cxx -c -o obj/intercept_symbol.o $(CXXFLAGS)

obj/decomposition_util_util.o: src/decomposition/util/util.cxx
	$(CXX) src/decomposition/util/util.cxx -c -o obj/decomposition_util_util.o $(CXXFLAGS)

obj/decomposition_record_record.o: src/decomposition/record/record.cxx
	$(CXX) src/decomposition/record/record.cxx -c -o obj/decomposition_record_record.o $(CXXFLAGS)

obj/decomposition_container_comm_tracker.o: src/decomposition/container/comm_tracker.cxx
	$(CXX) src/decomposition/container/comm_tracker.cxx -c -o obj/decomposition_container_comm_tracker.o $(CXXFLAGS)

obj/decomposition_container_symbol_tracker.o: src/decomposition/container/symbol_tracker.cxx
	$(CXX) src/decomposition/container/symbol_tracker.cxx -c -o obj/decomposition_container_symbol_tracker.o $(CXXFLAGS)

obj/decomposition_volumetric_volumetric.o: src/decomposition/volumetric/volumetric.cxx
	$(CXX) src/decomposition/volumetric/volumetric.cxx -c -o obj/decomposition_volumetric_volumetric.o $(CXXFLAGS)

obj/decomposition_path_path.o: src/decomposition/path/path.cxx
	$(CXX) src/decomposition/path/path.cxx -c -o obj/decomposition_path_path.o $(CXXFLAGS)

obj/discretization_util_util.o: src/discretization/util/util.cxx
	$(CXX) src/discretization/util/util.cxx -c -o obj/discretization_util_util.o $(CXXFLAGS)

obj/discretization_parameterization_parameterization.o: src/discretization/parameterization/parameterization.cxx
	$(CXX) src/discretization/parameterization/parameterization.cxx -c -o obj/discretization_parameterization_parameterization.o $(CXXFLAGS)

obj/discretization_record_record.o: src/discretization/record/record.cxx
	$(CXX) src/discretization/record/record.cxx -c -o obj/discretization_record_record.o $(CXXFLAGS)

obj/discretization_container_comm_tracker.o: src/discretization/container/comm_tracker.cxx
	$(CXX) src/discretization/container/comm_tracker.cxx -c -o obj/discretization_container_comm_tracker.o $(CXXFLAGS)

obj/discretization_container_symbol_tracker.o: src/discretization/container/symbol_tracker.cxx
	$(CXX) src/discretization/container/symbol_tracker.cxx -c -o obj/discretization_container_symbol_tracker.o $(CXXFLAGS)

obj/discretization_volumetric_volumetric.o: src/discretization/volumetric/volumetric.cxx
	$(CXX) src/discretization/volumetric/volumetric.cxx -c -o obj/discretization_volumetric_volumetric.o $(CXXFLAGS)

obj/discretization_path_path.o: src/discretization/path/path.cxx
	$(CXX) src/discretization/path/path.cxx -c -o obj/discretization_path_path.o $(CXXFLAGS)

obj/dispatch_dispatch.o: src/dispatch/dispatch.cxx
	$(CXX) src/dispatch/dispatch.cxx -c -o obj/dispatch_dispatch.o $(CXXFLAGS)

obj/replay_path_path.o: src/replay/path/path.cxx
	$(CXX) src/replay/path/path.cxx -c -o obj/replay_path_path.o $(CXXFLAGS)

clean:
	rm -f obj/*.o lib/libcritter.a lib/libcritter.so
