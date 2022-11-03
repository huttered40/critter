include config.mk

all: lib/libcritter.a

test: lib/libcritter.a

lib/libcritter.a:\
		obj/util_util.o\
		obj/intercept_comp.o\
		obj/intercept_comm.o\
		obj/intercept_symbol.o\
		obj/profile_util_util.o\
		obj/profile_record_record.o\
		obj/profile_container_comm_tracker.o\
		obj/profile_container_symbol_tracker.o\
		obj/profile_volumetric_volumetric.o\
		obj/profile_path_path.o\
		obj/accelerate_util_util.o\
		obj/accelerate_record_record.o\
		obj/accelerate_container_comm_tracker.o\
		obj/accelerate_container_symbol_tracker.o\
		obj/accelerate_volumetric_volumetric.o\
		obj/accelerate_path_path.o\
		obj/skeletonize_util_util.o\
		obj/skeletonize_record_record.o\
		obj/skeletonize_container_comm_tracker.o\
		obj/skeletonize_container_symbol_tracker.o\
		obj/skeletonize_volumetric_volumetric.o\
		obj/skeletonize_path_path.o\
		obj/dispatch_dispatch.o
	ar -crs lib/libcritter.a obj/util_util.o\
					obj/intercept_comp.o\
					obj/intercept_comm.o\
					obj/intercept_symbol.o\
					obj/profile_util_util.o\
					obj/profile_record_record.o\
					obj/profile_container_comm_tracker.o\
					obj/profile_container_symbol_tracker.o\
					obj/profile_volumetric_volumetric.o\
					obj/profile_path_path.o\
					obj/accelerate_util_util.o\
					obj/accelerate_record_record.o\
					obj/accelerate_container_comm_tracker.o\
					obj/accelerate_container_symbol_tracker.o\
					obj/accelerate_volumetric_volumetric.o\
					obj/accelerate_path_path.o\
					obj/skeletonize_util_util.o\
					obj/skeletonize_record_record.o\
					obj/skeletonize_container_comm_tracker.o\
					obj/skeletonize_container_symbol_tracker.o\
					obj/skeletonize_volumetric_volumetric.o\
					obj/skeletonize_path_path.o\
					obj/dispatch_dispatch.o

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

obj/profile_util_util.o: src/profile/util/util.cxx
	$(CXX) src/profile/util/util.cxx -c -o obj/profile_util_util.o $(CXXFLAGS)

obj/profile_record_record.o: src/profile/record/record.cxx
	$(CXX) src/profile/record/record.cxx -c -o obj/profile_record_record.o $(CXXFLAGS)

obj/profile_container_comm_tracker.o: src/profile/container/comm_tracker.cxx
	$(CXX) src/profile/container/comm_tracker.cxx -c -o obj/profile_container_comm_tracker.o $(CXXFLAGS)

obj/profile_container_symbol_tracker.o: src/profile/container/symbol_tracker.cxx
	$(CXX) src/profile/container/symbol_tracker.cxx -c -o obj/profile_container_symbol_tracker.o $(CXXFLAGS)

obj/profile_volumetric_volumetric.o: src/profile/volumetric/volumetric.cxx
	$(CXX) src/profile/volumetric/volumetric.cxx -c -o obj/profile_volumetric_volumetric.o $(CXXFLAGS)

obj/profile_path_path.o: src/profile/path/path.cxx
	$(CXX) src/profile/path/path.cxx -c -o obj/profile_path_path.o $(CXXFLAGS)

obj/accelerate_util_util.o: src/accelerate/util/util.cxx
	$(CXX) src/accelerate/util/util.cxx -c -o obj/accelerate_util_util.o $(CXXFLAGS)

obj/accelerate_record_record.o: src/accelerate/record/record.cxx
	$(CXX) src/accelerate/record/record.cxx -c -o obj/accelerate_record_record.o $(CXXFLAGS)

obj/accelerate_container_comm_tracker.o: src/accelerate/container/comm_tracker.cxx
	$(CXX) src/accelerate/container/comm_tracker.cxx -c -o obj/accelerate_container_comm_tracker.o $(CXXFLAGS)

obj/accelerate_container_symbol_tracker.o: src/accelerate/container/symbol_tracker.cxx
	$(CXX) src/accelerate/container/symbol_tracker.cxx -c -o obj/accelerate_container_symbol_tracker.o $(CXXFLAGS)

obj/accelerate_volumetric_volumetric.o: src/accelerate/volumetric/volumetric.cxx
	$(CXX) src/accelerate/volumetric/volumetric.cxx -c -o obj/accelerate_volumetric_volumetric.o $(CXXFLAGS)

obj/accelerate_path_path.o: src/accelerate/path/path.cxx
	$(CXX) src/accelerate/path/path.cxx -c -o obj/accelerate_path_path.o $(CXXFLAGS)

obj/skeletonize_util_util.o: src/skeletonize/util/util.cxx
	$(CXX) src/skeletonize/util/util.cxx -c -o obj/skeletonize_util_util.o $(CXXFLAGS)

obj/skeletonize_record_record.o: src/skeletonize/record/record.cxx
	$(CXX) src/skeletonize/record/record.cxx -c -o obj/skeletonize_record_record.o $(CXXFLAGS)

obj/skeletonize_container_comm_tracker.o: src/skeletonize/container/comm_tracker.cxx
	$(CXX) src/skeletonize/container/comm_tracker.cxx -c -o obj/skeletonize_container_comm_tracker.o $(CXXFLAGS)

obj/skeletonize_container_symbol_tracker.o: src/skeletonize/container/symbol_tracker.cxx
	$(CXX) src/skeletonize/container/symbol_tracker.cxx -c -o obj/skeletonize_container_symbol_tracker.o $(CXXFLAGS)

obj/skeletonize_volumetric_volumetric.o: src/skeletonize/volumetric/volumetric.cxx
	$(CXX) src/skeletonize/volumetric/volumetric.cxx -c -o obj/skeletonize_volumetric_volumetric.o $(CXXFLAGS)

obj/skeletonize_path_path.o: src/skeletonize/path/path.cxx
	$(CXX) src/skeletonize/path/path.cxx -c -o obj/skeletonize_path_path.o $(CXXFLAGS)

obj/dispatch_dispatch.o: src/dispatch/dispatch.cxx
	$(CXX) src/dispatch/dispatch.cxx -c -o obj/dispatch_dispatch.o $(CXXFLAGS)

obj/replay_path_path.o: src/replay/path/path.cxx
	$(CXX) src/replay/path/path.cxx -c -o obj/replay_path_path.o $(CXXFLAGS)

clean:
	rm -f obj/*.o lib/libcritter.a lib/libcritter.so
