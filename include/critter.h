#ifndef __CRITTER_H__
#define __CRITTER_H__

#include <mpi.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <utility>
#include <iomanip>
#include <vector>
#include <stack>
#include <stdint.h>
#include <functional>
#include <map>
#include <set>
#include <unordered_map>
#include <bitset>
#include <cmath>
#include <assert.h>

namespace critter{

	// alpha-beta butterfly, BSP
	// BSPcommCost,ABCommCost,BSPsynchCost,ABSynchCost,CommTime,SynchTime,DataMvtTime,CompTime,RunTime
void start(size_t mode = 1, size_t mechanism = 0, const char* _cost_models_ = "11", const char* _breakdown_ = "000000001", size_t _max_num_symbols_ = 75, size_t _max_timer_name_length_ = 40);
void stop(size_t mode = 1, size_t mechanism = 0, size_t factor = 1);
}

#include "../src/critter.hpp"

#endif /*CRITTER_H_*/
