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
#include <array>
#include <stdint.h>
#include <functional>
#include <map>
#include <unordered_map>
#include <bitset>
#include <cmath>
#include <assert.h>

namespace critter{

// *****************************************************************************************************************************************************************
// User functions
void start(size_t mode = 1);
void stop(size_t mode = 1, size_t factor = 1);

// User variables
// Note: `cost_model_size` must equal `cost_models.count()`. This will not be checked at compile time.
constexpr size_t cost_model_size  = 2;					// must match number of bits set in cost_model (below)
constexpr std::bitset<2> cost_models(0b11);				// alpha-beta butterfly, BSP
// Note: `breakdown_size` must equal `breakdown.count()`. This will not be checked at compile time.
constexpr size_t breakdown_size  			= 5;			// must match number of bits set in breakdown (below)
constexpr std::bitset<5+2*cost_model_size> breakdown(0b110010101);  		// RunTime,CompTime,DataMvtTime,SynchTime,CommTime,ABSynchCost,BSPsynchCost,ABCommCost,BSPcommCost
constexpr int internal_tag                      	= 1669220;		// arbitrary
constexpr int internal_tag1                     	= 1669221;		// arbitrary
constexpr int internal_tag2                     	= 1669222;		// arbitrary
constexpr int internal_tag3                     	= 1669223;		// arbitrary
constexpr int internal_tag4                     	= 1669224;		// arbitrary
constexpr int internal_tag5                     	= 1669225;		// arbitrary
constexpr size_t max_timer_name_length 			= 40;			// max length of a symbol defining a timer
constexpr size_t max_num_symbols       			= 75;			// max number of symbols to be tracked

// *****************************************************************************************************************************************************************
}

#include "../src/critter.hpp"

#endif /*CRITTER_H_*/
