#ifndef CRITTER__RECORD__RECORD_H_
#define CRITTER__RECORD__RECORD_H_

#include "../util.h"

namespace critter{
namespace internal{

void record(std::ofstream& Stream);
void record(std::ostream& Stream);

}
}

#endif /*CRITTER__RECORD__RECORD_H_*/
