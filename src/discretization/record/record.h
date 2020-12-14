#ifndef CRITTER__DISCRETIZATION__RECORD__RECORD_H_
#define CRITTER__DISCRETIZATION__RECORD__RECORD_H_

#include "../../util/util.h"

namespace critter{
namespace internal{
namespace discretization{

class record{
public:
  static void write_file(int variantID, int print_mode, float overhead_time);
  static void print(int variantID, int print_mode, float overhead_time);

private:
  static void set_tuning_statistics();
  static void set_kernel_statistics();
};

}
}
}

#endif /*CRITTER__DISCRETIZATION__RECORD__RECORD_H_*/
