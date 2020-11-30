#ifndef __CRITTER_SYMBOL_H__
#define __CRITTER_SYMBOL_H__

#include "../src/intercept/symbol.h"

#define CRITTER_START(ARG)\
  do {\
    critter::internal::symbol_start(#ARG);\
    } while (0);

#define CRITTER_STOP(ARG)\
  do {\
    critter::internal::symbol_stop(#ARG);\
  } while (0);

#define TAU_START(ARG)\
  do {\
    critter::internal::symbol_start(#ARG);\
    } while (0);

#define TAU_STOP(ARG)\
  do {\
    critter::internal::symbol_stop(#ARG);\
  } while (0);

#define TAU_FSTART(ARG)\
  do {\
    critter::internal::symbol_start(#ARG);\
    } while (0);

#define TAU_FSTOP(ARG)\
  do {\
    critter::internal::symbol_stop(#ARG);\
  } while (0);

#define CRITTER_CONDITIONAL_VALUE_CAPTURE_START(ARG)\
  auto _critter_lambda_##ARG = [=]{

#define CRITTER_CONDITIONAL_REFERENCE_CAPTURE_START(ARG)\
  auto _critter_lambda_##ARG = [&]{

#define CRITTER_CONDITIONAL_STOP(ARG)\
  }; critter::internal::symbol_conditional_start(#ARG, _critter_lambda_##ARG);

#endif /*CRITTER_SYMBOL_H_*/
