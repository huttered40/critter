#ifndef CRITTER__INTERCEPT__SYMBOL_H_
#define CRITTER__INTERCEPT__SYMBOL_H_

#include <functional>

namespace critter{

void symbol_invoke(const char* symbol, float flops=0, int param1=-1, int param2=-1, int param3=-1, int param4=-1, int param5=-1);

namespace internal{

void symbol_start(const char* symbol);
void symbol_stop(const char* symbol);

void symbol_conditional_start(const char* symbol, std::function<void(void)> f);

};
};

#endif /*CRITTER__INTERCEPT__SYMBOL_H_*/
