#ifndef TEST__UTILS_H
#define TEST__UTILS_H

#include <chrono>
#include <functional>

std::chrono::microseconds measurePerf(int rep, std::function<void()> fun);

#endif
