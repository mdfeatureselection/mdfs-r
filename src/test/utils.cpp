#include "utils.h"

std::chrono::microseconds measurePerf(int rep, std::function<void()> fun)
{
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    for (int i = 0; i < rep; ++i) {
        fun();
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    return (std::chrono::duration_cast<std::chrono::microseconds>(end - start));
}
