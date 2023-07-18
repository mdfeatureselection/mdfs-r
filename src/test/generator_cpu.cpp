#include "generator_cpu.h"

#include <random>

RawData* generateRandom(int v, int o, int d) {
    double* data = new double[v * o];
    int* decision = new int[o];

    std::mt19937 gen(0);
    std::normal_distribution<> norm(0.0f, 1.0f);

    for (int i = 0; i < v * o; ++i) {
        data[i] = norm(gen);
    }

    for (int oi = 0; oi < o; ++oi) {
        bool dec = false;
        for (int i = 0; i < d; ++i) {
            dec = dec ^ (data[i*o + oi] > 0.0f);
        }
        decision[oi] = dec;
    }

    RawData* df = new RawData(RawDataInfo(o,v), data, decision);
    return df;
}
