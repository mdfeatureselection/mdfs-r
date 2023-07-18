extern "C" {
    #include "tester_cpu.h"
}

#include "utils.h"
#include "generator_cpu.h"
#include "../cpu/mdfs.h"

#include <iostream>

double* run(
        int var_count,
        int obj_count,
        int dim,
        int dis,
        int div,
        uint32_t seed,
        float range,
        int with_decision,
        int stat_mode
) {
    int rep = 1;

    RawData *df = generateRandom(var_count, obj_count, dim);

    std::unique_ptr<const DiscretizationInfo> dfi(new DiscretizationInfo(seed, dis, div, range));

    MDFSInfo mdfs_info(
        dim,
        div,
        dis,
        0.25f,
        0.0f,
        NULL,
        0, // required for proper functioning
        true,
        nullptr,
        false
    );

    MDFSOutput out(MDFSOutputType::MaxIGs, dim, var_count, 0);

    // NOTE: done to avoid defining the function type
    auto mdfsFun = mdfs[dim-1];

    if (with_decision) {
        switch (stat_mode) {
            case StatMode::Entropy:
                mdfsFun = mdfsDecisionConditionalEntropy[dim-1];
                break;
            case StatMode::MutualInformation:
                mdfsFun = mdfs[dim-1];
                break;
            case StatMode::VariationOfInformation:
                mdfsFun = mdfsDecisionConditionalVariationOfInformation[dim-1];
                break;
        }
    } else {
        switch (stat_mode) {
            case StatMode::Entropy:
                mdfsFun = mdfsEntropy[dim-1];
                break;
            case StatMode::MutualInformation:
                mdfsFun = mdfsMutualInformation[dim-1];
                break;
            case StatMode::VariationOfInformation:
                mdfsFun = mdfsVariationOfInformation[dim-1];
                break;
        }
    }

    std::chrono::microseconds time = measurePerf(rep, [&]() -> void { mdfsFun(mdfs_info, df, nullptr, std::move(dfi), out); });;

    std::cout << "Procedure time " << time.count() << " us" << "\n";

    double* igs = new double[var_count];
    out.copyMaxIGsAsDouble(igs);

    delete df->data;
    delete df->decision;
    delete df;

    return igs;
}
