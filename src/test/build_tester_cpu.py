from pathlib import Path

from cffi import FFI

ffibuilder = FFI()

ffibuilder.cdef(
    """
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
);
"""
)

cpu_sources = list(Path("../cpu").glob("**/*.cpp"))
cpu_tester_sources = ["tester_cpu.cpp", "generator_cpu.cpp", "utils.cpp"]

ffibuilder.set_source(
    "_tester_cpu",
    """
    #include "tester_cpu.h"
    """,
    sources=cpu_tester_sources + cpu_sources,
    extra_compile_args=["-std=c++17", "-O3", "-Wall", "-Wextra", "-pedantic"],
)

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
