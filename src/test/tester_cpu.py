import numpy as np

from _tester_cpu import ffi, lib

stat_mode_map = {
    "H": 0,
    "MI": 1,
    "VI": 2,
}


def run_mdfs(
    variables_count,
    objects_count,
    dimensions,
    discretizations,
    divisions,
    seed,
    range,
    with_decision=True,
    stat_mode="MI",
):
    result = lib.run(
        variables_count,
        objects_count,
        dimensions,
        discretizations,
        divisions,
        seed,
        range,
        with_decision,
        stat_mode_map[stat_mode],
    )
    buf = ffi.buffer(result, variables_count * ffi.sizeof("double"))
    return np.frombuffer(buf, np.dtype("f8"))
