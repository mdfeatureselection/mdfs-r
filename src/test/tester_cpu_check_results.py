import numpy as np

from config import TESTS, get_results_path
from tester_cpu import run_mdfs


def load_results(name):
    return np.loadtxt(get_results_path(name), dtype=np.dtype("f8"))


def check_results(name, params):
    data = run_mdfs(**params)
    saved_data = load_results(name)

    if (abs(data - saved_data) / saved_data > 0.0002).any():
        print(name + ": ERROR: results differ")
        print(abs(data - saved_data) / saved_data)
    else:
        print(name + ": OK")


for name, params in TESTS.items():
    check_results(name, params)
