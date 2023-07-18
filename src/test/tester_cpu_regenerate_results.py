import numpy as np

from config import TESTS, get_results_path
from tester_cpu import run_mdfs


def save_results(name, data):
    np.savetxt(get_results_path(name), data)


for name in TESTS:
    data = run_mdfs(**TESTS[name])
    save_results(name, data)
