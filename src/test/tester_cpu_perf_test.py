from tester_cpu import run_mdfs

params = {
    "variables_count": 5000,
    "objects_count": 100,
    "dimensions": 2,
    "discretizations": 1,
    "divisions": 1,
    "seed": 12345,
    "range": 0.0,
}

run_mdfs(**params)
