#include "utils.h"
#include "../gpu/datafile.h"
#include "../gpu/cucubes.h"

#include <iostream>
#include <random>

const int SEED = 1234;

void genXOR(InputFile inf, int n, int* vars) {
	std::mt19937 gen(SEED);
	std::uniform_real_distribution<double> dis(-1.0, 1.0);

	for (int i = 0; i < inf.vars; i++) {
		for (int j = 0; j < inf.objs; j++) {
			inf.data[i * inf.objs + j] = dis(gen);
		}
	}

	for (int i = 0; i < inf.objs; i++) {
		inf.decision[i] = 0;
		for (int j = 0; j < n; j++) {
			inf.decision[i] ^= (inf.data[vars[j] * inf.objs + i] > 0.0) ? 1 : 0;
		}
	}
}

int main(int argc, char** argv) {
	if (argc < 6) {
		std::cerr << "Wrong number of arguments, expected 5: var obj dis dim div\n";
		return 1;
	}
	int VARS = atoi(argv[1]);
	int OBJS = atoi(argv[2]);
	int DISC = atoi(argv[3]);
	int DIM = atoi(argv[4]);
	int DIV = atoi(argv[5]);

	int n = 2;
	int vars[2] = {0, 1};

	InputFile inf;
	inf.vars = VARS;
	inf.objs = OBJS;
	inf.data = new double[VARS * OBJS];
	inf.decision = new int[OBJS];

	double* IGmax = new double[VARS];

	auto time1 = measurePerf(1, [&]() -> void {
		genXOR(inf, n, vars);
	});

	std::cout << "Generation time " << time1.count() << " us" << "\n";

	auto time2 = measurePerf(1, [&]() -> void {
		run_cucubes(OBJS, VARS, DIM, DIV, DISC, SEED, 0.5, 0.1, inf.data, inf.decision, IGmax);
	});

	std::cout << "Procedure time " << time2.count() << " us" << "\n";

	delete[] inf.data;
	delete[] inf.decision;
	delete[] IGmax;
	return 0;
}
