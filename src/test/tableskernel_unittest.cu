#include <stdio.h>
#include <math.h>
#include <stdint.h>

#include "../gpu/tableskernel.cuh"

#define TILE_SIZE 16
#define DIM 4
#define DIV 3
#define BITS 2

void cpuKernel(KernelParam param) {
	int dim = param.dim - 1;
	int varlen = param.packs[0] + param.packs[1];
	int ncnt = (int)pow(param.div + 1, dim);
	int ntup = (int)pow(param.tileSize, dim);
	int bits = (int)std::ceil(std::log2f((float)(param.div + 1)));

	uint32_t counters[ncnt][2];

	uint64_t* data[dim];
	for (int i = 0; i < dim; i++) {
		data[i] = param.data[i >= dim - param.index ? i + 1 : i];
	}

	// po dyskretyzacjach
	for (int r = 0; r < param.disc; r++) {
		printf("d%d\n", r);
		// po krotkach
		for (int t = 0; t < ntup; t++) {
			for (int cnt = 0; cnt < ncnt; cnt++) {
				counters[cnt][0] = counters[cnt][1] = 0;
			}

			// po decyzjach
			for (int d = 0; d < 2; d++) {
				// po paczkach
				int done = 0;
				for (int p = d ? param.packs[0] : 0;
					p < param.packs[d] + (d ? param.packs[0] : 0);
					p += bits) {
					// po obiektach
					for (int o = 0; done + o < param.objs[d] && o < 64; o++) {
						int pos = 0;
						int shift = 1;

						int tt = t;
						// po współrzędnych
						for (int v = 0; v < dim; v++) {
							int val = 0;
							int var = tt % param.tileSize;
							tt /= param.tileSize;
							// po bitach
							for (int b = 0; b < bits; b++) {
								val |= ((data[v][var * varlen + p + b] >> o) & 1) << b;
							}
							pos += shift * val;
							shift *= param.div + 1;
						}

						counters[pos][d]++;
					}
					done += 64;
				}
			}

			// memcpy? bigendian.
			for (int cnt = 0; cnt < ncnt; cnt++) {
				param.counters[param.index][t * ncnt + cnt] =
					(uint64_t)counters[cnt][0] << 32 | (uint64_t)counters[cnt][1];
				//printf(">>> %x %x %llx\n", (uint64_t)counters[cnt][0], (uint64_t)counters[cnt][1],
				//	param.counters[param.index][t * ncnt + cnt]);
			}
		}

		param.counters[param.index] += ntup * ncnt;
		for (int i = 0; i < dim; i++) {
			data[i] += param.tileSize * varlen;
		}
	}

}

void generateInput(uint64_t* data, int div, int bits, int disc, int vars, int obj0, int obj1) {
	int objs[2] = {obj0, obj1};
	int off[2] = {0, bits * ((obj0 + 63) / 64)};

	int varlen = bits * (((obj0 + 63) / 64) + ((obj1 + 63) / 64));
	int disclen = varlen * vars;

	for (int r = 0; r < disc; r++) {
		for (int v = 0; v < vars; v++) {
			for (int d = 0; d < 2; d++) {
				for (int o = 0; o < objs[d]; o++) {
					int val = rand() % (div + 1);
					for (int b = 0; b < bits; b++) {
						data[r * disclen + v * varlen + off[d] + bits * (o / 64) + b] |=
							(uint64_t)((val >> b) & 1) << o;
					}
				}
			}
		}
	}
}


int main() {
	KernelParam input(TILE_SIZE, DIM, DIV, RM_AVG, BF_SPLIT, true, 0);
	printf("input init... Done.\n");

	input.disc = 3;

	//input.vars = 12345; // WSZYSTKIE zmienne

	// Ile bitów na obiekt?
	int bits = (int)std::ceil(std::log2f((float)(input.div + 1)));
	printf("bits: %d, BITS: %d\n", bits, BITS);

	// Nieistotne:
	// Alokujemy dopełnione do tileSize
	//std::size_t len = (input.vars + input.tileSize - 1) / input.tileSize * input.tileSize;
	//cudaMalloc(&input.IG, len * sizeof(float));
	//cudaMemset(input.IG, 0, len * sizeof(float));

	// Nieistotne:
	//input.pseudo[0] = 10.0;
	//input.pseudo[1] = 10.0;

	input.objs[0] = 12342;
	input.objs[1] = 33312;

	for (int dec = 0; dec < 2; dec++) {
		input.packs[dec] = ((input.objs[dec] + 63) / 64) * bits;
		printf("input.packs[%d] = %d\n", dec, input.packs[dec]);
	}

	uint64_t vol = 1;
	for (int i = 1; i < input.dim; i++) {
		vol *= input.tileSize * (input.div + 1);
	}
	printf("vol: %d\n", vol);

	printf("cudaMalloc...");
	uint64_t varlen = input.packs[0] + input.packs[1];
	printf("varlen: %d\n", varlen);
	uint64_t vars = input.disc * input.tileSize;
	for (int i = 0; i < input.dim; i++) {
		//input.offset[i] = 0;
		cudaMalloc(&input.data[i], vars * varlen * sizeof(uint64_t));
		cudaMalloc(&input.counters[i], input.disc * vol * sizeof(uint64_t));
	}
	printf(" Done\n");

	printf("generateInput...");
	KernelParam refInput = input;
	for (int i = 0; i < refInput.dim; i++) {
		refInput.data[i] = new uint64_t[vars * varlen];
		memset(refInput.data[i], 0, vars * varlen * sizeof(uint64_t));
		refInput.counters[i] = new uint64_t[input.disc * vol];

		generateInput(refInput.data[i], input.div, bits, refInput.disc, refInput.tileSize, input.objs[0], input.objs[1]);
		cudaMemcpy(input.data[i], refInput.data[i], vars * varlen * sizeof(uint64_t), cudaMemcpyHostToDevice);
	}
	printf(" Done\n");

	// RUN
	uint64_t* buffer = new uint64_t[input.disc * vol];
	for (int ix = 0; ix < input.dim; ix++) {
		refInput.index = ix;
		input.index = ix;

		cpuKernel(refInput);
		tablesKernelWrapper<TILE_SIZE, DIM - 1, DIV, BITS, 1>(input, 0);

		printf("%d: ", ix);
		cudaMemcpy(buffer, input.counters[ix], input.disc * vol * sizeof(uint64_t), cudaMemcpyDeviceToHost);
		for (int i = 0; i < input.disc * vol; i++) {
			if (buffer[i] != refInput.counters[ix][i]) {
				printf("FAIL(%d %d) %lx %lx\n", ix, i, buffer[i], refInput.counters[ix][i]);
				exit(1);
			}
		}
		printf("OK\n");
	}

	return 0;
}
