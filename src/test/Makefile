NVCC = nvcc

CXXFLAGS = -O3 -std=c++17 -Wall -Wextra -pedantic -fpic
NVCCFLAGS = -std=c++14 -O3 -arch=compute_30 -ccbin '$(CXX)' -Xcompiler='-fpic'

SOURCES_CPU = $(wildcard ../cpu/*.cpp)
OBJS_CPU = $(SOURCES_CPU:.cpp=.o)
SOURCES_GPU_CXX = $(wildcard ../gpu/*.cpp)
SOURCES_GPU_NVCC = $(wildcard ../gpu/*.cu)
OBJS_GPU_CXX = $(SOURCES_GPU_CXX:.cpp=.o)
OBJS_GPU_NVCC = $(SOURCES_GPU_NVCC:.cu=.o)
OBJS_GPU = $(OBJS_GPU_CXX) $(OBJS_GPU_NVCC)

.PHONY: default all clean

default: tester_gpu

all: tester_gpu tableskernel_unittest

tester_gpu: $(OBJS_GPU) tester_gpu.o utils.o
	$(CXX) $^ -o $@ -lcudart

tableskernel_unittest: $(OBJS_GPU) tableskernel_unittest.o
	$(CXX) $^ -o $@ -lcudart

clean:
	rm -rf *.o $(OBJS_CPU) $(OBJS_GPU) _*

%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@
