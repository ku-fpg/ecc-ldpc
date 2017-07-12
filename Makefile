NVCC_COMPILER_OPTS=-std=c++98
ARCH_FLAGS=-arch=sm_35 -gencode=arch=compute_35,code=sm_35 -rdc=true -lcudadevrt

all: cudabits/cached_mult.ptx
	cabal build

cudabits/cached_mult.ptx: cudabits/cached_mult.cu
	(cd cudabits; nvcc cached_mult.cu ${ARCH_FLAGS} --compiler-options "${NVCC_COMPILER_OPTS}" --ptx)

clean:
	rm -f cudabits/cached_mult.ptx
	cabal clean

