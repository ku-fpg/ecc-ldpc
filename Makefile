NVCC_COMPILER_OPTS=#-std=c++98
ARCH_FLAGS=-arch=sm_37 -gencode=arch=compute_37,code=sm_37 -rdc=true -lcudadevrt
NVCC_OPTS= --ptx -g -dc

all: cudabits/cached_mult.ptx cudabits/arraylet1.ptx
	cabal build

cudabits/cached_mult.ptx: cudabits/cached_mult.cu cudabits/common.h
	(cd cudabits; nvcc cached_mult.cu ${ARCH_FLAGS} --compiler-options "${NVCC_COMPILER_OPTS}" ${NVCC_OPTS})

cudabits/arraylet1.ptx: cudabits/arraylet1.cu cudabits/common.h
	(cd cudabits; nvcc arraylet1.cu ${ARCH_FLAGS} --compiler-options "${NVCC_COMPILER_OPTS}" ${NVCC_OPTS})
clean:
	rm -f cudabits/cached_mult.ptx
	rm -f cudabits/arraylet1.ptx
	cabal clean

