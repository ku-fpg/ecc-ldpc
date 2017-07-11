NVCC_COMPILER_OPTS=-std=c++98

all: cudabits/cached_mult.ptx
	cabal build

cudabits/cached_mult.ptx: cudabits/cached_mult.cu
	(cd cudabits; nvcc cached_mult.cu  --compiler-options "${NVCC_COMPILER_OPTS}" --ptx)

clean:
	rm -f cudabits/cached_mult.ptx
	cabal clean

