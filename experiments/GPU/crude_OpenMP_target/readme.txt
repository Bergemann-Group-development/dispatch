Here be dragons!

This code is not well tested and no fitness for any purpose is
assured. It has been written as a playfield for experimentation with
Dispatch related solvers on modern GPUs.

The code assumes GCC 10.1 with GCCs OpenMP runtime and CUDA 10.2.

Compilation is generally done with:
gfortran -cpp -DNSIZE=16 -O3 -march=native -funsafe-math-optimizations -ffast-math -ffree-line-length-0 -ftree-vectorize -fopenmp -foffload=nvptx-none="-O3" -foffload="-ftree-vectorize" -foffload="-funsafe-math-optimizations" -foffload="-ffast-math" -foffload="-lm" -foffload="-flto" mockup2.f90 muscl.f90 riemann.f90 -lgomp

NSIZE can change to play with different problem sizes.

The number of timesteps is currently hardcoded, and large. It should
be turned into a pre-processor define as well.

With tuning the code can easily leverage all of a Turing class NVIDIA GPU.
