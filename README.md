## A further optimization of the Barnes-Hut algorithm using a multi-threaded approach (OpenMP)

The N-body gravitational problem can be naively solved by a simple double loop algorithm performing $$\frac{1}{2} N(N-1)$$ calculations (hence a $$O(n^2)$$ time complexity). A faster (but less accurate) algorithm was proposed by Barnes and Hut
that could solve the N-body problem with $$O(N\log{N})$$ complexity, drastically reducing the number of arithmetic operations required to compute the trajectories of the bodies.

In this repository, part of a final exam for a HPC course,
you will find a parallelized version of the Barnes-Hut algorithm that uses multithreading (OpenMP) to speed up performances


![Recording 2024-12-18 at 20 01 02](https://github.com/user-attachments/assets/f9faca39-266a-4b31-8b5d-0f307fbc838a)

### How to compile and run the program

(N.B. You must have SFML libraries installed in a standard path for the compilation to run smoothly)

It is recommended to build the program using CMake. In the terminal, type:

`$ cmake -S ./src -B ./build`

Once the CMake configuration has finished, build the code:

`$ cmake --build build`

This will automatically create in /build directory an executable file,`barnes-hut.out`. To run the code, simply digit:

`$ build/barnes-hut.out`
