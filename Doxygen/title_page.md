\mainpage GRChombo Doxygen
# About GRChombo
GRChombo is a new open-source code for numerical general relativity simulations.
It is developed and maintained by a collaboration of numerical relativists with a
wide range of research interests, from early universe cosmology to astrophysics
and mathematical general relativity, and has been used in many papers since its
first release in 2015.

GRChombo is written entirely in C++14, using hybrid MPI/OpenMP parallelism and
vector intrinsics to achieve good performance on the latest architectures.
Furthermore, it makes use of the Chombo library for adaptive mesh refinement
to allow automatic increasing and decreasing of the grid resolution in regions
of arbitrary shape and topology.

Please visit www.grchombo.org for the full list of developers and their
institutions, a list of publications using GRChombo, and some videos.

# About the GRChombo Doxygen
This doxygen documentation is not intended as standalone documentation of
GRChombo. Please refer to our
[wiki](https://github.com/GRChombo/GRChombo/wiki) for a more pedagogical
introduction to GRChombo. We currently do not guarantee that this doxygen is
kept completely up to date. For the latest version, please generate your own
documentation by following the instructions on 
[this wiki page](https://github.com/GRChombo/GRChombo/wiki/Doxygen-documentation)
