# GRChombo
GRChombo is an open-source code for numerical general relativity simulations.
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

## Getting started
Detailed installation instructions and usage examples are available in
our [wiki](https://github.com/GRChombo/GRChombo/wiki), with the home page giving guidance on where to start.

## Contributing
We welcome feedback, bug reports, and contributions. Please consult the [wiki](https://github.com/GRChombo/GRChombo/wiki)
for our coding style and testing policy before filing a pull request.

## License
GRChombo is licensed under the BSD 3-Clause License. Please see LICENSE for details.
