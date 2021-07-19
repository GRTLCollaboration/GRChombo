# Using TwoPunctures initial data with this example

To build this example with [TwoPunctures](https://github.com/GRChombo/TwoPunctures)
initial data, set the environment variable `TWOPUNCTURES_SOURCE` to the Source
subdirectory of the local clone of the [TwoPunctures repository](https://github.com/GRChombo/TwoPunctures)
e.g.
```bash
export TWOPUNCTURES_SOURCE=/path/to/TwoPunctures/Source
```
Then, simply build as normal e.g.
```bash
make all -j 4
```
To stop using TwoPunctures, undefine the `TWOPUNCTURES_SOURCE` environment variable:
```bash
unset TWOPUNCTURES_SOURCE
```
Alternatively, you can avoid defining an environment variable by defining it in
the make command:
```bash
make all -j 4 TWOPUNCTURES_SOURCE=/path/to/TwoPunctures/Source
```
Note that the parameter names for TwoPunctures initial data differ to that of
the vanilla example: see [params_two_punctures.txt](./params_two_punctures.txt)
for the parameter names.
