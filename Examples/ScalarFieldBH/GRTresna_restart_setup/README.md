## Restarting from GRTresna and convergence testing 

GRTresna is the GRTL Collaboration code for solving the initial data constraints of GR. See [the GRTresna repository](https://github.com/GRTLCollaboration/GRTresna).

These parameters relate to the GRTresna [ScalarFieldBH example](https://github.com/GRTLCollaboration/GRTresna/tree/main/Examples/ScalarFieldBH). One should run the solver there, and use the output as a restart file for this evolution example. We also provide an example of how to perform a convergence test on the outputs, which is essential for verifying that things are working correctly. For more details see the [GRTresna wiki page](https://github.com/GRTLCollaboration/GRTresna/wiki/Convergence-testing).