# GRChombo
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

## Getting started
Detailed installation instructions and usage examples are available in
our [wiki](https://github.com/GRChombo/GRChombo/wiki).

## Contributing
We welcome feedback, bug reports, and contributions. Please consult the [wiki](https://github.com/GRChombo/GRChombo/wiki)
for our coding style and testing policy before filing a pull request.

## License
GRChombo is licensed under the BSD 3-Clause License. Please see LICENSE for details.

## Fawcett Tutorial

The following tutorial has been created to be aimed at Masters/early PhD
students with some HPC experience. We also assume that you have an account set
up on fawcett, and that you are happy with using GitHub. Inevitably it will not be 
possible to cover everything, but if anything is unclear please do get in contact.

To run simulations with GRChombo, first we need to get a copy of both the
Chombo and GRChombo repositories.

### Getting Chombo and GRChombo

The GRChombo copy of the Chombo repository can be found
[here](https://github.com/GRChombo/Chombo). To get a copy, clone it using a
command such as
```bash
git clone https://github.com/GRChombo/Chombo.git ~/Chombo
```
e.g. to clone into your home directory (which we will assume below).

Next, clone the GRChombo repository (and checkout this branch) using a command
such as
```bash
git clone -b training/fawcett https://github.com/GRChombo/GRChombo.git ~/GRChombo
```

Now we can start building the Chombo libraries

### Building Chombo
Before we can build the Chombo libraries, we need to first set the configuration
variables that the build system uses. These are set in the file
`Chombo/lib/mk/Make.defs.local`. This will differ for each system you build on,
but we try to add this file for clusters we have used \[GR\]Chombo on before in
[`GRChombo/InstallNotes/MakeDefsLocalExamples`](./InstallNotes/MakeDefsLocalExamples)
There is a ready-to-use Make.defs.local for use on the `cosmosx` and `skylake`
partitions on fawcett (where `cosmosx` is a shared memory node with skylake
cores, see
[here](https://www.maths.cam.ac.uk/computing/faculty-hpc-system-fawcett)). Copy this to the right place in your copy of the Chombo
repository using the command:
```bash
cp ~/GRChombo/InstallNotes/MakeDefsLocalExamples/fawcett-Skylake.Make.defs.local \
~/Chombo/lib/mk/Make.defs.local
```
For more details on what the variables are in the Make.defs.local file, have a
look at
[this wiki page](https://github.com/GRChombo/GRChombo/wiki/Compiling-Chombo).

Next, we need to load the necessary modules. We will use the Intel compiler with
Intel MPI and a compatible HDF5 module. This can either be done using the
script in
[`GRChombo/InstallNotes/MakeDefsLocalExamples`](./InstallNotes/MakeDefsLocalExamples)
with the command:
```bash
source ~/GRChombo/InstallNotes/MakeDefsLocalExamples/fawcett-modules.sh
```
or by loading the modules explicitly with the commands:
```bash
module load intel/compilers/2018.3
module load intel/impi/2018.3/intel
module load hdf5-intel/1.10.4
```
Note that you should also check using 
```bash 
module list
```
that no conflicting
modules are loaded, as this can cause problems with compilation. If they are,
you can remove them using 
```bash
module unload <module_name>
```

Finally we can build the Chombo libraries with
```bash
cd ~/Chombo/lib
make lib -j <num tasks>
```
Note that a sensible number of tasks will depend on whether you are building
on the login node (in which case choose a small number such as 4) or if you are
in an interactive job (which case feel free to use up to about twice
as many tasks as you have cores available).

For more details on building Chombo (including what the variables are in the
`Make.defs.local` file) have a look at
[this wiki page](https://github.com/GRChombo/GRChombo/wiki/Compiling-Chombo).

Now we are ready to build a GRChombo example.

### Building the GRChombo BinaryBH example
First we need to let GRChombo know where to look for the Chombo libraries.
This is done with the `CHOMBO_HOME` environment variable and can be set using
the command
```bash
export CHOMBO_HOME=${HOME}/Chombo/lib
```
You might want to add this line to your `~/.bashrc` file so that every time you
start a new shell, the environment variable is automatically set.

Now we are ready to build with the commands
```bash
cd ~/GRChombo/Examples/BinaryBH
make all -j <num tasks>
```
Assuming everything compiled and linked sucessfully, you should have an
executable called something like (the name will change depending on the Chombo
configuration variables)
```
Main_BinaryBH3d.Linux.64.mpiicpc.ifort.OPTHIGH.MPI.OPENMPCC.Skylake.ICC2018.3.ex
```

### Running the example
In general, on most computing clusters, you should not run jobs from your home
directory.
Usually there is some separate parallel storage to be used for simulations
(often called the "scratch" partition) and you should consult the documentation
of the particular cluster on where to find it. On fawcett, this is
the `/nfs/st01/hpc-gr-epss` directory. First, let's create a subdirectory
for our simulation
```bash
mkdir -p /nfs/st01/hpc-gr-epss/${USER}/BinaryBH/MyFirstSim
```
and let's copy our binary across to the parent of this new directory using
a command such as
```bash
cp ~/GRChombo/Examples/BinaryBH/Main_BinaryBH3d.Linux.64.mpiicpc.ifort.OPTHIGH.MPI.OPENMPCC.Skylake.ICC2018.3.ex \
/nfs/st01/hpc-gr-epss/${USER}/BinaryBH/
```
I have created a parameter file to be used for today. It contains initial data
for the quasicircular inspiral of an equal-mass BH binary which should do
about 10 orbits (though due to low resolution and approximate initial data,
it will probably be a fair number fewer - it would also take around 2 days to
finish). Let's copy it across to our simulation
directory that we created:
```bash
cp ~/GRChombo/Examples/BinaryBH/params_fawcett.txt /nfs/st01/hpc-gr-epss/${USER}/BinaryBH/MyFirstSim/params.txt
```
We will need to create a jobscript to submit to the scheduler. I made an example
template in `/nfs/st01/hpc-gr-epss/GRChombo-Training/example-jobscript` which
you can copy over to your simulation directory with the command:
```bash
cp /nfs/st01/hpc-gr-epss/GRChombo-Training/example-jobscript /nfs/st01/hpc-gr-epss/${USER}/BinaryBH/MyFirstSim/jobscript
```
You will need to edit a few lines in this template before you can submit it.
Use your favourite available text editor to provide the full path to the
executable i.e. replace the line
```bash
EXEC="/full/path/to/binary.ex"
```
with
```bash
EXEC="/nfs/st01/hpc-gr-epss/<crsdid>/BinaryBH/Main_BinaryBH3d.Linux.64.mpiicpc.ifort.OPTHIGH.MPI.OPENMPCC.Skylake.ICC2018.3.ex"
```
where `<crsid>` is your username. We may also need to change the reservation
name if it has changed. Feel free to change the time limit. Once you are happy,
save and exit your editor and submit it to the scheduler using
```bash
cd /nfs/st01/hpc-gr-epss/${USER}/BinaryBH/MyFirstSim/ && sbatch jobscript
```
Assuming there are available resources, the job should start shortly.
You can query the status of the job with the command
```bash
squeue -u $USER
```
Within the first minute of the job starting, you should be able to see the
output from each MPI rank in `pout/pout.n` where `n` is the rank number.
To monitor the output of this, you can use the command
```bash
tail -f pout/pout.0
```
Note that there isn't much output during the initial setup after the loading of
parameters so don't panic if it appears nothing is happening for a bit.
You should also check the main job output in the `slurm-<jobid>.out` file
to check that it hasn't crashed.

Congratulations, you have just run your first GRChombo simulation!

### Visualisation

In this section, we will outline how to get set up with visualisation using
Paraview on Fawcett. It is a bit fiddly, so we will go over the bare minimum
to get up and running. For more detailed information on what Paraview can do,
have a look at the
[tutorial](https://www.paraview.org/Wiki/images/b/bc/ParaViewTutorial56.pdf)
(or ask Amelia!). Adapted from instructions provided by Kacper Kornet.

First, we need to set up ssh access to fawcett. Information on how to do this
is given [here](https://www.maths.cam.ac.uk/computing/fawcett-ssh-access). We
have sent this to you in advance, so you should have all of the steps set up
prior to the tutorial.

Second, we need you to download Paraview v5.6.0 (note that is important that
you download this specific version and not a different one) onto your local machine i.e. the
laptop/desktop you are using to access fawcett. This can be done
[here](https://www.paraview.org/download/). Please do this in advance of the tutorial.

Before we tunnel to fawcett, we need to add a server on the Paraview
client. To do this, open Paraview, and click the connect icon near the top left corner,
then 'Add Server'. Change the name to whatever you like e.g. 'localhost'. All
other settings should be the default:  

Server type: Client/Server  
Host: localhost  
Port: 11111  

Click Configure. On the next screen, set 

Startup Type: Manual

Accept by clicking Save and Close. Once you have set this up, you won't need
to do it again.

Now we need to start the Paraview server on fawcett. Ssh into fawcett - as
outlined above, this should already be configured with ssh forwarding, and so
shouldn't require explicit login via another machine. Load Paraview using

```bash
module load paraview/5.6.0/upstream
```

Now run the command
```bash
pvserver
```

It should print something like:

```bash
Waiting for client...
Connection URL: cs://fawcett:11111
Accepting connection(s): fawcett:11111
```

Note that the number 11111 can be different if someone else
is already running pvserver. In what follows, we will refer to it as `<fawcett_port>`.

(For large datasets, it is a good idea to submit a job to the queue and run
pvserver with a modified environment

```bash
OSPRAY_THREADS=8 KNOB_MAX_WORKER_THREADS=8 pvserver
```

where 8 is the number of cores reserved for the job. However for the purposes
of this tutorial we should be fine without.)

We next need to set up the tunnel to fawcett. On your local machine, run

```bash
ssh -NL 11111:localhost:cosmos_port username@fawcett.maths.cam.ac.uk
```

It should ask for authentication, then look as though it has 'hung up'.

On your local Paraview client, again click 'connect' and choose the connection
configured in the earlier step. This should now have opened a connection to
fawcett that you can use to access your data using `File->Open`. We will go
over any smaller points in the tutorial if necessary. You should now be able
to see your two black holes!

For more detailed information, consult the
[wiki](https://github.com/GRChombo/GRChombo/wiki/).
