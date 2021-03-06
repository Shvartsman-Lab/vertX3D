# vertX3D

## Requirements

a) C++ compiler (g++ 5.5.0 was used for development) supporting C++11 <br />
b) CMAKE (version 3.5.1 or newer) <br />
c) ParaView (Interactive 3D viewer for VTK output files generated by vertX3D)

## Getting Started

Follow the following steps to install vertX3D

a) Clone the repository:

    git clone https://github.com/Shvartsman-Lab/vertX3D.git

b) Enter the build directory and use CMAKE to generate a MakeFile:

    cd vertX3D
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release

If instead you are debugging the code, replace "Release" with "Debug". Be warned the code will run slower in debug mode.

c) Compile the examples:
    
    make

## Running Examples

The general form of running the code is as follows

    ./name_of_executable [array_max] [tmax] [seed] [num_threads]

If a segmentation fault is encountered when running the code, it is likely that [array_max] is too small, and an array requires more entries than allocated. [seed] is used to seed the random number generator which for example is utilized in simulations with T1 transformations. [tmax] specifies how long the simulation will run and if [num_threads] is greater than 1 the code will run in parallel.

For the included example simulations, please use the following commands to run them:

### Spreading on a flat surface

    ./vertX3D_spreading_flat 100000 600 1 4

### Spreading on a sphere

    ./vertX3D_spreading_sphere 100000 600 1 4

### Homeostasis

    ./vertX3D_homeostasis 100000 20 1 4

By default, the homeostasis example will run using the overcrowded initial configuration. If you wish to run using the stretched or normal initial configurations, in homeostasis.cpp replace the line

    set_initial_fromFile("./initial/overcrowded.vt3d");

with

    set_initial_fromFile("./initial/stretched.vt3d");

for the stretched configuration or

    set_initial_fromFile("./initial/normal.vt3d");

for the normal configuration. After modifying the code run make in the build folder to re-compile.

### Fluidization

    ./vertX3D_fluidization 200000 2000 1 4

By default, the fluidization example will run the using annealing on the rate of T1 transformations. To run using a fixed rate of T1s, in fluidization.cpp comment out the line:

    T1_spont_act(0.1,300-300*Time/tmax);

and uncomment the line:

    T1_spont_act(0.1,kT1);

kT1 can be set at the top of fluidization.cpp. After modifying the code run make in the build folder to re-compile.

### Sphere crumpling

    ./vertX3D_sphere 300000 100 1 4

## Visualization

The code will output VTK files which can be easily read by ParaView, an interactive visualization tool. ParaView version 5.6.1 was used during development, however newer versions will likely work too. ParaView can be downloaded from https://www.paraview.org/download/. If you are using Linux, make sure to download “ParaView-5.6.1-MPI-Linux-64bit.tar.gz” not “ParaView-5.6.1-osmesa-MPI-Linux-64bit.tar.gz”.

Once you have ParaView launched, output files can be imported using the upper left folder icon or from the File>Open dropdown menu. The output is separated into 3 layers: the apical, basal and lateral sides which can be viewed separately or all at once. Use the animation controls to play your simulation.
