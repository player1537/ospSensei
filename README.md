# OSPRay + SENSEI Integration

This repository contains a library that adds [OSPRay][] scientific visualization
support to the [SENSEI][] library. The primary goal of this library is to enable
developers to easily integrate multiprocess simulations with multiprocess
rendering using OSPRay's MPI module.


## Building

This library expects to have MPICH, OSPRay, VTK, and SENSEI installed already
and expects them to be findable via CMake's `find_package()` command.

For convenience, a script at the root level helps manage the development
environment this code expects. Running the following commands should yield a
suitable development environment, however these commands aren't necessary. Any
method of installing the required dependencies will work.

```console
$ # Download dependencies
$ ./go.sh ospray git clone
$ ./go.sh vtk git clone
$ ./go.sh sensei git clone
$
$ # Setup Docker environment
$ ./go.sh docker build
$ ./go.sh docker start
$
$ # Build and install OSPRay
$ ./go.sh -docker ospray cmake configure
$ ./go.sh -docker ospray cmake parbuild
$ ./go.sh -docker ospray cmake install
$
$ # Build and install VTK
$ ./go.sh -docker -ospray vtk cmake configure
$ ./go.sh -docker -ospray vtk cmake parbuild
$ ./go.sh -docker -ospray vtk cmake install
$
$ # Build and install SENSEI
$ ./go.sh -docker -ospray -vtk sensei cmake configure
$ ./go.sh -docker -ospray -vtk sensei cmake parbuild
$ ./go.sh -docker -ospray -vtk sensei cmake install
$
$ # Build and install this library
$ ./go.sh -docker -ospray -vtk -sensei cmake configure
$ ./go.sh -docker -ospray -vtk -sensei cmake build
$ ./go.sh -docker -ospray -vtk -sensei cmake install
```


[SENSEI]: https://github.com/SENSEI-insitu/SENSEI
[OSPRay]: https://github.com/ospray/ospray
