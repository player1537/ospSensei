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


## Building: On Deadpool

```console
$ # Download dependencies
$ ./go.sh ospray git clone
$ ./go.sh vtk git clone
$ ./go.sh sensei git clone
$
$ # Build and install OSPRay
$ ./go.sh ospray cmake configure
$ ./go.sh ospray cmake parbuild
$ ./go.sh ospray cmake install
$
$ # Build and install VTK
$ ./go.sh -ospray vtk cmake configure
$ ./go.sh -ospray vtk cmake parbuild
$ ./go.sh -ospray vtk cmake install
$
$ # Build and install SENSEI
$ ./go.sh -ospray -vtk sensei cmake configure
$ ./go.sh -ospray -vtk sensei cmake parbuild
$ ./go.sh -ospray -vtk sensei cmake install
$
$ # Build and install this library
$ ./go.sh -ospray -vtk -sensei cmake configure
$ ./go.sh -ospray -vtk -sensei cmake build
$ ./go.sh -ospray -vtk -sensei cmake install
```

**Testing**:

```console
$ # Run nbody example app
$ rm -f ospSensei.*.ppm
$ ./go.sh -ospray -vtk -sensei exec mpirun -np 1 nbody
$ # Might segfault due to faulty MPI finalization
$ ls ospSensei.*.ppm
$ # Expect: 10 files
$
$ # Run mandelbrot example app
$ rm -f ospSensei.*.ppm
$ ./go.sh -ospray -vtk -sensei exec mpirun -np 1 mandelbrot
$ # Might segfault due to faulty MPI finalization
$ ls ospSensei.*.ppm
$ # Expect: 1 file
```

## Usage: Examples

Two example apps are included within the repository. These are in
[nbody.cpp][app/nbody/nbody.cpp] and
[mandelbrot.cpp][app/mandelbrot/mandelbrot.cpp].

In contrast with SENSEI, this library does not (yet) provide an equivalent to
SENSEI's `ConfigurableAnalysis` that parses an XML file and changes the
analysis/data adaptors at runtime. This means that the specific `ospSensei`
class needs to be referenced directly.

In the simplest case, the analysis adaptors are used as follows:

```c++
#include <vtkNew.h>
#include <vtkPolyData.h>
#include <VTKDataAdaptor.h>
#include <ospSensei/OSPRayParticleVisualization.h>

int main() {
    // ...

    using AnalysisAdaptor = ospSensei::OSPRayParticleVisualization;
    vtkNew<AnalysisAdaptor> analysis;
    analysis->SetCommunicator(MPI_COMM_WORLD);
    analysis->SetMeshName("bodies");
    analysis->SetArrayName("position");
    analysis->SetWidth(512);
    analysis->SetHeight(512);
    analysis->Initialize();

    for (size_t iter=0; iter<n_iters; ++iter) {
        // ...

        vtkNew<vtkPolyData> polyData;

        using DataAdaptor = sensei::VTKDataAdaptork;
        vtkNew<DataAdaptor> data;
        data->SetDataObject("bodies", polyData);

        analysis->Execute(data);

        data->ReleaseData();
    }

    analysis->Finalize();

    // ...
}
```

```c++
#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>
#include <VTKDataAdaptor.h>
#include <ospSensei/OSPRayUnstructuredVolumeVisualization.h>

int main() {
    // ...

    using AnalysisAdaptor = ospSensei::OSPRayUnstructuredVolumeVisualization;
    vtkNew<AnalysisAdaptor> analysis;
    analysis->SetCommunicator(MPI_COMM_WORLD);
    analysis->SetMeshName("mandelbrot");
    analysis->SetArrayName("nsteps");
    analysis->SetWidth(512);
    analysis->SetHeight(512);
    analysis->Initialize();

    for (size_t iter=0; iter<n_iters; ++iter) {
        // ...

        vtkNew<vtkUnstructuredGrid> ugrid;

        using DataAdaptor = sensei::VTKDataAdaptor;
        vtkNew<DataAdaptor> data;
        data->SetDataObject("mandelbrot", ugrid);

        analysis->SetUseD3(true); // or false
        analysis->Execute(data);

        data->ReleaseData();
    }

    analysis->Finalize();

    // ...
```


[SENSEI]: https://github.com/SENSEI-insitu/SENSEI
[OSPRay]: https://github.com/ospray/ospray
