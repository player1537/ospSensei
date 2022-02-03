# Sensei Tests

This repository contains code to experiment with interfacing [SENSEI][] with
[OSPRay Studio][] using simple simulations and datasets.

Currently, one simulation and one dataset are planned:
- N-Body Simulation, based on [mini-nbody][].
- HACC data, based on the [Hardware/Hybrid Accelerated Cosmology Code][].

The goal of this repository is to outlive being a "test" repository and provide
an OSPRay / OSPRay Studio implementation of a SENSEI "Analysis Adapter." At that
time, a SENSEI XML configuration file like the following could be expected to
work:

```xml
<sensei>
  <!-- Description:
    - analyse data from a mesh named "mesh" and the array within named "xyz";
    - render it as a point cloud;
    - visualize it in an interactive window.  -->

  <analysis
    type="ospray_studio"
    subtype="interactive"
    mesh="mesh"
    array="xyz"
    association="point"
    enabled="1" />
</sensei>
```

This isn't set in stone.



[SENSEI]: https://github.com/SENSEI-insitu/SENSEI
[OSPRay Studio]: https://github.com/ospray/ospray_studio
[mini-nbody]: https://github.com/harrism/mini-nbody
[Hardware/Hybrid Accelerated Cosmology Code]: https://cpac.hep.anl.gov/projects/hacc/
