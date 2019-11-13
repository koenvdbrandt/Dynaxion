# TCAD DF-ISE mesh converter
This code takes the `.grd` and `.dat` files of the DF-ISE format from TCAD simulations as input. The `.grd` file contains the vertex coordinates (3D or 2D) of each mesh node and the `.dat` file contains the value of each electric field vector component for each mesh node, grouped by model regions (such as silicon bulk or metal contacts). The regions are defined in the `.grd` file by grouping vertices into edges, faces and, consecutively, volumes or elements.

A new regular mesh is created by scanning the model volume in regular X Y and Z steps (not necessarily coinciding with original mesh nodes) and using a barycentric interpolation method to calculate the respective electric field vector on the new point. The interpolation uses the four closest, no-coplanar, neighbor vertex nodes such, that the respective tetrahedron encloses the query point. For the neighbors search, the software uses the Octree implementation [@octree].

The output `.init` or `.apf` file can be imported into Allpix Squared. The INIT file is an ASCII text file with a header followed by a list of columns organized as
```bash
node.x node.y node.z observable.x observable.y observable.z
```
The APF (Allpix Squared Field) data format contains the field data in binary form and is therefore a bit more compact and can be read much faster. Whenever possible, this format should be preferred.

### Compilation

When compiling the Allpix Squared framework, the TCAD DF-ISE mesh converter is automatically compiled and installed in the Allpix Squared installation directory.

It is also possible to compile the converter separately as stand-alone tool within this directory:
```bash
$ mkdir build
$ cd build
$ cmake ..
$ make
```

It should be noted that the TCAD DF-ISE mesh converter depends on the core utilities of the Allpix Squared framework found in the directory `src/core/utils`. Thus, it is discouraged to move the converter code outside the repository as this directory would have to be copied and included in the code as well. Furthermore, updates are only distributed through the repository and new release versions of the Allpix Squared framework.

### Features
- TCAD DF-ISE file format reader.
- Fast radius neighbor search for three-dimensional point clouds.
- Barycentric interpolation between non-regular mesh points.
- Several cuts available on the interpolation algorithm variables.
- Interpolated data visualization tool.

### Parameters
* `model`: Field file format to use, can be **INIT** or **APF**, defaults to **APF** (binary format).
* `dimension`: Specify mesh dimensionality (defaults to 3).
* `region`: Region name or list of region names to be meshed (defaults to `bulk`).
* `observable`: Observable to be interpolated (defaults to `ElectricField`).
* `initial_radius`: Initial node neighbors search radius in micro meters. Defaults to the minimal cell dimension of the final interpolated mesh.
* `radius_step`: Radius step if no neighbor is found (defaults to `0.5um`).
* `max_radius`: Maximum search radius (default is `10um`).
* `volume_cut`: Minimum volume for tetrahedron for non-coplanar vertices (defaults to minimum double value).
* `divisions`: Number of divisions of the new regular mesh for each dimension, 2D or 3D vector depending on the `dimension` setting. Defaults to 100 bins in each dimension.
* `xyz`: Array to replace the system coordinates of the mesh. A detailed description of how to use this parameter is given below.
* `mesh_tree`: Boolean to enable creation of a root file with the TCAD mesh nodes stored in a `ROOT::TTree`. This setting is deactivated by default.
* `workers`: Number of worker threads to be used for the interpolation. Defaults to the available number of cores on the machine (hardware concurrency).

### Usage
To run the program, the following command should be executed from the installation folder:
```bash
mesh_converter -f <file_prefix> [<options>] [<arguments>]
```
The converter will look for a configuration file with `<file_prefix>` and `.conf` extension. This default configuration file name can be replaced with the `-c` option.
The list with options can be accessed using the `-h` option.
Possible options and their default values are:
```
-f <file_prefix>       common prefix of DF-ISE grid (.grd) and data (.dat) files
-c <config_file>	   configuration file setting mesh conversion parameters
-h                     display this help text
-l <file>              file to log to besides standard output (disabled by default)
-o <init_file_prefix>  output file prefix without .init (defaults to file name of <file_prefix>)
-v <level>             verbosity level (default reporting level is INFO)
```

Observables currently implemented for interpolation are: `ElectrostaticPotential`, `ElectricField`, `DopingConcentration`, `DonorConcentration` and `AcceptorConcentration`.
The output INIT/APF file will be saved with the same file_prefix as the `.grd` and `.dat` files and the additional name suffix `_<observable>_interpolated` and the appropriate file extension, where `<observable>` is replaced with the selected quantity.

The new coordinate system of the mesh can be changed by providing an array for the *xyz* keyword in the configuration file. The first entry of the array, representing the new mesh *x* coordinate, should indicate the TCAD original mesh coordinate (*x*, *y* or *z*), and so on for the second (*y*) and third (*z*) array entry. For example, if one wants to have the TCAD *x*, *y* and *z* mesh coordinates mapped into the *y*, *z* and *x* coordinates of the new mesh, respectively, the configuration file should have `xyz = z x y`. If one wants to flip one of the coordinates, the minus symbol (`-`) can be used in front of one of the coordinates (such as `xyz = z x -y`).

The program can be used to convert 3D and 2D TCAD mesh files. Note that when converting 2D meshes, the *x* coordinate will be fixed to 1 and the interpolation will happen over the *yz* plane.
The keyword mesh_tree can be used as a switch to enable or disable the creation of a root file with the original TCAD mesh points stored as a `ROOT::TTree` for later, fast, inspection.

In addition, the `mesh_plotter` tool can be used, in order to visualize the new mesh interpolation results, from the installation folder as follows:
```bash
mesh_plotter -f <file_name> [<options>] [<arguments>]
```
The following command-line options are supported:
```
-f <file_name>         name of the interpolated file in APF or INIT format
-c <cut>               projection height index (default is mesh_pitch / 2)
-h                     display this help text
-l                     plot with logarithmic scale if set
-o <output_file_name>  name of the file to output (default is efield.png)
-p <plane>             plane to be plotted. xy, yz or zx (default is yz)
```

The list with options and defaults is displayed with the `-h` option.
In a 3D mesh, the plane to be plotted must be identified by using the option `-p` with argument *xy*, *yz* or *zx*, defaulting to *yz*.
The data to be plotted can be selected with the `-d` option, the arguments are *ex*, *ey*, *ez* for the vector components or the default value *n* for the norm of the electric field.
The number of mesh divisions in each dimension is automatically read from the `init`/`apf` file, by default the cut in the third dimension is done in the center but can be shifted using the `-c` option described above.

### Octree
J. Behley, V. Steinhage, A.B. Cremers. *Efficient Radius Neighbor Search in Three-dimensional Point Clouds*, Proc. of the IEEE International Conference on Robotics and Automation (ICRA), 2015 [@octree].

Copyright 2015 Jens Behley, University of Bonn.
This project is free software made available under the MIT License. For details see the LICENSE.md file.

[@octree]: http://jbehley.github.io/papers/behley2015icra.pdf
