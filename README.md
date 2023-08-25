# Cardiac Geometries

This is a library for creating idealized cardiac geometries in FEniCS using `gmsh`.
Currently you can create bi-ventricular (BiV) and left-ventricular (LV) ellipsoidal geometries as well as slab geometries. There is also support for creating bi-ventricular and slab geometries embedded in a torso / bath.


This package can also output analytic fiber orientations for LV and Slab and integrated with the [ldrb](https://github.com/finsberg/ldrb) algorithm BiV geometries.

## Install
User are encourage to use the the [provided docker image](https://github.com/ComputationalPhysiology/cardiac_geometries/pkgs/container/cardiac_geometries)
```
docker pull ghcr.io/computationalphysiology/cardiac_geometries:latest
```
which comes pre-installed with FEniCs and gmsh.

You can also install `cardiac-geometries` using `pip`
```
python3 -m pip install cardiac-geometries
```
but this requires FEniCS and gmsh to be installed in other ways. See more at https://computationalphysiology.github.io/cardiac_geometries/install.html

## Getting started
`cardiac-geometries` comes with a command-line interface, and to e.g create an LV mesh you can do
```
cardiac-geometries create-lv-ellipsoid lv-mesh
```
The same functionality can be accessed through the Python API, e.g
```python
import cardiac_geometries

geo = cardiac_geometries.mesh.create_lv_ellipsoid(outdir="lv-mesh")
```
See https://computationalphysiology.github.io/cardiac_geometries/quickstart.html for more info.


## Documentation
Please read the documentation at http://computationalphysiology.github.io/cardiac_geometries for more info.

## Contributing
See https://computationalphysiology.github.io/cardiac_geometries/CONTRIBUTING.html

## License
MIT
