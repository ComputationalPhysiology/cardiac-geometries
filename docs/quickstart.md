# Getting started

If you have [installed](install.md) the packed with `pip` then you should have a command called `cardiac-geometries` available. You can see the help menu by executing the command `cardiac-geometries --help`

You can create a left ellipsoidal ventricular geometry by using the following command

```
cardiac-geometries create-lv-ellipsoid lv-mesh
```
or using the docker container
```
docker run --rm -v $PWD:/home/shared -w /home/shared -t ghcr.io/computationalphysiology/cardiac_geometries:latest cardiac-geometries create-lv-ellipsoid lv-mesh
```
This will create a new folder `lv-mesh` with the following content
```
lv-mesh
├── line_mesh.h5
├── line_mesh.xdmf
├── lv_ellipsoid.msh
├── markers.json
├── mesh.h5
├── mesh.xdmf
├── triangle_mesh.h5
├── triangle_mesh.xdmf
├── vertex_mesh.h5
└── vertex_mesh.xdmf
```
Here you will find both the `.msh` file that is created by `gmsh` and all the files created by `fenics`. You can also run the following command
```
cardiac-geometries create-lv-ellipsoid lv-mesh --create-fibers
```
and then one additional file called `microstructure.h5` will be saved to the output folder containing the fiber, sheet and sheet normal vector fields. To load these files using FEniCS you can use the following code snippet

```python
from pathlib import Path
import dolfin


def read_meshfunction(fname, obj):
    with dolfin.XDMFFile(Path(fname).as_posix()) as f:
        f.read(obj, "name_to_read")

outdir = Path("lv-mesh")
vertex_mesh_name = outdir / "vertex_mesh.xdmf"
line_mesh_name = outdir / "line_mesh.xdmf"
triangle_mesh_name = outdir / "triangle_mesh.xdmf"
tetra_mesh_name = outdir / "mesh.xdmf"
microstructure_path = outdir / "microstructure.h5"

mesh = dolfin.Mesh()

with dolfin.XDMFFile(tetra_mesh_name.as_posix()) as infile:
    infile.read(mesh)

cfun = dolfin.MeshFunction("size_t", mesh, 3)
read_meshfunction(tetra_mesh_name, cfun)

ffun_val = dolfin.MeshValueCollection("size_t", mesh, 2)
read_meshfunction(triangle_mesh_name, ffun_val)
ffun = dolfin.MeshFunction("size_t", mesh, ffun_val)

efun_val = dolfin.MeshValueCollection("size_t", mesh, 1)
read_meshfunction(line_mesh_name, efun_val)
efun = dolfin.MeshFunction("size_t", mesh, efun_val)

vfun_val = dolfin.MeshValueCollection("size_t", mesh, 0)
read_meshfunction(vertex_mesh_name, vfun_val)
vfun = dolfin.MeshFunction("size_t", mesh, vfun_val)

V = dolfin.VectorFunctionSpace(mesh, "P", 1)
f0 = dolfin.Function(V)
s0 = dolfin.Function(V)
n0 = dolfin.Function(V)

with dolfin.HDF5File(mesh.mpi_comm(), microstructure_path.as_posix(), "r") as h5file:
    h5file.read(f0, "f0")
    h5file.read(s0, "s0")
    h5file.read(n0, "n0")
```

Note that here the fibers as assumed to be in a first order Lagrange space ($\mathbb{P}_1$). If you want another function space, e.g fourth order quadrature space, you can use the flag `--fiber-space`, e.g

```
cardiac-geometries create-lv-ellipsoid lv-mesh --create-fibers --fiber-space=Quadrature_4
```
in which case the function space should be
```python
V = dolfin.VectorFunctionSpace(
    mesh,
    dolfin.FiniteElement(
        family="Quadrature",
        cell=mesh.ufl_cell(),
        degree=4,
        quad_scheme="default",
    )
)
```

## Creating the mesh from the Land Cardiac Mechanics Benchmark paper

> Land S, Gurev V, Arens S, Augustin CM, Baron L, Blake R, Bradley C, Castro S,
 Crozier A, Favino M, Fastl TE. Verification of cardiac mechanics software:
 benchmark problems and solutions for testing active and passive material
 behaviour. Proc. R. Soc. A. 2015 Dec 8;471(2184):20150641.

To create the geometry from the Land Cardiac Mechanics Benchmark paper you can use the following command

```
cardiac-geometries create-lv-ellipsoid benchmark-mesh --create-fibers --r-short-endo=7 --r-short-epi=10 --r-long-endo=17 --r-long-epi=20  --mu-apex-endo=0 --mu-apex-epi=0 --mu-base-endo=1.869328527979773 --mu-base-epi=1.8234765819369754 --psize-ref=3
```
Note that this mesh will have a flat base located at $z_{\text{base}} = 5$. The `mu-base-endo` and `mu-base-epi` is computed using the following formula

$$
\mu_{\text{base}}^{\text{endo}} = \mathrm{arccos} \left( \frac{z_{\text{base}}}{r\_ long \_ endo} \right)
$$

and

$$
\mu_{\text{base}}^{\text{epi}} = \mathrm{arccos} \left( \frac{z_{\text{base}}}{r\_ long \_ epi} \right)
$$

which means that if you want the base to be located at $z = 0$ you can set $\mu_{\text{base}}^{\text{endo}} = \mu_{\text{base}}^{\text{epi}} = \frac{\pi}{2} \approx 1.5707963267948966$, i.e

```
cardiac-geometries create-lv-ellipsoid mymesh --create-fibers --r-short-endo=7 --r-short-epi=10 --r-long-endo=17 --r-long-epi=20  --mu-apex-endo=0 --mu-apex-epi=0 --mu-base-endo=1.5707963267948966 --mu-base-epi=1.5707963267948966 --psize-ref=3
```

## Creating mesh for a slab

If it also possible to create a slab tissue using the command `create-slab` (note that this require `gmsh` to be installed with with [OpenCascade kernel](https://www.opencascade.com) which is installed in the docker image). In this case you can specify the length in the $x, y$ and $z$ direction through the arguments `lx`, `ly` and `lz` respectively. In this case the domain will be the box $[0, l_x] \times [0, l_y] \times [0, l_z]$ where is assumed that the apex is located at the place $z = 0$, the base at $z = l_z$, the endocardium at $y = 0$ and the epicardium at $y = l_y$. Also here you can generate fibers with the same rule based approach.
