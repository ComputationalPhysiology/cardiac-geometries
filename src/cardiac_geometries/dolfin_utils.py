import typing
from pathlib import Path

import dolfin
import meshio

from . import calculus


class MarkerFunctions(typing.NamedTuple):
    vfun: dolfin.MeshFunction
    efun: dolfin.MeshFunction
    ffun: dolfin.MeshFunction
    cfun: dolfin.MeshFunction


class Geometry(typing.NamedTuple):
    mesh: dolfin.Mesh
    markers: typing.Dict[str, typing.Tuple[int, int]]
    marker_functions: MarkerFunctions


def create_mesh(mesh, cell_type):
    # From http://jsdokken.com/converted_files/tutorial_pygmsh.html
    cells = mesh.get_cells_type(cell_type)
    if cells.size == 0:
        return None

    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    out_mesh = meshio.Mesh(
        points=mesh.points,
        cells={cell_type: cells},
        cell_data={"name_to_read": [cell_data]},
    )
    return out_mesh


def read_meshfunction(fname, obj):
    try:
        with dolfin.XDMFFile(Path(fname).as_posix()) as f:
            f.read(obj, "name_to_read")
    except RuntimeError:
        pass


def gmsh2dolfin(
    msh_file: typing.Union[Path, str],
    outdir: typing.Optional[typing.Union[Path, str]] = None,
    unlink: bool = False,
) -> Geometry:

    msh = meshio.gmsh.read(msh_file)
    if outdir is None:
        outdir = Path(msh_file).absolute().parent
    else:
        outdir = Path(outdir)

    outdir.mkdir(exist_ok=True, parents=True)

    vertex_mesh = create_mesh(msh, "vertex")
    line_mesh = create_mesh(msh, "line")
    triangle_mesh = create_mesh(msh, "triangle")
    tetra_mesh = create_mesh(msh, "tetra")

    vertex_mesh_name = outdir / "vertex_mesh.xdmf"
    if vertex_mesh is not None:
        meshio.write(vertex_mesh_name, vertex_mesh)

    line_mesh_name = outdir / "line_mesh.xdmf"
    if line_mesh is not None:
        meshio.write(line_mesh_name, line_mesh)

    triangle_mesh_name = outdir / "triangle_mesh.xdmf"
    if triangle_mesh is not None:
        meshio.write(triangle_mesh_name, triangle_mesh)

    if tetra_mesh is None:
        raise RuntimeError("Unable to create mesh")

    tetra_mesh_name = outdir / "mesh.xdmf"
    meshio.write(
        tetra_mesh_name,
        tetra_mesh,
    )

    mesh = dolfin.Mesh()

    with dolfin.XDMFFile(tetra_mesh_name.as_posix()) as infile:
        infile.read(mesh)

    cfun = dolfin.MeshFunction("size_t", mesh, 3)
    read_meshfunction(tetra_mesh_name, cfun)
    if unlink:
        tetra_mesh_name.unlink(missing_ok=True)
        tetra_mesh_name.with_suffix(".h5").unlink(missing_ok=True)

    ffun_val = dolfin.MeshValueCollection("size_t", mesh, 2)
    read_meshfunction(triangle_mesh_name, ffun_val)
    ffun = dolfin.MeshFunction("size_t", mesh, ffun_val)
    for value in ffun_val.values():
        mesh.domains().set_marker(value, 2)
    ffun.array()[ffun.array() == max(ffun.array())] = 0
    if unlink:
        triangle_mesh_name.unlink(missing_ok=True)
        triangle_mesh_name.with_suffix(".h5").unlink(missing_ok=True)
    else:
        ffun_path = outdir / "ffun.xdmf"
        with dolfin.XDMFFile(ffun_path.as_posix()) as infile:
            infile.write(ffun)

    efun_val = dolfin.MeshValueCollection("size_t", mesh, 1)
    read_meshfunction(line_mesh_name, efun_val)
    efun = dolfin.MeshFunction("size_t", mesh, efun_val)
    efun.array()[efun.array() == max(efun.array())] = 0
    if unlink:
        line_mesh_name.unlink(missing_ok=True)
        line_mesh_name.with_suffix(".h5").unlink(missing_ok=True)

    vfun_val = dolfin.MeshValueCollection("size_t", mesh, 0)
    read_meshfunction(vertex_mesh_name, vfun_val)
    vfun = dolfin.MeshFunction("size_t", mesh, vfun_val)
    vfun.array()[vfun.array() == max(vfun.array())] = 0
    if unlink:
        vertex_mesh_name.unlink(missing_ok=True)
        vertex_mesh_name.with_suffix(".h5").unlink(missing_ok=True)

    markers = msh.field_data
    marker_functions = MarkerFunctions(vfun=vfun, efun=efun, ffun=ffun, cfun=cfun)

    geo = Geometry(
        mesh=mesh,
        markers=markers,
        marker_functions=marker_functions,
    )
    return geo


def mark_cell_function(fun, mesh, foc, regions):
    """
    Iterates over the mesh and stores the
    region number in a meshfunction
    """

    if foc is None:
        foc = calculus.estimate_focal_point(mesh)

    for cell in dolfin.cells(mesh):

        # Get coordinates to cell midpoint
        x = cell.midpoint().x()
        y = cell.midpoint().y()
        z = cell.midpoint().z()

        T = calculus.cartesian_to_prolate_ellipsoidal(x, y, z, foc)

        fun[cell] = calculus.strain_region_number(T, regions)

    return fun
