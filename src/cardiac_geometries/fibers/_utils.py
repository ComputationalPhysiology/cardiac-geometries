from pathlib import Path
from typing import NamedTuple
from typing import Optional
from typing import Union

import dolfin
import mpi4py


class Microstructure(NamedTuple):
    f0: dolfin.Function
    s0: dolfin.Function
    n0: dolfin.Function


def save_microstructure(
    system: Microstructure,
    outdir: Optional[Union[Path, str]],
    comm: mpi4py.MPI.Intracomm,
) -> None:
    if outdir is None:
        return
    path = Path(outdir) / "microstructure.h5"
    with dolfin.HDF5File(comm, Path(path).as_posix(), "w") as h5file:
        h5file.write(system.f0, "f0")
        h5file.write(system.s0, "s0")
        h5file.write(system.n0, "n0")

    with dolfin.XDMFFile(
        comm,
        (Path(outdir) / "microstructure_viz.xdmf").as_posix(),
    ) as xdmf:
        xdmf.write_checkpoint(system.f0, "f0", 0, dolfin.XDMFFile.Encoding.HDF5, False)
        xdmf.write_checkpoint(system.s0, "s0", 0, dolfin.XDMFFile.Encoding.HDF5, True)
        xdmf.write_checkpoint(system.n0, "n0", 0, dolfin.XDMFFile.Encoding.HDF5, True)


def facet_function_from_heart_mesh(
    ffun: dolfin.MeshFunction,
    heart_mesh: dolfin.Mesh,
) -> dolfin.MeshFunction:
    try:
        from scipy.spatial import KDTree
    except ImportError as e:
        raise ImportError("Please install scipy - pip install scipy") from e

    assert ffun.mesh().id() in heart_mesh.topology().mapping()

    # Create a KDTree for the facet midpoints
    tree = KDTree([f.midpoint().array() for f in dolfin.facets(ffun.mesh())])

    # New facet function
    D = ffun.mesh().topology().dim()
    new_ffun = dolfin.MeshFunction("size_t", heart_mesh, D - 1, 0)

    for heart_cell in dolfin.cells(heart_mesh):
        for facet in dolfin.facets(heart_cell):
            global_index = tree.query(facet.midpoint().array())[1]
            new_ffun[facet] = ffun[global_index]

    return new_ffun
