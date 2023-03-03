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
