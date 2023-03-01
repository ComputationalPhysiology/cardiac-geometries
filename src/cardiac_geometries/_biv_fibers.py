from pathlib import Path
from typing import Dict
from typing import Optional
from typing import Tuple
from typing import Union

import dolfin

from ._lv_ellipsoid_fibers import Microstructure


def create_biv_fibers(
    mesh: dolfin.Mesh,
    ffun: dolfin.MeshFunction,
    markers: Dict[str, Tuple[int, int]],
    fiber_space: str,
    alpha_endo=60,
    alpha_epi=-60,
    outdir: Optional[Union[str, Path]] = None,
) -> Microstructure:
    try:
        import ldrb
    except ImportError as e:
        msg = "Need ldrb to create fibers for BiV geometry - pip install ldrb"
        raise ImportError(msg) from e

    ldrb_markers = {
        "base": markers["BASE"][0],
        "lv": markers["ENDO_LV"][0],
        "rv": markers["ENDO_RV"][0],
        "epi": markers["EPI"][0],
    }

    f0, s0, n0 = ldrb.dolfin_ldrb(
        mesh=mesh,
        fiber_space=fiber_space,
        ffun=ffun,
        markers=ldrb_markers,
        alpha_endo_lv=alpha_endo,  # Fiber angle on the endocardium
        alpha_epi_lv=alpha_epi,  # Fiber angle on the epicardium
    )

    if outdir is not None:
        path = Path(outdir) / "microstructure.h5"
        with dolfin.HDF5File(mesh.mpi_comm(), path.as_posix(), "w") as h5file:
            h5file.write(f0, "f0")
            h5file.write(s0, "s0")
            h5file.write(n0, "n0")

    return Microstructure(f0=f0, s0=s0, n0=n0)
