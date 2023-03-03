from pathlib import Path
from typing import Dict
from typing import Optional
from typing import Tuple
from typing import Union

import dolfin

from .._import_checks import has_ldrb
from ._utils import Microstructure
from ._utils import save_microstructure


def create_biv_fibers(
    mesh: dolfin.Mesh,
    ffun: dolfin.MeshFunction,
    markers: Dict[str, Tuple[int, int]],
    fiber_space: str,
    alpha_endo=60,
    alpha_epi=-60,
    outdir: Optional[Union[str, Path]] = None,
) -> Microstructure:

    if not has_ldrb():
        msg = "Need ldrb to create fibers for BiV geometry - pip install ldrb"
        raise ImportError(msg)

    import ldrb

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
    system = Microstructure(f0=f0, s0=s0, n0=n0)

    save_microstructure(system=system, outdir=outdir, comm=mesh.mpi_comm())

    return system
