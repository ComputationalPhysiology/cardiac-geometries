from pathlib import Path
from typing import Dict
from typing import Optional
from typing import Tuple
from typing import Union

import dolfin
import numpy as np
import ufl

from .._import_checks import has_ldrb
from ._utils import Microstructure
from ._utils import save_microstructure, facet_function_from_heart_mesh


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


def setup_assigner(vs, index):
    # Set-up separate potential function for post processing
    VS0 = vs.function_space().sub(index)
    V = VS0.collapse()
    v = dolfin.Function(V)
    # Set-up object to optimize assignment from a function to subfunction
    v_assigner = dolfin.FunctionAssigner(VS0, V)
    v_assigner.assign(vs.sub(index), v)
    return v, v_assigner


def sub_function(vs, index):
    sub_vs = dolfin.Function(dolfin.cpp.function.Function(vs._cpp_object, index))
    ufl.Coefficient.__init__(
        sub_vs,
        sub_vs.function_space().ufl_function_space(),
        count=sub_vs._cpp_object.id(),
    )
    return sub_vs


def create_biv_in_torso_fibers(
    mesh: dolfin.Mesh,
    ffun: dolfin.MeshFunction,
    cfun: dolfin.MeshFunction,
    markers: Dict[str, Tuple[int, int]],
    fiber_space: str = "P_1",
    alpha_endo=60,
    alpha_epi=-60,
    outdir: Optional[Union[str, Path]] = None,
) -> Microstructure:
    if not has_ldrb():
        msg = "Need ldrb to create fibers for BiV geometry - pip install ldrb"
        raise ImportError(msg)

    # Extract heart
    heart = dolfin.MeshView.create(cfun, markers["HEART"][0])

    # And vertex map
    vmap = np.array(heart.topology().mapping()[mesh.id()].vertex_map(), dtype="uint64")

    # Extract facet function
    new_ffun = facet_function_from_heart_mesh(ffun, heart)

    # Only works for P 1 because we need vertex map
    assert fiber_space == "P_1"

    import ldrb

    ldrb_markers = {
        "base": markers["BASE"][0],
        "lv": markers["ENDO_LV"][0],
        "rv": markers["ENDO_RV"][0],
        "epi": markers["EPI"][0],
    }

    f0, s0, n0 = ldrb.dolfin_ldrb(
        mesh=heart,
        fiber_space=fiber_space,
        ffun=new_ffun,
        markers=ldrb_markers,
        alpha_endo_lv=alpha_endo,  # Fiber angle on the endocardium
        alpha_epi_lv=alpha_epi,  # Fiber angle on the epicardium
    )

    Vh = f0.function_space()
    V = dolfin.FunctionSpace(mesh, Vh.ufl_element())

    f0_torso = dolfin.Function(V)
    s0_torso = dolfin.Function(V)
    n0_torso = dolfin.Function(V)
    # We need to take dimension by dimension because the vertex
    # map only works on coordinates (not on vertices)
    for dim in range(3):
        V_ = V.sub(dim)
        Vh_ = Vh.sub(dim)
        torso_dofs = np.array(V_.dofmap().dofs(), dtype=int)
        heart_dofs = np.array(Vh_.dofmap().dofs(), dtype=int)

        if dim == 0:
            f0_torso.vector()[torso_dofs] = 1.0
        elif dim == 1:
            s0_torso.vector()[torso_dofs] = 1.0
        else:
            n0_torso.vector()[torso_dofs] = 1.0

        torso_inds = torso_dofs[dolfin.vertex_to_dof_map(V_.collapse())[vmap]]
        heart_inds = heart_dofs[dolfin.vertex_to_dof_map(Vh_.collapse())]
        f0_torso.vector()[torso_inds] = f0.vector()[heart_inds]
        s0_torso.vector()[torso_inds] = s0.vector()[heart_inds]
        n0_torso.vector()[torso_inds] = n0.vector()[heart_inds]

    system = Microstructure(f0=f0_torso, s0=s0_torso, n0=n0_torso)

    save_microstructure(system=system, outdir=outdir, comm=mesh.mpi_comm())

    return system
