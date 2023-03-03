from pathlib import Path
from typing import Dict
from typing import Optional
from typing import Tuple
from typing import Union

import dolfin
import numpy as np

from ._utils import Microstructure
from ._utils import save_microstructure


def laplace(
    mesh: dolfin.Mesh,
    ffun: dolfin.MeshFunction,
    markers: Dict[str, Tuple[int, int]],
    function_space: str = "P_1",
):

    endo_marker = markers["ENDO"][0]
    epi_marker = markers["EPI"][0]

    V = dolfin.FunctionSpace(mesh, "CG", 1)

    u = dolfin.TrialFunction(V)
    v = dolfin.TestFunction(V)
    a = dolfin.dot(dolfin.grad(u), dolfin.grad(v)) * dolfin.dx
    L = v * dolfin.Constant(0) * dolfin.dx

    t = dolfin.Function(V)

    endo_bc = dolfin.DirichletBC(V, 0, ffun, endo_marker, "topological")
    epi_bc = dolfin.DirichletBC(V, 1, ffun, epi_marker, "topological")
    bcs = [endo_bc, epi_bc]
    dolfin.solve(
        a == L,
        t,
        bcs,
    )
    if function_space != "P_1":
        family, degree = function_space.split("_")
        W = dolfin.FunctionSpace(
            mesh,
            dolfin.FiniteElement(
                family=family,
                cell=mesh.ufl_cell(),
                degree=int(degree),
                quad_scheme="default",
            ),
        )
        t = dolfin.interpolate(t, W)
    return t


def normalize(u):
    return u / np.linalg.norm(u, axis=0)


def compute_system(
    t_func: dolfin.Function,
    r_short_endo,
    r_short_epi,
    r_long_endo,
    r_long_epi,
    alpha_endo: float = -60,
    alpha_epi: float = 60,
    **kwargs,
):

    V = t_func.function_space()
    element = V.ufl_element()
    mesh = V.mesh()

    dof_coordinates = V.tabulate_dof_coordinates()

    alpha = lambda x: (alpha_endo + (alpha_epi - alpha_endo) * x) * (np.pi / 180)
    r_long = lambda x: r_long_endo + (r_long_epi - r_long_endo) * x
    r_short = lambda x: r_short_endo + (r_short_epi - r_short_endo) * x

    drl_dt = r_long_epi - r_long_endo
    drs_dt = r_short_epi - r_short_endo

    t = t_func.vector().get_local()

    rl = r_long(t)
    rs = r_short(t)
    al = alpha(t)

    x = dof_coordinates[:, 0]
    y = dof_coordinates[:, 1]
    z = dof_coordinates[:, 2]

    a = np.sqrt(y**2 + z**2) / rs
    b = x / rl
    mu = np.arctan2(a, b)
    theta = np.pi - np.arctan2(z, -y)
    theta[mu < 1e-7] = 0.0

    e_t = np.array(
        [
            drl_dt * np.cos(mu),
            drs_dt * np.sin(mu) * np.cos(theta),
            drs_dt * np.sin(mu) * np.sin(theta),
        ],
    )
    e_t = normalize(e_t)

    e_mu = np.array(
        [
            -rl * np.sin(mu),
            rs * np.cos(mu) * np.cos(theta),
            rs * np.cos(mu) * np.sin(theta),
        ],
    )
    e_mu = normalize(e_mu)

    e_theta = np.array(
        [
            np.zeros_like(t),
            -rs * np.sin(mu) * np.sin(theta),
            rs * np.sin(mu) * np.cos(theta),
        ],
    )
    e_theta = normalize(e_theta)

    f0 = np.sin(al) * e_mu + np.cos(al) * e_theta
    f0 = normalize(f0)

    s0 = np.cross(e_mu, e_theta, axis=0)
    s0 = normalize(s0)

    n0 = np.cross(f0, s0, axis=0)
    n0 = normalize(n0)

    Vv = dolfin.FunctionSpace(
        mesh,
        dolfin.VectorElement(
            family=element.family(),
            degree=element.degree(),
            cell=element.cell(),
            quad_scheme="default",
        ),
    )

    dim = Vv.mesh().geometry().dim()
    start, end = Vv.sub(0).dofmap().ownership_range()
    x_dofs = np.arange(0, end - start, dim)
    y_dofs = np.arange(1, end - start, dim)
    z_dofs = np.arange(2, end - start, dim)

    start, end = V.dofmap().ownership_range()
    scalar_dofs = [
        dof
        for dof in range(end - start)
        if V.dofmap().local_to_global_index(dof)
        not in V.dofmap().local_to_global_unowned()
    ]

    fiber = dolfin.Function(Vv)
    f = np.zeros_like(fiber.vector().get_local())

    f[x_dofs] = f0[0, scalar_dofs] / np.linalg.norm(f0, axis=0)
    f[y_dofs] = f0[1, scalar_dofs] / np.linalg.norm(f0, axis=0)
    f[z_dofs] = f0[2, scalar_dofs] / np.linalg.norm(f0, axis=0)

    fiber.vector().set_local(f)
    fiber.vector().apply("insert")
    fiber.rename("fiber", "fibers")

    sheet = dolfin.Function(Vv)
    s = np.zeros_like(f)
    s[x_dofs] = s0[0, scalar_dofs]
    s[y_dofs] = s0[1, scalar_dofs]
    s[z_dofs] = s0[2, scalar_dofs]
    sheet.vector().set_local(s)
    sheet.vector().apply("insert")
    sheet.rename("sheet", "fibers")

    sheet_normal = dolfin.Function(Vv)
    n = np.zeros_like(f)
    n[x_dofs] = n0[0, scalar_dofs]
    n[y_dofs] = n0[1, scalar_dofs]
    n[z_dofs] = n0[2, scalar_dofs]
    sheet_normal.vector().set_local(n)
    sheet_normal.vector().apply("insert")
    sheet_normal.rename("sheet_normal", "fibers")

    return Microstructure(f0=fiber, s0=sheet, n0=sheet_normal)


def create_microstructure(
    mesh,
    ffun,
    markers,
    r_short_endo,
    r_short_epi,
    r_long_endo,
    r_long_epi,
    alpha_endo,
    alpha_epi,
    function_space,
    outdir: Optional[Union[str, Path]] = None,
) -> Microstructure:

    t = laplace(mesh, ffun, markers, function_space=function_space)
    system = compute_system(
        t,
        r_short_endo=r_short_endo,
        r_short_epi=r_short_epi,
        r_long_endo=r_long_endo,
        r_long_epi=r_long_epi,
        alpha_endo=alpha_endo,
        alpha_epi=alpha_epi,
    )

    save_microstructure(system=system, outdir=outdir, comm=mesh.mpi_comm())

    return system
