import json
from enum import Enum
from pathlib import Path
from typing import Any
from typing import Dict
from typing import NamedTuple
from typing import Optional
from typing import Tuple
from typing import Union

import dolfin
from dolfin import FiniteElement  # noqa: F401
from dolfin import tetrahedron  # noqa: F401
from dolfin import VectorElement  # noqa: F401

from ._dolfin_utils import read_meshfunction
from .viz import dict_to_h5
from .viz import h5_to_dict
from .viz import h5pyfile


class H5Paths(str, Enum):
    mesh = "/mesh"
    cfun = "/meshfunctions/cfun"
    ffun = "/meshfunctions/ffun"
    efun = "/meshfunctions/efun"
    vfun = "/meshfunctions/vfun"
    f0 = "/microstructure/f0"
    s0 = "/microstructure/s0"
    n0 = "/microstructure/n0"
    markers = "/markers"
    info = "/info"


class Paths(NamedTuple):
    outdir: Path

    @property
    def tetra(self) -> Path:
        return self.outdir / "mesh.xdmf"

    @property
    def triangle(self) -> Path:
        return self.outdir / "triangle_mesh.xdmf"

    @property
    def line(self) -> Path:
        return self.outdir / "line_mesh.xdmf"

    @property
    def vertex(self) -> Path:
        return self.outdir / "vertex_mesh.xdmf"

    @property
    def microstructure(self) -> Path:
        return self.outdir / "microstructure.h5"

    @property
    def info(self) -> Path:
        return self.outdir / "info.json"

    @property
    def markers(self) -> Path:
        return self.outdir / "markers.json"


class Microstructure(NamedTuple):
    f0: dolfin.Function
    s0: dolfin.Function
    n0: dolfin.Function


def load_microstructure(mesh, microstructure_path) -> Microstructure:
    # Get signature
    with h5pyfile(microstructure_path) as h5file:
        signature = h5file["f0"].attrs["signature"].decode()

    V = dolfin.FunctionSpace(mesh, eval(signature))
    f0 = dolfin.Function(V)
    s0 = dolfin.Function(V)
    n0 = dolfin.Function(V)

    with dolfin.HDF5File(
        mesh.mpi_comm(),
        microstructure_path.as_posix(),
        "r",
    ) as h5file:
        h5file.read(f0, "f0")
        h5file.read(s0, "s0")
        h5file.read(n0, "n0")

    return Microstructure(f0=f0, s0=s0, n0=n0)


class Geometry(NamedTuple):
    mesh: dolfin.Mesh
    cfun: Optional[dolfin.MeshFunction] = None
    ffun: Optional[dolfin.MeshFunction] = None
    efun: Optional[dolfin.MeshFunction] = None
    vfun: Optional[dolfin.MeshFunction] = None
    markers: Optional[Dict[str, Tuple[int, int]]] = None
    info: Optional[Dict[str, Any]] = None
    f0: Optional[dolfin.Function] = None
    s0: Optional[dolfin.Function] = None
    n0: Optional[dolfin.Function] = None

    def save(self, fname: Union[str, Path]) -> None:
        path = Path(fname).with_suffix(".h5")
        path.unlink(missing_ok=True)

        with dolfin.HDF5File(self.mesh.mpi_comm(), path.as_posix(), "w") as h5file:

            h5file.write(self.mesh, H5Paths.mesh.value)
            for obj, p in [
                (self.cfun, H5Paths.cfun.value),
                (self.ffun, H5Paths.ffun.value),
                (self.efun, H5Paths.efun.value),
                (self.vfun, H5Paths.vfun.value),
                (self.f0, H5Paths.f0.value),
                (self.s0, H5Paths.s0.value),
                (self.n0, H5Paths.n0.value),
            ]:
                if obj is not None:
                    h5file.write(obj, p)

        if self.info is not None:
            dict_to_h5(self.info, path, H5Paths.info.value)
        if self.markers is not None:
            dict_to_h5(self.markers, path, H5Paths.markers.value)

    @classmethod
    def from_file(cls, fname: Union[str, Path]):
        path = Path(fname)
        if not path.is_file():
            msg = f"File {path} does not exist"
            raise FileNotFoundError(msg)

        groups = {}
        signatures = {"f0": None, "s0": None, "n0": None}
        with h5pyfile(path, "r") as h5file:

            if "info" in h5file:
                info = h5_to_dict(h5file["info"])
            else:
                info = None

            if "markers" in h5file:
                markers = h5_to_dict(h5file["markers"])
            else:
                markers = None

            for p in H5Paths:
                groups[p.name] = p.value in h5file

            for name in ["f0", "s0", "n0"]:
                if groups[name]:
                    signatures[name] = (
                        h5file[H5Paths[name].value].attrs["signature"].decode()
                    )

        if not groups.pop("mesh"):
            msg = f"Missing 'mesh' in {fname}"
            raise RuntimeError(msg)

        microstructure: Dict[str, dolfin.Function] = {}

        mesh = dolfin.Mesh()
        with dolfin.HDF5File(mesh.mpi_comm(), path.as_posix(), "r") as h5file:
            h5file.read(mesh, H5Paths.mesh, True)
            meshfunctions = {
                "cfun": dolfin.MeshFunction("size_t", mesh, 3),
                "ffun": dolfin.MeshFunction("size_t", mesh, 2),
                "efun": dolfin.MeshFunction("size_t", mesh, 1),
                "vfun": dolfin.MeshFunction("size_t", mesh, 0),
            }
            for name, exists in groups.items():
                if not exists:
                    continue

                if name in meshfunctions:
                    h5file.read(meshfunctions[name], H5Paths[name].value)
                    continue

                if name in signatures:
                    signature = signatures[name]
                    if signature is None:
                        continue
                    V = dolfin.FunctionSpace(mesh, eval(signature))
                    microstructure[name] = dolfin.Function(V)
                    h5file.read(microstructure[name], H5Paths[name].value)
                    continue

        return cls(
            mesh=mesh, info=info, markers=markers, **microstructure, **meshfunctions
        )

    @classmethod
    def from_folder(cls, folder):
        paths = Paths(folder)

        if not paths.tetra.is_file():
            msg = f"File {paths.tetra} not found"
            raise RuntimeError(msg)

        mesh = dolfin.Mesh()
        with dolfin.XDMFFile(paths.tetra.as_posix()) as h5file:
            h5file.read(mesh)

        cfun = dolfin.MeshFunction("size_t", mesh, 3)
        read_meshfunction(paths.tetra, cfun)

        ffun_val = dolfin.MeshValueCollection("size_t", mesh, 2)
        if paths.triangle.is_file():
            read_meshfunction(paths.triangle, ffun_val)
        ffun = dolfin.MeshFunction("size_t", mesh, ffun_val)

        efun_val = dolfin.MeshValueCollection("size_t", mesh, 1)
        if paths.line.is_file():
            read_meshfunction(paths.line, efun_val)
        efun = dolfin.MeshFunction("size_t", mesh, efun_val)

        vfun_val = dolfin.MeshValueCollection("size_t", mesh, 0)
        if paths.vertex.is_file():
            read_meshfunction(paths.vertex, vfun_val)
        vfun = dolfin.MeshFunction("size_t", mesh, vfun_val)

        if paths.microstructure.is_file():
            f0, s0, n0 = load_microstructure(
                mesh=mesh,
                microstructure_path=paths.microstructure,
            )
        else:
            f0 = s0 = n0 = None

        if paths.markers.is_file():
            markers = json.loads(paths.markers.read_text())
        else:
            markers: Dict[str, Tuple[int, int]] = {}

        if paths.info.is_file():
            info = json.loads(paths.info.read_text())
        else:
            info: Dict[str, Any] = {}

        return cls(
            mesh=mesh,
            cfun=cfun,
            ffun=ffun,
            efun=efun,
            vfun=vfun,
            f0=f0,
            s0=s0,
            n0=n0,
            markers=markers,
            info=info,
        )
