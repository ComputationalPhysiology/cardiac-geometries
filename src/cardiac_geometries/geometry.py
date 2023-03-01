import json
import warnings
from enum import auto
from enum import Enum
from pathlib import Path
from typing import Any
from typing import Dict
from typing import NamedTuple
from typing import Optional
from typing import Sequence
from typing import Tuple
from typing import Union

import dolfin
from dolfin import FiniteElement  # noqa: F401
from dolfin import tetrahedron  # noqa: F401
from dolfin import VectorElement  # noqa: F401

from .viz import dict_to_h5
from .viz import h5_to_dict
from .viz import h5pyfile


class MeshTypes(Enum):
    slab = auto()
    lv_ellipsoid = auto()
    biv_ellipsoid = auto()


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


class FileNames(str, Enum):
    mesh = "mesh.xdmf"
    cfun = "mesh.xdmf:name_to_read"
    ffun = "triangle_mesh.xdmf:name_to_read"
    efun = "line_mesh.xdmf:name_to_read"
    vfun = "vertex_mesh.xdmf:name_to_read"
    f0 = "microstructure.h5:f0"
    s0 = "microstructure.h5:s0"
    n0 = "microstructure.h5:n0"
    markers = "markers.json"
    info = "info.json"


class Microstructure(NamedTuple):
    f0: dolfin.Function
    s0: dolfin.Function
    n0: dolfin.Function


def load_microstructure(mesh, microstructure_path) -> Microstructure:
    # Get signature
    with h5pyfile(microstructure_path) as h5file:
        signature = h5file["f0"].attrs["signature"].decode()

    sig = eval(signature)
    sig._quad_scheme = "default"
    V = dolfin.FunctionSpace(mesh, sig)
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


def read_signature(fname, group):
    try:
        with h5pyfile(fname) as h5file:
            signature = h5file[group].attrs["signature"].decode()
    except Exception:
        return None
    return signature


class H5Path(NamedTuple):
    h5group: str = ""
    fname: str = ""
    mesh_key: str = ""
    is_dolfin: bool = True
    is_mesh: bool = False
    is_meshfunction: bool = False
    is_function: bool = False
    dim: int = -1

    def to_dict(self):
        return self._asdict()


def load_schema(path: Path) -> Optional[Dict[str, H5Path]]:
    if not path.is_file():
        return None

    data = json.loads(Path(path).read_text())

    # Remove invalid keys
    def filter_values(d: Dict[str, Any]):
        return {k: v for k, v in d.items() if k in H5Path._fields}

    return {k: H5Path(**filter_values(v)) for k, v in data.items()}


def extract_mesh_keys(d: Dict[str, H5Path]) -> Sequence[Tuple[str, str]]:
    required_mesh_keys = []
    for name, p in d.items():
        if p.is_mesh:
            continue
        if not p.is_dolfin:
            continue
        required_mesh_keys.append((name, p.mesh_key))
    return list(required_mesh_keys)


def read_xdmf(fname, obj, group=None):
    try:
        with dolfin.XDMFFile(Path(fname).as_posix()) as f:
            if group is None:
                f.read(obj)
            else:
                f.read(obj, group)
    except RuntimeError as e:
        msg = f"Unable to read {fname} and for group {group}. Got:\n" + e.args[0]
        warnings.warn(msg, category=UserWarning, stacklevel=3)


def read_h5(fname, comm, obj, group=None):
    try:
        with dolfin.HDF5File(comm, Path(fname).as_posix(), "r") as f:
            if group is None:
                f.read(obj)
            else:
                f.read(obj, group)
    except RuntimeError as e:
        msg = f"Unable to read {fname} and for group {group}. Got:\n" + e.args[0]
        warnings.warn(msg, category=UserWarning, stacklevel=3)


def read(
    fname: Path,
    obj: Any,
    group: Optional[str],
    current_mesh: dolfin.Mesh,
) -> None:
    if fname.suffix == ".xdmf":
        read_xdmf(fname, obj, group)
    elif fname.suffix == ".h5":
        read_h5(fname, current_mesh.mpi_comm(), obj, group)
    else:
        raise RuntimeError(f"Unknown file format for {fname}")


def extract_fname_group(fname: str, folder=".") -> Tuple[Path, Optional[str]]:
    fg = fname.split(":")
    if len(fg) == 1:
        return Path(folder) / fg[0], None
    if len(fg) > 2:
        print(f"Warning: Invalid fname {fname}")
    return Path(folder) / fg[0], fg[1]


def dump_schema(path: Union[Path, str], schema: Dict[str, H5Path]) -> None:
    Path(path).write_text(
        json.dumps({k: v._asdict() for k, v in schema.items()}, indent=2),
    )


class Geometry:
    def __init__(self, **kwargs):

        self._fields = ["schema"]
        schema = kwargs.pop("schema")
        if schema is None:
            schema = type(self).default_schema()

        self.schema = {}
        missing_schema_entries = []
        for k, v in kwargs.items():
            s = schema.get(k)
            if s is None:
                missing_schema_entries.append(k)
                continue
            self._fields.append(k)
            setattr(self, k, v)
            self.schema[k] = s

        if len(missing_schema_entries) > 0:
            msg = (
                f"Missing schema entry for keys {missing_schema_entries!r}. "
                "Objects will not be set as geometry attributes"
            )
            warnings.warn(UserWarning(msg), stacklevel=2)

    def __repr__(self) -> str:
        fields = ", ".join(self._fields)
        return f"{type(self).__name__}({fields})"

    @staticmethod
    def default_schema() -> Dict[str, H5Path]:
        return {
            "mesh": H5Path(
                h5group=H5Paths.mesh.value,
                is_mesh=True,
                fname=FileNames.mesh.value,
            ),
            "cfun": H5Path(
                h5group=H5Paths.cfun.value,
                is_meshfunction=True,
                dim=3,
                mesh_key="mesh",
                fname=FileNames.cfun.value,
            ),
            "ffun": H5Path(
                h5group=H5Paths.ffun.value,
                is_meshfunction=True,
                dim=2,
                mesh_key="mesh",
                fname=FileNames.ffun.value,
            ),
            "efun": H5Path(
                h5group=H5Paths.efun.value,
                is_meshfunction=True,
                dim=1,
                mesh_key="mesh",
                fname=FileNames.efun.value,
            ),
            "vfun": H5Path(
                h5group=H5Paths.vfun.value,
                is_meshfunction=True,
                dim=0,
                mesh_key="mesh",
                fname=FileNames.vfun.value,
            ),
            "f0": H5Path(
                h5group=H5Paths.f0.value,
                is_function=True,
                mesh_key="mesh",
                fname=FileNames.f0.value,
            ),
            "s0": H5Path(
                h5group=H5Paths.s0.value,
                is_function=True,
                mesh_key="mesh",
                fname=FileNames.s0.value,
            ),
            "n0": H5Path(
                h5group=H5Paths.n0.value,
                is_function=True,
                mesh_key="mesh",
                fname=FileNames.n0.value,
            ),
            "info": H5Path(
                h5group=H5Paths.info.value,
                is_dolfin=False,
                fname=FileNames.info.value,
            ),
            "markers": H5Path(
                h5group=H5Paths.markers.value,
                is_dolfin=False,
                fname=FileNames.markers.value,
            ),
        }

    def save(
        self,
        fname: Union[str, Path],
        schema_path: Optional[Union[str, Path]] = None,
        unlink: bool = True,
    ) -> None:
        path = Path(fname).with_suffix(".h5")
        file_mode = "w"
        if path.is_file():
            if unlink:
                path.unlink()
            else:
                file_mode = "a"

        with dolfin.HDF5File(
            dolfin.MPI.comm_world,
            path.as_posix(),
            file_mode,
        ) as h5file:
            for name, p in self.schema.items():
                obj = getattr(self, name, None)
                if obj is not None and p.is_dolfin:
                    if p.h5group == "":
                        raise RuntimeError("Cannot write object with empty path")
                    h5file.write(obj, p.h5group)

        for name, p in self.schema.items():
            obj = getattr(self, name, None)
            if obj is not None and not p.is_dolfin:
                if p.h5group == "":
                    raise RuntimeError("Cannot write object with empty path")
                dict_to_h5(obj, path, p.h5group)

        if schema_path is None:
            schema_path = path.with_suffix(".json")
        dump_schema(schema_path, schema=self.schema)

    @classmethod
    def from_file(
        cls,
        fname: Union[str, Path],
        schema_path: Optional[Union[str, Path]] = None,
        schema: Optional[Dict[str, H5Path]] = None,
    ):
        path = Path(fname)
        if not path.is_file():
            msg = f"File {path} does not exist"
            raise FileNotFoundError(msg)

        if schema_path is not None:
            schema = load_schema(Path(schema_path))

        if schema is None:
            schema = cls.default_schema()

        groups = {}
        data = {}
        signatures = {}
        with h5pyfile(path, "r") as h5file:

            for name, p in schema.items():
                if p.h5group == "":
                    continue
                groups[name] = p.h5group in h5file

                if not p.is_dolfin:
                    if groups[name]:
                        data[name] = h5_to_dict(h5file[p.h5group])
                    else:
                        data[name] = None

                if p.is_function and groups[name]:
                    signatures[name] = h5file[p.h5group].attrs["signature"].decode()

        required_mesh_keys = extract_mesh_keys(schema)
        for name, mesh_key in required_mesh_keys:
            if not groups.get(mesh_key, False):
                msg = f"Missing mesh key '{mesh_key}' for key {name} in {fname}"
                raise RuntimeError(msg)

        mesh = dolfin.Mesh()
        with dolfin.HDF5File(mesh.mpi_comm(), path.as_posix(), "r") as h5file:

            # Meshes
            for name, p in schema.items():
                if p.is_mesh and groups[name]:
                    data[name] = dolfin.Mesh()

                    h5file.read(data[name], p.h5group, True)

            # Meshfunctions
            for name, p in schema.items():
                if p.is_meshfunction and groups[name]:
                    current_mesh = data[p.mesh_key]
                    data[name] = dolfin.MeshFunction("size_t", current_mesh, p.dim)
                    h5file.read(data[name], p.h5group)
                    data[name].array()[data[name].array() == 2**64 - 1] = 0
                    continue

                if p.is_function and groups[name]:
                    current_mesh = data[p.mesh_key]
                    signature = signatures.get(name)
                    if signature is None:
                        continue
                    sig = eval(signature)
                    sig._quad_scheme = "default"
                    V = dolfin.FunctionSpace(current_mesh, sig)
                    data[name] = dolfin.Function(V)
                    h5file.read(data[name], p.h5group)
                    continue

        return cls(**data, schema=schema)

    @classmethod
    def from_folder(cls, folder, schema: Optional[Dict[str, H5Path]] = None):
        folder = Path(folder)

        if schema is None:
            schema = cls.default_schema()

        # Load mesh first
        data = {}

        for name, p in schema.items():
            if p.fname == "":
                continue
            if not p.is_mesh:
                continue
            fname, group = extract_fname_group(p.fname, folder)
            if not fname.is_file():
                continue
            data[name] = current_mesh = dolfin.Mesh()
            read(fname, data[name], group=group, current_mesh=current_mesh)

        # Now load rest of the dolfin stuff
        for name, p in schema.items():
            if p.fname == "":
                continue
            if p.is_meshfunction:
                fname, group = extract_fname_group(p.fname, folder)
                if not fname.is_file():
                    continue
                current_mesh = data[p.mesh_key]

                val = dolfin.MeshValueCollection("size_t", current_mesh, p.dim)
                read(fname, val, group=group, current_mesh=current_mesh)
                data[name] = dolfin.MeshFunction("size_t", current_mesh, val)
                data[name].array()[data[name].array() == 2**64 - 1] = 0
                continue

            if p.is_function:
                fname, group = extract_fname_group(p.fname, folder)
                if not fname.is_file():
                    continue
                signature = read_signature(fname, group)
                if signature is None:
                    continue
                current_mesh = data[p.mesh_key]
                sig = eval(signature)
                sig._quad_scheme = "default"
                V = dolfin.FunctionSpace(current_mesh, sig)
                data[name] = dolfin.Function(V)
                read(fname, data[name], group=group, current_mesh=current_mesh)
                continue

        # Finally read json
        for name, p in schema.items():
            if p.is_dolfin:
                continue
            fname = folder / p.fname
            if fname.suffix == ".json":
                data[name] = json.loads(fname.read_text())
            else:
                raise RuntimeError(f"Unknown file format for {fname}")

        return cls(**data, schema=schema)
