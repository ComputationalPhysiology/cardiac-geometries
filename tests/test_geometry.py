import json
from typing import Any
from typing import Dict
from typing import NamedTuple

import dolfin
import pytest
from cardiac_geometries.geometry import Geometry
from cardiac_geometries.geometry import H5Path
from cardiac_geometries.geometry import load_scheme


class ExampleData(NamedTuple):
    info: Dict[str, Any]
    mesh: dolfin.Mesh
    ffun: dolfin.MeshFunction
    f0: dolfin.Function
    info2: Dict[str, Any]
    mesh2: dolfin.Mesh
    ffun2: dolfin.MeshFunction
    f02: dolfin.Function


@pytest.fixture(scope="session")
def example_data():
    info = {"param1": 13, "param2": 20}
    mesh = dolfin.UnitCubeMesh(2, 2, 2)
    ffun = dolfin.MeshFunction("size_t", mesh, 2)
    ffun.array()[0] = 42
    V = dolfin.VectorFunctionSpace(mesh, "CG", 1)
    f0 = dolfin.Function(V)
    f0.vector()[0] = 3

    info2 = {"param1": 12, "param2": 21}
    mesh2 = dolfin.UnitCubeMesh(3, 3, 3)
    ffun2 = dolfin.MeshFunction("size_t", mesh2, 2)
    ffun2.array()[0] = 41
    V = dolfin.VectorFunctionSpace(mesh2, "CG", 1)
    f02 = dolfin.Function(V)
    f02.vector()[0] = 2

    return ExampleData(
        info=info,
        mesh=mesh,
        ffun=ffun,
        f0=f0,
        info2=info2,
        mesh2=mesh2,
        ffun2=ffun2,
        f02=f02,
    )


def test_load_invalid_scheme(tmp_path):
    scheme = {
        "mesh": dict(h5group="/mesh", is_mesh=True, invalid_key=42),
        "ffun": dict(h5group="/ffun", is_meshfunction=True, dim=2),
        "f0": dict(h5group="/f0", is_function=True),
    }
    path = tmp_path / "scheme.json"
    path.write_text(json.dumps(scheme, indent=2))
    new_scheme = load_scheme(path)
    for name, d in scheme.items():
        for k, v in d.items():
            if k == "invalid_key":
                continue
            assert getattr(new_scheme[name], k) == v


def test_save_load_simple(tmp_path, example_data):
    scheme = {
        "info": H5Path(h5group="/info", is_dolfin=False),
        "mesh": H5Path(h5group="/mesh", is_mesh=True),
        "ffun": H5Path(h5group="/ffun", is_meshfunction=True, dim=2, mesh_key="mesh"),
        "f0": H5Path(h5group="/f0", is_function=True, mesh_key="mesh"),
    }
    geo = Geometry(**example_data._asdict(), scheme=scheme)
    path = tmp_path / "geo.h5"
    scheme_path = path.with_suffix(".json")
    geo.save(path, scheme_path=scheme_path)
    assert path.is_file()

    new_geo = Geometry.from_file(path, scheme_path=scheme_path)
    assert new_geo.scheme == geo.scheme
    assert new_geo.info == geo.info
    assert (new_geo.mesh.coordinates() == geo.mesh.coordinates()).all()
    assert (new_geo.ffun.array() == geo.ffun.array()).all()
    assert (new_geo.f0.vector().get_local() == geo.f0.vector().get_local()).all()


def test_save_load_multiple_meshes(tmp_path, example_data):
    scheme = {
        "info": H5Path(h5group="/group1/info", is_dolfin=False),
        "mesh": H5Path(h5group="/group1/mesh", is_mesh=True),
        "ffun": H5Path(
            h5group="/group1/ffun",
            is_meshfunction=True,
            dim=2,
            mesh_key="mesh",
        ),
        "f0": H5Path(h5group="/group1/f0", is_function=True, mesh_key="mesh"),
        "info2": H5Path(h5group="/group2/info", is_dolfin=False),
        "mesh2": H5Path(h5group="/group2/mesh", is_mesh=True),
        "ffun2": H5Path(
            h5group="/group2/ffun",
            is_meshfunction=True,
            dim=2,
            mesh_key="mesh2",
        ),
        "f02": H5Path(h5group="/group2/f0", is_function=True, mesh_key="mesh2"),
    }
    geo = Geometry(**example_data._asdict(), scheme=scheme)
    path = tmp_path / "geo.h5"
    scheme_path = path.with_suffix(".json")
    geo.save(path, scheme_path=scheme_path)
    assert path.is_file()

    new_geo = Geometry.from_file(path, scheme_path=scheme_path)
    assert new_geo.scheme == geo.scheme
    assert new_geo.info == geo.info
    assert (new_geo.mesh.coordinates() == geo.mesh.coordinates()).all()
    assert (new_geo.ffun.array() == geo.ffun.array()).all()
    assert (new_geo.f0.vector().get_local() == geo.f0.vector().get_local()).all()

    assert new_geo.info2 == geo.info2
    assert (new_geo.mesh2.coordinates() == geo.mesh2.coordinates()).all()
    assert (new_geo.ffun2.array() == geo.ffun2.array()).all()
    assert (new_geo.f02.vector().get_local() == geo.f02.vector().get_local()).all()