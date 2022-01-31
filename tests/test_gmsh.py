from pathlib import Path
from cardiac_geometries import gmsh


def test_lv_prolate_flat_base():
    mesh_path = Path("prolate_lv_ellipsoid_flat_base.msh")
    mesh_path.unlink(missing_ok=True)
    gmsh.prolate_lv_ellipsoid_flat_base(mesh_path)
    mesh_path.unlink(missing_ok=False)


def test_lv_prolate():
    mesh_path = Path("prolate_lv_ellipsoid.msh")
    mesh_path.unlink(missing_ok=True)
    gmsh.prolate_lv_ellipsoid(mesh_path)
    mesh_path.unlink(missing_ok=False)


def test_lv_flat_base():
    mesh_path = Path("lv_ellipsoid_flat_base.msh")
    mesh_path.unlink(missing_ok=True)
    gmsh.lv_ellipsoid_flat_base(mesh_path)
    mesh_path.unlink(missing_ok=False)


def test_lv_simple():
    mesh_path = Path("lv_ellipsoid.msh")
    mesh_path.unlink(missing_ok=True)
    gmsh.lv_ellipsoid(mesh_path, psize_ref=0.05)
    mesh_path.unlink(missing_ok=False)


def test_create_benchmark_geometry_land15():
    path = gmsh.create_benchmark_geometry_land15()
    path.unlink(missing_ok=False)
