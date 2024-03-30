from pathlib import Path

from cardiac_geometries import cli
from cardiac_geometries.geometry import Geometry
from click.testing import CliRunner


def test_create_slab(tmp_path: Path):
    runner = CliRunner()
    res1 = runner.invoke(cli.create_slab, ["--create-fibers", tmp_path.as_posix()])
    assert res1.exit_code == 0
    outfile = tmp_path / "geo.h5"
    res2 = runner.invoke(cli.folder2h5, [tmp_path.as_posix(), "--outfile", outfile])
    assert res2.exit_code == 0
    assert outfile.is_file()
    geo = Geometry.from_file(outfile)
    assert geo.f0 is not None


def test_create_slab_in_bath(tmp_path: Path):
    runner = CliRunner()
    res1 = runner.invoke(cli.create_slab_in_bath, [tmp_path.as_posix()])
    assert res1.exit_code == 0
    outfile = tmp_path / "geo.h5"
    res2 = runner.invoke(cli.folder2h5, [tmp_path.as_posix(), "--outfile", outfile])
    assert res2.exit_code == 0
    assert outfile.is_file()
    geo = Geometry.from_file(outfile)

    assert "Bath" in geo.markers
    bath_marker = geo.markers["Bath"][0]
    assert bath_marker in geo.cfun.array()


def test_create_lv_ellipsoid(tmp_path: Path):
    runner = CliRunner()
    res1 = runner.invoke(
        cli.create_lv_ellipsoid,
        ["--create-fibers", tmp_path.as_posix(), "--aha"],
    )
    assert res1.exit_code == 0
    outfile = tmp_path / "geo.h5"
    res2 = runner.invoke(cli.folder2h5, [tmp_path.as_posix(), "--outfile", outfile])
    assert res2.exit_code == 0
    assert outfile.is_file()
    geo = Geometry.from_file(outfile)
    assert geo.mesh.topology().dim() == 3
    assert geo.f0 is not None

    # Make sure all aha regions have nonzero volume
    import dolfin

    dx = dolfin.dx(domain=geo.mesh, subdomain_data=geo.cfun)
    for k, v in geo.markers.items():
        if v[1] != 3:
            continue
        vol = dolfin.assemble(dolfin.Constant(1.0) * dx(v[0]))
        assert vol > 0, k


def test_create_lv_ellipsoid_2D(tmp_path: Path):
    runner = CliRunner()
    res1 = runner.invoke(
        cli.create_lv_ellipsoid,
        [tmp_path.as_posix(), "--axisymmetric", "--create-fibers"],
    )
    assert res1.exit_code == 0, res1.stderr
    outfile = tmp_path / "geo.h5"
    res2 = runner.invoke(cli.folder2h5, [tmp_path.as_posix(), "--outfile", outfile])
    assert res2.exit_code == 0
    assert outfile.is_file()
    geo = Geometry.from_file(outfile)
    assert geo.mesh.topology().dim() == 2
    assert geo.f0 is not None


def test_create_biv_ellipsoid(tmp_path: Path):
    runner = CliRunner()
    res1 = runner.invoke(
        cli.create_biv_ellipsoid,
        ["--create-fibers", tmp_path.as_posix()],
    )
    assert res1.exit_code == 0
    outfile = tmp_path / "geo.h5"
    res2 = runner.invoke(cli.folder2h5, [tmp_path.as_posix(), "--outfile", outfile])
    assert res2.exit_code == 0
    assert outfile.is_file()
    geo = Geometry.from_file(outfile)
    assert geo.f0 is not None


def test_create_biv_ellipsoid_torso(tmp_path: Path):
    runner = CliRunner()
    res1 = runner.invoke(
        cli.create_biv_ellipsoid_torso,
        ["--create-fibers", tmp_path.as_posix()],
    )
    assert res1.exit_code == 0
    outfile = tmp_path / "geo.h5"
    res2 = runner.invoke(cli.folder2h5, [tmp_path.as_posix(), "--outfile", outfile])
    assert res2.exit_code == 0
    assert outfile.is_file()
    geo = Geometry.from_file(outfile)
    assert geo.mesh is not None
