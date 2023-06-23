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


def test_create_lv_ellipsoid(tmp_path: Path):
    runner = CliRunner()
    res1 = runner.invoke(
        cli.create_lv_ellipsoid,
        ["--create-fibers", tmp_path.as_posix()],
    )
    assert res1.exit_code == 0
    outfile = tmp_path / "geo.h5"
    res2 = runner.invoke(cli.folder2h5, [tmp_path.as_posix(), "--outfile", outfile])
    assert res2.exit_code == 0
    assert outfile.is_file()
    geo = Geometry.from_file(outfile)
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
        [tmp_path.as_posix()],
    )
    assert res1.exit_code == 0
    outfile = tmp_path / "geo.h5"
    res2 = runner.invoke(cli.folder2h5, [tmp_path.as_posix(), "--outfile", outfile])
    assert res2.exit_code == 0
    assert outfile.is_file()
    geo = Geometry.from_file(outfile)
    assert geo.mesh is not None
