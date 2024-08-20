from textwrap import dedent
import json
import shutil
from pathlib import Path
from collections import Counter
import logging
import argparse
import subprocess
import meshio

STLFILE = "surface_mesh.stl"
MSHFILE = "surface_mesh.msh"
GEOFILE = "surface_mesh.geo"
COARSE_FOLDER = "coarse"
ORIGINAL_FOLDER = "original"


geofile = dedent(
    """
Merge "surface_mesh.stl";
//+


DefineConstant[
  // Angle between two triangles above which an edge is considered as sharp
  angle = {40, Min 20, Max 120, Step 1,
    Name "Parameters/Angle for surface detection"},
  // For complex geometries, patches can be too complex, too elongated or too
  // large to be parametrized; setting the following option will force the
  // creation of patches that are amenable to reparametrization:
  forceParametrizablePatches = {1, Choices{0,1},
    Name "Parameters/Create surfaces guaranteed to be parametrizable"},
  // For open surfaces include the boundary edges in the classification process:
  includeBoundary = 0,
  // Force curves to be split on given angle:
  curveAngle = 180
];
ClassifySurfaces{angle * Pi/180, includeBoundary, forceParametrizablePatches,
                 curveAngle * Pi / 180};

// Create a geometry for all the discrete curves and surfaces in the mesh, by
// computing a parametrization for each one
CreateGeometry;

// Create a volume as usual
Surface Loop(1) = Surface{:};
Volume(1) = {1};

Physical Surface("LV", 1) = {1};
Physical Surface("RV", 2) = {2};
Physical Surface("EPI", 3) = {3};
Physical Surface("MV", 4) = {4};
Physical Volume("Wall", 8) = {1};

// We specify element sizes imposed by a size field, just because we can :-)
funny = DefineNumber[0, Choices{0,1},
  Name "Parameters/Apply funny mesh size field?" ];

Field[1] = MathEval;
If(funny)
  Field[1].F = "2*Sin((x+y)/5) + 3";
Else
  Field[1].F = "4";
EndIf
Background Field = 1;


OptimizeMesh "Gmsh";
Mesh.OptimizeThreshold = 0.5;
Mesh.AngleToleranceFacetOverlap = 0.04;
Mesh.Algorithm3D = 1;
Coherence;
Mesh.MshFileVersion = 2.2;
"""
)

here = Path(__file__).parent


logger = logging.getLogger(__name__)


def create_facet_function(mesh, u):
    logger.debug("Creating facet function")
    import dolfin

    tdim = mesh.topology().dim()
    mesh.init(tdim - 1, tdim)
    ffun = dolfin.MeshFunction("size_t", mesh, 2)
    for facet in dolfin.facets(mesh):
        if facet.exterior():
            midpoint = u(facet.midpoint())
            round_midpoint = int(round(midpoint, 1))
            corners = []
            neighbors = []
            for v in dolfin.vertices(facet):
                corners.append(u(v.point()))
                for f in dolfin.facets(v):
                    if f.exterior():
                        neighbors.append(round(u(f.midpoint())))
            count = Counter(neighbors).most_common()
            possible_values = [count[0][0]]
            if len(count) > 1:
                possible_values.append(count[1][0])

            if round_midpoint not in possible_values:
                round_midpoint = count[0][0]
            ffun[facet] = round_midpoint
        else:
            ffun[facet] = 0
    return ffun


def convert_surface_mesh_to_dolfin(vtp: Path, output: Path) -> None:
    """Convert a vtp file to a dolfin xdmf file using the surface mesh
    representation.

    Parameters
    ----------
    vtp : Path
        Path to the vtp file
    output : Path
        Path to the output folder
    """
    logger.info(f"Converting {vtp} to dolfin xdmf")
    import dolfin

    m = meshio.read(vtp)
    logger.debug(f"Writing {output / ORIGINAL_FOLDER / 'surface_mesh.xdmf'}")
    meshio.write(output / ORIGINAL_FOLDER / "surface_mesh.xdmf", m)
    mesh = dolfin.Mesh()
    with dolfin.XDMFFile((output / ORIGINAL_FOLDER / "surface_mesh.xdmf").as_posix()) as xdmf:
        xdmf.read(mesh)

    V = dolfin.FunctionSpace(mesh, "CG", 1)
    arr = m.point_data["class"][dolfin.dof_to_vertex_map(V)]
    u = dolfin.Function(V)
    u.set_allow_extrapolation(True)
    u.vector().set_local(arr)
    logger.debug(f"Writing {output / ORIGINAL_FOLDER / 'surface_data.xdmf'}")
    with dolfin.XDMFFile((output / ORIGINAL_FOLDER / "surface_data.xdmf").as_posix()) as xdmf:
        xdmf.write_checkpoint(u, "u", 0, dolfin.XDMFFile.Encoding.HDF5, False)

    import cardiac_geometries

    msh_path = output / MSHFILE
    if not msh_path.is_file():
        raise FileNotFoundError(f"File not found: {msh_path}")

    geo = cardiac_geometries.dolfin_utils.gmsh2dolfin(msh_path, output / COARSE_FOLDER)
    ffun = create_facet_function(geo.mesh, u)
    logger.debug(f"Writing {output / COARSE_FOLDER / 'ffun.xdmf'}")
    with dolfin.XDMFFile((output / COARSE_FOLDER / "ffun.xdmf").as_posix()) as xdmf:
        xdmf.write(ffun)

    (output / COARSE_FOLDER / "markers.json").write_text(
        json.dumps(
            {
                "BASE": [1, 2],
                "EPI": [2, 2],
                "ENDO_LV": [3, 2],
                "ENDO_RV": [4, 2],
            }
        )
    )


def convert_volume_mesh_to_dolfin(vtu: Path, output: Path) -> None:
    """Convert a vtu file to a dolfin xdmf file using the volume mesh
    representation.

    Parameters
    ----------
    vtu : Path
        Path to the vtu file
    output : Path
        Path to the output folder
    """
    logger.info(f"Converting {vtu} to dolfin xdmf")
    import dolfin

    m = meshio.read(vtu)
    logger.debug(f"Writing {output / ORIGINAL_FOLDER / 'mesh.xdmf'}")
    meshio.write(output / ORIGINAL_FOLDER / "mesh.xdmf", m)
    mesh = dolfin.Mesh()
    with dolfin.XDMFFile((output / ORIGINAL_FOLDER / "mesh.xdmf").as_posix()) as xdmf:
        xdmf.read(mesh)

    W = dolfin.FunctionSpace(mesh, "CG", 1)
    w = dolfin.Function(W)
    w.set_allow_extrapolation(True)

    coarsen_mesh = dolfin.Mesh()
    with dolfin.XDMFFile((output / COARSE_FOLDER / "mesh.xdmf").as_posix()) as xdmf:
        xdmf.read(coarsen_mesh)

    W_coarse = dolfin.FunctionSpace(coarsen_mesh, "CG", 1)
    w_coarse = dolfin.Function(W_coarse)

    logger.debug(f"Writing {output / ORIGINAL_FOLDER / 'volume_data.xdmf'}")
    logger.debug(f"Writing {output / COARSE_FOLDER / 'volume_data.xdmf'}")
    xdmf_orig = dolfin.XDMFFile((output / ORIGINAL_FOLDER / "volume_data.xdmf").as_posix())
    xdmf_coarse = dolfin.XDMFFile((output / COARSE_FOLDER / "volume_data.xdmf").as_posix())
    for k, v in m.point_data.items():
        w.vector().set_local(v[dolfin.dof_to_vertex_map(W)])
        xdmf_orig.write_checkpoint(w, k, 0, dolfin.XDMFFile.Encoding.HDF5, True)
        logger.debug(f"Interpolating {k} to coarse mesh")
        w_coarse.interpolate(w)
        xdmf_coarse.write_checkpoint(w_coarse, k, 0, dolfin.XDMFFile.Encoding.HDF5, True)

    xdmf_coarse.close()
    xdmf_orig.close()


def create_fine_facet_function(vtp: Path, output: Path) -> None:
    logger.info("Creating fine facet function")
    import dolfin

    m = meshio.read(vtp)

    mesh = dolfin.Mesh()
    with dolfin.XDMFFile((output / ORIGINAL_FOLDER / "surface_mesh.xdmf").as_posix()) as xdmf:
        xdmf.read(mesh)

    vmesh = dolfin.Mesh()
    with dolfin.XDMFFile((output / ORIGINAL_FOLDER / "mesh.xdmf").as_posix()) as xdmf:
        xdmf.read(vmesh)

    V = dolfin.FunctionSpace(mesh, "CG", 1)
    arr = m.point_data["class"][dolfin.dof_to_vertex_map(V)]
    u = dolfin.Function(V)
    u.set_allow_extrapolation(True)
    u.vector().set_local(arr)

    ffun = create_facet_function(vmesh, u)
    logger.debug(f"Writing {output / ORIGINAL_FOLDER / 'ffun.xdmf'}")
    with dolfin.XDMFFile((output / ORIGINAL_FOLDER / "ffun.xdmf").as_posix()) as xdmf:
        xdmf.write(ffun)

    (output / ORIGINAL_FOLDER / "markers.json").write_text(
        json.dumps(
            {
                "BASE": [1, 2],
                "EPI": [2, 2],
                "ENDO_LV": [3, 2],
                "ENDO_RV": [4, 2],
            }
        )
    )


def coarsen_mesh(vtp: Path, output: Path) -> None:
    """Convert a vtp file to a gmsh mesh file using the surface mesh
    representation. The surface mesh is coarsened using the gmsh
    algorithm.

    Parameters
    ----------
    vtp : Path
        Path to the vtp file
    output : Path
        Path to the output folder
    """
    logger.info(f"Converting {vtp} to gmsh mesh")
    m = meshio.read(vtp)
    logger.debug(f"Writing {output / STLFILE}")
    meshio.write(output / STLFILE, m)
    logger.debug(f"Writing {output / GEOFILE}")
    (output / GEOFILE).write_text(geofile)
    logger.debug(f"Coarsening using gmsh. Converting {STLFILE} to {MSHFILE}")
    subprocess.run(
        [
            "gmsh",
            GEOFILE,
            "-3",
            "-o",
            MSHFILE,
        ],
        cwd=output,
    )
    logger.debug("Finished running gmsh")


def create_fibers_ldrb(output: Path, folder: str) -> None:
    import dolfin

    try:
        import ldrb
    except ImportError:
        logger.error("ldrb library not found not found. Install with 'pip install ldrb'")
        return

    logger.info(f"Creating fibers for mesh in {folder}")

    # Load mesh
    mesh = dolfin.Mesh()
    with dolfin.XDMFFile((output / folder / "mesh.xdmf").as_posix()) as xdmf:
        xdmf.read(mesh)

    ffun = dolfin.MeshFunction("size_t", mesh, 2)
    with dolfin.XDMFFile((output / folder / "ffun.xdmf").as_posix()) as xdmf:
        xdmf.read(ffun)

    markers = json.loads((output / folder / "markers.json").read_text())

    # Create fibers
    ldrb_markers = {
        "base": markers["BASE"][0],
        "lv": markers["ENDO_LV"][0],
        "rv": markers["ENDO_RV"][0],
        "epi": markers["EPI"][0],
    }

    f0, s0, n0 = ldrb.dolfin_ldrb(
        mesh=mesh,
        fiber_space="CG_1",
        ffun=ffun,
        markers=ldrb_markers,
        alpha_endo_lv=60,
        alpha_epi_lv=-60,
        beta_endo_lv=0,
        beta_epi_lv=0,
    )

    # Save fibers
    with dolfin.XDMFFile((output / folder / "microstructure.xdmf").as_posix()) as xdmf:
        xdmf.write_checkpoint(f0, "f0", 0, dolfin.XDMFFile.Encoding.HDF5, True)
        xdmf.write_checkpoint(s0, "s0", 0, dolfin.XDMFFile.Encoding.HDF5, True)
        xdmf.write_checkpoint(n0, "n0", 0, dolfin.XDMFFile.Encoding.HDF5, True)


def get_parser():
    parser = argparse.ArgumentParser(description="Convert vtp and vtu file to dolfin xdmf")
    parser.add_argument("vtp", type=Path, help="Path to vtp file")
    parser.add_argument("vtu", type=Path, help="Path to vtu file")
    parser.add_argument("output", type=Path, default=None, help="Path to output folder")
    parser.add_argument("-f", "--force", action="store_true", help="Force overwrite")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    return parser


def check_args(vtp: Path, vtu: Path, output: Path) -> None:
    if not vtp.is_file():
        raise FileNotFoundError(f"File not found: {vtp}")
    if not vtu.is_file():
        raise FileNotFoundError(f"File not found: {vtu}")

    if not output.is_dir():
        logger.info(f"Creating output directory: {output}")
        output.mkdir(parents=True, exist_ok=True)

    logger.debug(f"Creating subdirectories in {output / COARSE_FOLDER}")
    (output / COARSE_FOLDER).mkdir(exist_ok=True)
    logger.debug(f"Creating subdirectories in {output / ORIGINAL_FOLDER}")
    (output / ORIGINAL_FOLDER).mkdir(exist_ok=True)


def main(
    vtp: Path,
    vtu: Path,
    output: Path,
    force: bool = False,
    verbose: bool = False,
    copy_original: bool = False,
    create_fibers: bool = False,
) -> int:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level)
    logging.getLogger("h5py").setLevel(logging.INFO)
    logging.getLogger("numba").setLevel(logging.INFO)
    logging.getLogger("FFC").setLevel(logging.WARNING)
    logging.getLogger("UFL_LEGACY").setLevel(logging.WARNING)

    check_args(vtp, vtu, output)

    if not (output / MSHFILE).is_file() or force:
        coarsen_mesh(vtp, output)

    if not (output / ORIGINAL_FOLDER / "surface_mesh.xdmf").is_file() or force:
        convert_surface_mesh_to_dolfin(vtp, output)

    if not (output / ORIGINAL_FOLDER / "mesh.xdmf").is_file() or force:
        convert_volume_mesh_to_dolfin(vtu, output)

    if not (output / ORIGINAL_FOLDER / "ffun.xdmf").is_file() or force:
        create_fine_facet_function(vtp, output)

    if create_fibers:
        if not (output / ORIGINAL_FOLDER / "microstructure.xdmf").is_file() or force:
            create_fibers_ldrb(output, ORIGINAL_FOLDER)

        if not (output / COARSE_FOLDER / "microstructure.xdmf").is_file() or force:
            create_fibers_ldrb(output, COARSE_FOLDER)

    if copy_original:
        logger.info("Copying original files")
        shutil.copy(vtp, output / "surface_mesh.vtp")
        shutil.copy(vtu, output / "mesh.vtu")

    return 0


if __name__ == "__main__":
    parser = get_parser()
    args = vars(parser.parse_args())
    raise SystemExit(main(**args))
