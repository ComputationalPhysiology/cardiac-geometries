import gmsh
import numpy as np

from . import utils


def slab(mesh_name: str = "", lx=20.0, ly=7.0, lz=3.0, dx=1.0):

    path = utils.handle_mesh_name(mesh_name=mesh_name)
    # Initialize gmsh:
    gmsh.initialize()
    gmsh.model.add("Slab")
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", dx)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", dx)

    gmsh.model.occ.addBox(0, 0, 0, lx, ly, lz)
    gmsh.model.occ.synchronize()

    volumes = gmsh.model.getEntities(dim=3)

    myo_marker = 7
    gmsh.model.addPhysicalGroup(volumes[0][0], [volumes[0][1]], myo_marker)
    gmsh.model.setPhysicalName(volumes[0][0], myo_marker, "Myocardium")

    surfaces = gmsh.model.occ.getEntities(dim=2)
    x0, x1, y0, y1, z0, z1 = 1, 2, 3, 4, 5, 6

    for surface in surfaces:
        com = gmsh.model.occ.getCenterOfMass(surface[0], surface[1])

        if np.allclose(com, [0, ly / 2, lz / 2]):
            # X = 0
            gmsh.model.addPhysicalGroup(surface[0], [surface[1]], x0)
            gmsh.model.setPhysicalName(surface[0], x0, "X0")
        elif np.allclose(com, [lx, ly / 2, lz / 2]):
            # X = lx
            gmsh.model.addPhysicalGroup(surface[0], [surface[1]], x1)
            gmsh.model.setPhysicalName(surface[0], x1, "X1")
        elif np.allclose(com, [lx / 2, 0, lz / 2]):
            # Y = 0
            gmsh.model.addPhysicalGroup(surface[0], [surface[1]], y0)
            gmsh.model.setPhysicalName(surface[0], y0, "Y0")
        elif np.allclose(com, [lx / 2, ly, lz / 2]):
            # Y = ly
            gmsh.model.addPhysicalGroup(surface[0], [surface[1]], y1)
            gmsh.model.setPhysicalName(surface[0], y1, "Y1")
        elif np.allclose(com, [lx / 2, ly / 2, 0]):
            # Z = 0
            gmsh.model.addPhysicalGroup(surface[0], [surface[1]], z0)
            gmsh.model.setPhysicalName(surface[0], z0, "Z0")
        elif np.allclose(com, [lx / 2, ly / 2, lz]):
            # Z = lz
            gmsh.model.addPhysicalGroup(surface[0], [surface[1]], z1)
            gmsh.model.setPhysicalName(surface[0], z1, "Z1")

        else:
            print("Wtf!")

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)

    gmsh.write(path.as_posix())

    gmsh.finalize()
    return path


def slab_in_bath(
    mesh_name: str = "",
    lx=1.0,
    ly=0.01,
    lz=0.5,
    bx=0.0,
    by=0.0,
    bz=0.1,
    dx=0.001,
):
    """Create slab inside a bath

    Parameters
    ----------
    mesh_name : str, optional
        Name of o, by default ""
    lx : float, optional
        _description_, by default 1.0
    ly : float, optional
        _description_, by default 0.01
    lz : float, optional
        _description_, by default 0.5
    bx : float, optional
        _description_, by default 0.0
    by : float, optional
        _description_, by default 0.0
    bz : float, optional
        _description_, by default 0.1
    dx : float, optional
        _description_, by default 0.001

    Returns
    -------
    _type_
        _description_
    """

    path = utils.handle_mesh_name(mesh_name=mesh_name)
    # Initialize gmsh:
    gmsh.initialize()
    gmsh.model.add("SlabInBath")
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", dx)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", dx)

    slab_id = gmsh.model.occ.addBox(0, 0, 0, lx, ly, lz)
    bath_id = gmsh.model.occ.addBox(
        -bx,
        -by,
        -bz,
        lx + 2 * bx,
        ly + 2 * by,
        lz + 2 * bz,
    )

    vol = gmsh.model.occ.fragment([(3, slab_id)], [(3, bath_id)])

    # gmsh.model.occ.fragment([(3, vol[0][-1])], [(3, bath_id)])

    gmsh.model.occ.synchronize()

    volumes = gmsh.model.getEntities(dim=3)

    for vol in volumes:
        if np.allclose(
            gmsh.model.occ.getCenterOfMass(vol[0], vol[1]),
            [lx / 2, ly / 2, lz / 2],
        ):
            myo_id = vol[1]
            break

    myo_marker = 13
    gmsh.model.addPhysicalGroup(3, [myo_id], myo_marker)
    gmsh.model.setPhysicalName(3, myo_marker, "Myocardium")

    bath_marker = 14
    bath_ids = [vol[1] for vol in volumes if vol[1] != myo_id]
    gmsh.model.addPhysicalGroup(3, bath_ids, bath_marker)
    gmsh.model.setPhysicalName(3, bath_marker, "Bath")

    surfaces = gmsh.model.occ.getEntities(dim=2)

    x0, x1, y0, y1, z0, z1 = 1, 2, 3, 4, 5, 6
    bx0, bx1, by0, by1, bz0, bz1 = 7, 8, 9, 10, 11, 12

    for surface in surfaces:
        com = gmsh.model.occ.getCenterOfMass(surface[0], surface[1])

        if np.allclose(com, [0, ly / 2, lz / 2]):
            marker = x0
            name = "X0"
        elif np.allclose(com, [lx, ly / 2, lz / 2]):
            marker = x1
            name = "X1"
        elif np.allclose(com, [lx / 2, 0, lz / 2]):
            marker = y0
            name = "y0"
        elif np.allclose(com, [lx / 2, ly, lz / 2]):
            marker = y1
            name = "Y1"
        elif np.allclose(com, [lx / 2, ly / 2, 0]):
            marker = z0
            name = "Z0"
        elif np.allclose(com, [lx / 2, ly / 2, lz]):
            marker = z1
            name = "Z1"
        elif np.allclose(com, [-bx, ly / 2, lz / 2]):
            marker = bx0
            name = "BathX0"
        elif np.allclose(com, [lx + bx, ly / 2, lz / 2]):
            marker = bx1
            name = "BathX1"
        elif np.allclose(com, [lx / 2, -by, lz / 2]):
            marker = by0
            name = "BathY0"
        elif np.allclose(com, [lx / 2, ly + by, lz / 2]):
            marker = by1
            name = "BathY1"
        elif np.allclose(com, [lx / 2, ly / 2, -bz]):
            marker = bz0
            name = "BathZ0"
        elif np.allclose(com, [lx / 2, ly / 2, lz + bz]):
            marker = bz1
            name = "BathZ1"

        else:
            print("Wtf!")
            continue

        print(surface, com, marker, name)
        gmsh.model.addPhysicalGroup(surface[0], [surface[1]], marker)
        gmsh.model.setPhysicalName(surface[0], marker, name)

    gmsh.model.geo.synchronize()

    # gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.model.mesh.generate(3)

    gmsh.write(path.as_posix())

    gmsh.finalize()
    return path
