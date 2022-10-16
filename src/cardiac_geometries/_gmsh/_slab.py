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
