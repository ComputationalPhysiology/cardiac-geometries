import math
from pathlib import Path

import gmsh

from . import utils


def biv_ellipsoid(
    mesh_name: str = "",
    center_lv=(0.0, 0.0, 0.0),
    radius_lv=1.0,
    radius_lv_endo=1.0,
    a_endo_lv: float = 2.5,
    b_endo_lv: float = 1.0,
    c_endo_lv: float = 1.0,
    a_epi_lv: float = 3.0,
    b_epi_lv: float = 1.5,
    c_epi_lv: float = 1.5,
    center_rv=(0.0, 0.5, 0.0),
    radius_rv=1.0,
    radius_rv_endo=1.0,
    a_endo_rv: float = 3.0,
    b_endo_rv: float = 1.5,
    c_endo_rv: float = 1.5,
    a_epi_rv: float = 4.0,
    b_epi_rv: float = 2.5,
    c_epi_rv: float = 2.0,
    angle1=-math.pi / 2,
    angle2=math.pi / 2,
    angle3=2 * math.pi,
    base_x: float = 0.0,
    char_length: float = 0.5,
) -> Path:
    path = utils.handle_mesh_name(mesh_name=mesh_name)
    gmsh.initialize()
    gmsh.model.add("biv")

    # Create a box for cutting the base
    a_epi = max(a_epi_lv, a_epi_rv)
    diam = -5 * a_epi  # Just make it sufficiently big
    box_id = gmsh.model.occ.addBox(base_x, a_epi, a_epi, diam, diam, diam)
    dim_tag_box = [(3, box_id)]

    # LV epicardium
    lv_id = gmsh.model.occ.addSphere(
        *center_lv, radius=radius_lv, angle1=angle1, angle2=angle2, angle3=angle3
    )
    dim_tag_lv = [(3, lv_id)]

    gmsh.model.occ.dilate(dim_tag_lv, *center_lv, a=a_epi_lv, b=b_epi_lv, c=c_epi_lv)

    # LV endocardium
    lv_endo_id = gmsh.model.occ.addSphere(
        *center_lv, radius=radius_lv_endo, angle1=angle1, angle2=angle2, angle3=angle3
    )
    dim_tag_lv_endo = [(3, lv_endo_id)]
    gmsh.model.occ.dilate(
        dim_tag_lv_endo, *center_lv, a=a_endo_lv, b=b_endo_lv, c=c_endo_lv
    )

    # RV epicardium
    rv_id = gmsh.model.occ.addSphere(
        *center_lv, radius=radius_rv, angle1=angle1, angle2=angle2, angle3=angle3
    )
    dim_tag_rv = [(3, rv_id)]

    gmsh.model.occ.dilate(dim_tag_rv, *center_lv, a=a_epi_rv, b=b_epi_rv, c=c_epi_rv)
    gmsh.model.occ.translate(dim_tag_rv, *center_rv)

    # RV endocardium
    rv_endo_id = gmsh.model.occ.addSphere(
        *center_rv, radius=radius_rv_endo, angle1=angle1, angle2=angle2, angle3=angle3
    )

    dim_tag_rv_endo = [(3, rv_endo_id)]
    gmsh.model.occ.dilate(
        dim_tag_rv_endo, *center_lv, a=a_endo_rv, b=b_endo_rv, c=c_endo_rv
    )
    gmsh.model.occ.translate(dim_tag_rv_endo, *center_rv)

    # Subtract LV epicardium from RV endo
    rv_endo, _ = gmsh.model.occ.cut(dim_tag_rv_endo, dim_tag_lv, removeTool=False)

    # Combine the two endocardial volumes
    endo, _ = gmsh.model.occ.fuse(rv_endo, dim_tag_lv_endo)
    # Combine the two epicardial volumes
    epi, _ = gmsh.model.occ.fuse(dim_tag_lv, dim_tag_rv)

    # Subtract the endocardial surfaces from the epicardial volumes
    ov, _ = gmsh.model.occ.cut(epi, endo)
    # And finally cut the base with the box
    ov2, _ = gmsh.model.occ.cut(ov, dim_tag_box)

    # Now we mark the different surfaces
    surfaces = gmsh.model.occ.getEntities(dim=2)

    gmsh.model.occ.synchronize()

    # Use the GUI to find out which surfaces are which
    # Alternatively we could do this programmatically using
    # the following functions, but lets assume that the order
    # remains the same.
    # for surf in surfaces:
    #     print(gmsh.model.occ.getCenterOfMass(surf[0], surf[1]))
    #     print(gmsh.model.occ.getBoundingBox(surf[0], surf[1]))
    #     print()

    base_marker = 1
    endo_lv_marker = 2
    endo_rv_marker = 3
    epi_marker = 4

    gmsh.model.addPhysicalGroup(
        dim=surfaces[1][0],
        tags=[surfaces[1][1]],
        tag=base_marker,
        name="BASE",
    )
    gmsh.model.addPhysicalGroup(
        dim=surfaces[6][0],
        tags=[surfaces[6][1], surfaces[7][1]],
        tag=endo_lv_marker,
        name="ENDO_LV",
    )
    gmsh.model.addPhysicalGroup(
        dim=surfaces[3][0],
        tags=[surfaces[3][1], surfaces[4][1], surfaces[5][1]],
        tag=endo_rv_marker,
        name="ENDO_RV",
    )
    gmsh.model.addPhysicalGroup(
        dim=surfaces[0][0],
        tags=[surfaces[0][1], surfaces[2][1]],
        tag=epi_marker,
        name="EPI",
    )

    gmsh.model.add_physical_group(dim=3, tags=[t[1] for t in ov2], tag=1)

    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", char_length)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", char_length)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")

    gmsh.write(path.as_posix())
    gmsh.finalize()

    return path


def biv_ellipsoid_torso(
    mesh_name: str = "",
    heart_as_surface: bool = False,
    torso_length: float = 20.0,
    torso_width: float = 20.0,
    torso_height: float = 20.0,
    rotation_angle: float = math.pi / 6,
    center_lv=(0.0, 0.0, 0.0),
    radius_lv=1.0,
    radius_lv_endo=1.0,
    a_endo_lv: float = 2.5,
    b_endo_lv: float = 1.0,
    c_endo_lv: float = 1.0,
    a_epi_lv: float = 3.0,
    b_epi_lv: float = 1.5,
    c_epi_lv: float = 1.5,
    center_rv=(0.0, 0.5, 0.0),
    radius_rv=1.0,
    radius_rv_endo=1.0,
    a_endo_rv: float = 3.0,
    b_endo_rv: float = 1.5,
    c_endo_rv: float = 1.5,
    a_epi_rv: float = 4.0,
    b_epi_rv: float = 2.5,
    c_epi_rv: float = 2.0,
    angle1=-math.pi / 2,
    angle2=math.pi / 2,
    angle3=2 * math.pi,
    base_x: float = 0.0,
    char_length: float = 0.5,
) -> Path:
    path = utils.handle_mesh_name(mesh_name=mesh_name)
    gmsh.initialize()
    gmsh.model.add("biv")

    # Create a box for cutting the base
    a_epi = max(a_epi_lv, a_epi_rv)
    diam = -5 * a_epi  # Just make it sufficiently big
    box_id = gmsh.model.occ.addBox(base_x, a_epi, a_epi, diam, diam, diam)
    dim_tag_box = [(3, box_id)]

    # LV epicardium
    lv_id = gmsh.model.occ.addSphere(
        *center_lv, radius=radius_lv, angle1=angle1, angle2=angle2, angle3=angle3
    )
    dim_tag_lv = [(3, lv_id)]

    gmsh.model.occ.dilate(dim_tag_lv, *center_lv, a=a_epi_lv, b=b_epi_lv, c=c_epi_lv)

    # LV endocardium
    lv_endo_id = gmsh.model.occ.addSphere(
        *center_lv, radius=radius_lv_endo, angle1=angle1, angle2=angle2, angle3=angle3
    )
    dim_tag_lv_endo = [(3, lv_endo_id)]
    gmsh.model.occ.dilate(
        dim_tag_lv_endo, *center_lv, a=a_endo_lv, b=b_endo_lv, c=c_endo_lv
    )

    # RV epicardium
    rv_id = gmsh.model.occ.addSphere(
        *center_lv, radius=radius_rv, angle1=angle1, angle2=angle2, angle3=angle3
    )
    dim_tag_rv = [(3, rv_id)]

    gmsh.model.occ.dilate(dim_tag_rv, *center_lv, a=a_epi_rv, b=b_epi_rv, c=c_epi_rv)
    gmsh.model.occ.translate(dim_tag_rv, *center_rv)

    # RV endocardium
    rv_endo_id = gmsh.model.occ.addSphere(
        *center_rv, radius=radius_rv_endo, angle1=angle1, angle2=angle2, angle3=angle3
    )

    dim_tag_rv_endo = [(3, rv_endo_id)]
    gmsh.model.occ.dilate(
        dim_tag_rv_endo, *center_lv, a=a_endo_rv, b=b_endo_rv, c=c_endo_rv
    )
    gmsh.model.occ.translate(dim_tag_rv_endo, *center_rv)

    # Subtract LV epicardium from RV endo
    rv_endo, _ = gmsh.model.occ.cut(dim_tag_rv_endo, dim_tag_lv, removeTool=False)

    # Combine the two endocardial volumes
    endo, _ = gmsh.model.occ.fuse(rv_endo, dim_tag_lv_endo)
    # Combine the two epicardial volumes
    epi, _ = gmsh.model.occ.fuse(dim_tag_lv, dim_tag_rv)

    # Subtract the endocardial surfaces from the epicardial volumes
    ov, _ = gmsh.model.occ.cut(epi, endo)
    # And finally cut the base with the box
    ov2, _ = gmsh.model.occ.cut(ov, dim_tag_box)

    torso_tag = gmsh.model.occ.addBox(
        -torso_length / 2.0,
        -torso_width / 2.0,
        -torso_height / 2.0,
        torso_length,
        torso_width,
        torso_height,
        3,
    )

    # Rotate torso around this point to align with realistic heart in torso
    gmsh.model.occ.rotate([(3, torso_tag)], 0, 0, 0, 0, 0, 1, -rotation_angle)
    gmsh.model.occ.rotate([(3, torso_tag)], 0, 0, 0, 0, 1, 0, rotation_angle)

    if heart_as_surface:
        mark_heart_as_surface(ov2, torso_tag)
    else:
        mark_heart_as_volume(ov2, torso_tag)

    # gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", char_length)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", char_length)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")

    gmsh.write(path.as_posix())
    gmsh.finalize()

    return path


def mark_heart_as_surface(ov2, torso_tag):
    surfaces = gmsh.model.occ.getEntities(dim=2)
    volumes = gmsh.model.occ.getEntities(dim=3)

    gmsh.model.occ.cut(
        [(3, torso_tag)],
        ov2,
    )

    gmsh.model.occ.synchronize()

    base_marker = 1
    endo_lv_marker = 2
    endo_rv_marker = 3
    epi_marker = 4

    side1_marker = 5
    side2_marker = 6
    side3_marker = 7
    side4_marker = 8
    side5_marker = 9
    side6_marker = 10

    tissue_marker = 11

    gmsh.model.addPhysicalGroup(
        dim=surfaces[0][0],
        tags=[surfaces[0][1]],
        tag=side1_marker,
        name="TOP",
    )
    gmsh.model.addPhysicalGroup(
        dim=surfaces[1][0],
        tags=[surfaces[1][1]],
        tag=side2_marker,
        name="LEFT",
    )
    gmsh.model.addPhysicalGroup(
        dim=surfaces[2][0],
        tags=[surfaces[2][1]],
        tag=side3_marker,
        name="FRONT",
    )
    gmsh.model.addPhysicalGroup(
        dim=surfaces[3][0],
        tags=[surfaces[3][1]],
        tag=side4_marker,
        name="RIGHT",
    )
    gmsh.model.addPhysicalGroup(
        dim=surfaces[4][0],
        tags=[surfaces[4][1]],
        tag=side5_marker,
        name="BACK",
    )
    gmsh.model.addPhysicalGroup(
        dim=surfaces[5][0],
        tags=[surfaces[5][1]],
        tag=side6_marker,
        name="BOTTOM",
    )
    gmsh.model.addPhysicalGroup(
        dim=surfaces[7][0],
        tags=[surfaces[7][1]],
        tag=base_marker,
        name="BASE",
    )
    gmsh.model.addPhysicalGroup(
        dim=surfaces[8][0],
        tags=[surfaces[6][1], surfaces[8][1]],
        tag=epi_marker,
        name="EPI",
    )
    gmsh.model.addPhysicalGroup(
        dim=surfaces[9][0],
        tags=[surfaces[9][1], surfaces[10][1], surfaces[11][1]],
        tag=endo_rv_marker,
        name="ENDO_RV",
    )

    gmsh.model.addPhysicalGroup(
        dim=surfaces[12][0],
        tags=[surfaces[12][1], surfaces[13][1]],
        tag=endo_lv_marker,
        name="ENDO_LV",
    )

    gmsh.model.add_physical_group(
        dim=3,
        tags=[volumes[0][1]],
        tag=tissue_marker,
        name="TISSUE",
    )


def mark_heart_as_volume(ov2, torso_tag):
    surfaces = gmsh.model.occ.getEntities(dim=2)

    gmsh.model.occ.fuse([(3, torso_tag)], ov2, removeTool=False)

    volumes = gmsh.model.occ.getEntities(dim=3)

    gmsh.model.occ.synchronize()

    base_marker = 1
    endo_lv_marker = 2
    endo_rv_marker = 3
    epi_marker = 4

    side1_marker = 5
    side2_marker = 6
    side3_marker = 7
    side4_marker = 8
    side5_marker = 9
    side6_marker = 10

    tissue_marker = 11
    heart_marker = 12

    gmsh.model.addPhysicalGroup(
        dim=surfaces[8][0],
        tags=[surfaces[8][1]],
        tag=side1_marker,
        name="TOP",
    )
    gmsh.model.addPhysicalGroup(
        dim=surfaces[9][0],
        tags=[surfaces[9][1]],
        tag=side2_marker,
        name="LEFT",
    )
    gmsh.model.addPhysicalGroup(
        dim=surfaces[10][0],
        tags=[surfaces[10][1]],
        tag=side3_marker,
        name="FRONT",
    )
    gmsh.model.addPhysicalGroup(
        dim=surfaces[11][0],
        tags=[surfaces[11][1]],
        tag=side4_marker,
        name="RIGHT",
    )
    gmsh.model.addPhysicalGroup(
        dim=surfaces[12][0],
        tags=[surfaces[12][1]],
        tag=side5_marker,
        name="BACK",
    )
    gmsh.model.addPhysicalGroup(
        dim=surfaces[13][0],
        tags=[surfaces[13][1]],
        tag=side6_marker,
        name="BOTTOM",
    )

    gmsh.model.addPhysicalGroup(
        dim=surfaces[0][0],
        tags=[surfaces[0][1], surfaces[2][1]],
        tag=epi_marker,
        name="EPI",
    )
    gmsh.model.addPhysicalGroup(
        dim=surfaces[1][0],
        tags=[surfaces[1][1]],
        tag=base_marker,
        name="BASE",
    )
    gmsh.model.addPhysicalGroup(
        dim=surfaces[3][0],
        tags=[surfaces[3][1], surfaces[4][1], surfaces[5][1]],
        tag=endo_rv_marker,
        name="ENDO_RV",
    )
    gmsh.model.addPhysicalGroup(
        dim=surfaces[7][0],
        tags=[surfaces[6][1], surfaces[7][1]],
        tag=endo_lv_marker,
        name="ENDO_LV",
    )
    gmsh.model.add_physical_group(
        dim=3,
        tags=[volumes[1][1]],
        tag=tissue_marker,
        name="TISSUE",
    )
    gmsh.model.add_physical_group(
        dim=3,
        tags=[volumes[0][1]],
        tag=heart_marker,
        name="HEART",
    )
