import contextlib
import warnings
from pathlib import Path
from textwrap import dedent
from typing import List
from typing import Tuple
from typing import Union

import dolfin
import numpy as np
import ufl


def value_size(obj: ufl.Coefficient) -> Union[List[int], int]:
    value_shape = obj.value_shape()
    if len(value_shape) == 0:
        return 1
    return [0]


body = dedent(
    """<?xml version="1.0"?>
    <Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">
        <Domain>
            {body}
        </Domain>
    </Xdmf>""",
)

series = dedent(
    """
    <Grid Name="{name}" GridType="Collection" CollectionType="Temporal">
        <Time TimeType="List">
            <DataItem Format="XML" Dimensions="{N}"> {lst}</DataItem>
        </Time>
    {entry}
    </Grid>
    """,
)


entry_single = dedent(
    """
    <Grid Name="time_{iter}" GridType="Uniform">
        {frame}
    </Grid>
    """,
)

entry = dedent(
    """
    <Grid Name="time_{iter}" GridType="Uniform">
        {frame}
    </Grid>
    """,
)


topology = dedent(
    """
    <Topology NumberOfElements="{ncells}" TopologyType="{cell}">
        <DataItem Dimensions="{ncells} {dim}" Format="HDF">{h5name}:/{h5group}</DataItem>
    </Topology>
    """,
)


topology_polyvert = dedent(
    """
    <Topology TopologyType="Polyvertex" NodesPerElement="{nverts}">
    </Topology>
    """,
)

geometry = dedent(
    """
    <Geometry GeometryType="{coords}">
        <DataItem Dimensions="{nverts} {dim}" Format="HDF">{h5name}:/{h5group}</DataItem>
    </Geometry>
    """,
)

vector_attribute = dedent(
    """
    <Attribute Name="{name}" AttributeType="Vector" Center="{center}">
        <DataItem Format="HDF" Dimensions="{nverts} {dim}">{h5name}:/{h5group}</DataItem>
    </Attribute>
    """,
)

scalar_attribute = dedent(
    """
    <Attribute Name="{name}" AttributeType="Scalar" Center="{center}">
        <DataItem Format="HDF" Dimensions="{nverts} {dim}">{h5name}:/{h5group}</DataItem>
    </Attribute>
    """,
)


@contextlib.contextmanager
def h5pyfile(h5name, filemode="r", comm=None):
    try:
        import h5py
    except ImportError:
        print("Please install h5py")
        print("python3 -m pip install h5py --no-binary=h5py")
        raise

    if comm is None:
        comm = dolfin.MPI.comm_world

    if h5py.h5.get_config().mpi and dolfin.MPI.size(comm) > 1:
        h5file = h5py.File(h5name, filemode, driver="mpio", comm=comm)
    else:
        if dolfin.MPI.size(comm) > 1:
            warnings.warn("h5py is not installed with MPI support")
        h5file = h5py.File(h5name, filemode)
    yield h5file

    h5file.close()


def dict_to_h5(data, h5name, h5group, comm=None):
    if comm is None:
        comm = dolfin.MPI.comm_world
    if comm.rank == 0:
        import h5py

        with h5py.File(h5name, "a") as h5file:
            if h5group == "":
                group = h5file
            else:
                group = h5file.create_group(h5group)
            for k, v in data.items():
                group.create_dataset(k, data=v)

    # Wait for root process to finish
    dolfin.MPI.barrier(comm)


def decode(x):
    if isinstance(x, bytes):
        return x.decode()
    elif isinstance(x, list):
        return [xi for xi in map(decode, x)]
    return x


def h5_to_dict(h5group):
    import h5py

    group = {}
    for key, value in h5group.items():
        if isinstance(value, h5py.Dataset):
            v = decode(value[...].tolist())
            group[key] = v

        # elif isinstance(value, h5py.Group):
        #     group[key] = h5_to_dict(h5group[key])

        # else:
        #     raise ValueError(f"Unknown value type {type(value)} for key {key}.")

    return group


def dolfin_to_hd5(obj: dolfin.Function, h5name, time="", name=None):
    """
    Save object to and HDF file.

    Parameters
    ----------
    obj : dolfin.Function
        The object you want to save
    name : str
        Name of the object
    h5group : str
        The folder you want to save the object
        withing the HDF file. Default: ''

    """
    assert isinstance(obj, dolfin.Function)
    name = obj.name() if name is None else name
    dolfin.info("Save {0} to {1}:{0}/{2}".format(name, h5name, time))

    group = name if time == "" else "/".join([name, str(time)])
    file_mode = "a" if Path(h5name).is_file() else "w"

    if value_size(obj) == 1:
        return save_scalar_function(obj, h5name, group, file_mode)
    else:
        return save_vector_function(obj, h5name, group, file_mode)


def save_scalar_function(obj, h5name, h5group="", file_mode="w"):
    V = obj.function_space()

    dim = V.mesh().geometry().dim()
    # TODO: gather
    coords_tmp = V.tabulate_dof_coordinates()
    coords = coords_tmp.reshape((-1, dim))
    # TODO: gather
    obj_arr = obj.vector().get_local()
    vecs = np.array(obj_arr).T

    coord_group = "/".join([h5group, "coordinates"])
    vector_group = "/".join([h5group, "vector"])
    with h5pyfile(h5name, file_mode) as h5file:
        if h5group in h5file:
            del h5file[h5group]

        h5file.create_dataset(coord_group, data=coords)
        h5file.create_dataset(vector_group, data=vecs)

    element = obj.ufl_element()

    if dim == 3:
        coords = "XYZ"
    elif dim == 2:
        coords = "XY"
    else:
        coords = "X"

    return {
        "h5group": h5group,
        "coordinates": coord_group,
        "vector": vector_group,
        "nverts": obj.vector().size(),
        "dim": 1,
        "family": element.family(),
        "geo_dim": dim,
        "coords": coords,
        "degree": element.degree(),
        "type": "scalar",
    }


def save_vector_function(obj, h5name, h5group="", file_mode="w"):
    V = obj.function_space()
    gs = obj.split(deepcopy=True)

    W = V.sub(0).collapse()
    dim = V.mesh().geometry().dim()
    # TODO: gather
    coords_tmp = W.tabulate_dof_coordinates()
    coords = coords_tmp.reshape((-1, dim))
    # TODO: gather
    us = [g.vector().get_local() for g in gs]
    vecs = np.array(us).T

    coord_group = "/".join([h5group, "coordinates"])
    vector_group = "/".join([h5group, "vector"])

    with h5pyfile(h5name, file_mode) as h5file:
        if h5group in h5file:
            del h5file[h5group]
        h5file.create_dataset(coord_group, data=coords)
        h5file.create_dataset(vector_group, data=vecs)

    element = obj.ufl_element()

    if dim == 3:
        coords = "XYZ"
    elif dim == 2:
        coords = "XY"
    else:
        coords = "X"

    return {
        "h5group": h5group,
        "coordinates": coord_group,
        "vector": vector_group,
        "nverts": obj.vector().size() / dim,
        "dim": dim,
        "family": element.family(),
        "geo_dim": dim,
        "coords": coords,
        "degree": element.degree(),
        "type": "vector",
    }


def fun_to_xdmf(fun, fname, name="function"):
    h5name = Path(fname).with_suffix(".h5")
    dolfin_to_hd5(fun, h5name, name=name)

    dim = fun.function_space().mesh().geometry().dim()

    if value_size(fun) == 1:
        nverts = len(fun.vector())
        fun_str = scalar_attribute.format(
            name=name,
            nverts=nverts,
            center="Node",
            h5group="/".join([name, "vector"]),
            dim=1,
            h5name=h5name.name,
        )
    else:
        nverts = int(len(fun.vector()) / dim)
        fun_str = vector_attribute.format(
            name=name,
            nverts=nverts,
            dim=dim,
            h5group="/".join([name, "vector"]),
            center="Node",
            h5name=h5name.name,
        )

    fun_top = topology_polyvert.format(nverts=nverts)
    fun_geo = geometry.format(
        nverts=nverts,
        dim=dim,
        coords="XYZ",
        h5group="/".join([name, "coordinates"]),
        h5name=h5name.name,
    )

    fun_entry = entry.format(frame=fun_geo + fun_top + fun_str, iter=0)
    T = body.format(body=fun_entry, name="Visualzation of {}".format(name))

    Path(fname).with_suffix(".xdmf").write_text(T)


def get_sign_direction(base_direction: str) -> Tuple[int, str]:
    if len(base_direction) == 1:
        return 1, base_direction
    if len(base_direction) == 2:
        if base_direction[0] == "+":
            return 1, base_direction[1]
        elif base_direction[0] == "-":
            return -1, base_direction[1]
    raise ValueError(f"Invalid base direction {base_direction!r}")


def fiber_to_xdmf(fun: dolfin.Function, fname, base_direction="z"):
    h5name = Path(fname).with_suffix(".h5")
    dolfin_to_hd5(fun, h5name, name="fiber")

    sign, direction = get_sign_direction(base_direction=base_direction)
    index = {"x": 0, "y": 1, "z": 2}[direction]

    f = fun.split(deepcopy=True)[index]
    f_arr = f.vector().get_local()
    scalar = np.arcsin(sign * f_arr) * 180 / np.pi

    with h5pyfile(h5name, "a") as h5file:
        h5file.create_dataset("fiber/scalar", data=scalar)

    dim = fun.function_space().mesh().geometry().dim()
    nverts = int(fun.vector().size() / dim)
    name = "fiber"

    fun_scal = scalar_attribute.format(
        name="angle",
        nverts=nverts,
        center="Node",
        h5group="/".join([name, "scalar"]),
        dim=1,
        h5name=h5name.name,
    )

    fun_vec = vector_attribute.format(
        name=name,
        nverts=nverts,
        dim=dim,
        center="Node",
        h5group="/".join([name, "vector"]),
        h5name=h5name.name,
    )

    fun_top = topology_polyvert.format(nverts=nverts)
    fun_geo = geometry.format(
        nverts=nverts,
        dim=dim,
        coords="XYZ",
        h5group="/".join([name, "coordinates"]),
        h5name=h5name.name,
    )

    fun_entry = entry_single.format(
        frame=fun_geo + fun_top + fun_scal + fun_vec,
        iter=0,
    )
    T = body.format(body=fun_entry, name="Visualzation of {}".format(name))

    with open("{}.xdmf".format(fname), "w") as f:
        f.write(T)
