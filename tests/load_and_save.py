"""Simple script for loading an existing geometry
and saving it to a different file and returning a
non-zero exit code if this fails. This is useful
for making sure that this works in parallel.

To run this, first create a geometry in serial, e.g

.. code::
    cardiac-geometries create-biv-ellipsoid  --create-fibers biv

and then run this script in parallel

.. code::
    mpirun -n 2 python3 tests/load_and_save.py biv/biv_ellipsoid.h5

"""
import logging
import sys
from pathlib import Path

import dolfin
from cardiac_geometries.geometry import Geometry


def unlink_geofile(comm, path: Path) -> None:
    if dolfin.MPI.rank(comm) == 0:
        logging.info("Delete existing file")
        path.unlink(missing_ok=True)
        path.with_suffix(".json").unlink(missing_ok=True)
    else:
        logging.info("Wait for root process to delete file")

    dolfin.MPI.barrier(comm)


def main(path: Path) -> int:
    logging.info(f"Start. Check if {path} exist")
    if not path.is_file():
        logging.info("No!")
        return 1

    logging.info("Loading file")
    geo = Geometry.from_file(path)
    outpath = Path("out.h5")

    comm = geo.comm

    unlink_geofile(comm, path=outpath)
    logging.info("Save geometry to outpath")
    geo.save(outpath)

    if outpath.is_file() and outpath.with_suffix(".json").is_file():
        logging.info("Success! File exist")
        unlink_geofile(comm, path=outpath)
        return 0  # Success
    logging.info("Failure")
    return 1  # Failure


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="[Process %(process)d]: %(message)s")
    raise SystemExit(main(Path.cwd() / Path(sys.argv[1])))
