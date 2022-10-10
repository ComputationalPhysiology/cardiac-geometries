FROM finsberg/fenics-gmsh


RUN RUN python3 -m pip install ".[fibers]"
