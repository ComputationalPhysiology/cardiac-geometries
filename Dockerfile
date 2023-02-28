FROM ghcr.io/scientificcomputing/fenics-gmsh:2023-02-20

COPY . /app
WORKDIR /app

RUN python3 -m pip install --upgrade pip && \
    python3 -m pip install .
