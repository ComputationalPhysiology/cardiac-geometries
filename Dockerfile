FROM finsberg/fenics-gmsh:latest

COPY . /app
WORKDIR /app

RUN python3 -m pip install --upgrade pip && \
    python3 -m pip install .
