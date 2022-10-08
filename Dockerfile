FROM finsberg/fenics

ARG GMSH_VERSION=4_10_5

# Install gmsh
RUN git clone -b gmsh_${GMSH_VERSION} --single-branch --depth 1 https://gitlab.onelab.info/gmsh/gmsh.git && \
    cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -DENABLE_BUILD_DYNAMIC=1 -B build-dir -S gmsh && \
    cmake --build build-dir && \
    cmake --install build-dir && \
    rm -rf /tmp/*

# GMSH installs python library in /usr/local/lib, see: https://gitlab.onelab.info/gmsh/gmsh/-/issues/1414
ENV PYTHONPATH=/usr/local/lib:$PYTHONPATH

RUN python -m pip install scipy matplotlib meshio h5py --no-binary=h5py

ENV LD_PRELOAD /usr/lib/aarch64-linux-gnu/libgomp.so.1
