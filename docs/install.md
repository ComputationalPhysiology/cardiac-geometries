# Installation

You can install the software with pip
```
python3 -m pip install cardiac-geometries
```

Alternatively you can install it directly from Github
```
python3 -m pip install git+https://github.com/ComputationalPhysiology/cardiac_geometries
```
However, to actually use the software you need to also install [`gmsh`](http://gmsh.info) with python bindings and to convert the mesh into FEniCS format you also need to install [FEniCS](https://fenicsproject.org/download/archive/).

## Install using Docker
Docker is a simple alternative if you don't want to go through the hazel of installing both FEniCS and Gmsh. We provide a pre-build Docker image with FEniCS and Gmsh installed at https://github.com/orgs/scientificcomputing/packages/container/package/fenics-gmsh, which should be compatible with both AMD64 and ARM64 architectures. You can run the image interactively using the command
```
docker run --rm -v $PWD:/home/shared -w /home/shared -it ghcr.io/scientificcomputing/fenics-gmsh:2023-03-01
```
If you want to get an image with `cardiac-gemometries` installed you can use the Docker image from the GitHub repo
```
docker run --rm -v $PWD:/home/shared -w /home/shared -it ghcr.io/computationalphysiology/cardiac_geometries:latest
```


## Installing python bindings for gmsh
If you are using Mac and have installed `gmsh` with Homebrew, i.e
```
brew install gmsh
```
then it will also install python bindings. In my case there is a file called `gmsh.py` located in `/opt/homebrew/lib`. In order to make python see this file, you can just update the `PYTHONPATH` environment variable, e.g
```
export PYTHONPATH=/opt/homebrew/lib:$PYTHONPATH
```
Similar hacks might be possible for Linux as well
