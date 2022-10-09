# Cardiac Geometries

This is a library for creating idealized cardiac geometries in FEniCS using `gmsh`. There is also some functionality for doing this with `mshr`, however, `mshr` is not really supported anymore.

Install with `pip`
```
python3 -m pip install cardiac-geometries
```
To run the code you need to have both FEniCS and `gmsh` installed with python bindings. The easiest way to achieve this is by using the docker image `finsberg/fenics-gmsh`. Use the command
```
docker run --rm -v $PWD:/home/shared -w /home/shared -it finsberg/fenics-gmsh
```
to started a Docker container with FEniCS (development version) and gmsh installed. This docker image is also compatible with both AMD64 and ARM64 architectures. The [Dockerfile](Dockerfile) in this repository contains instructions on how to install gmsh with python bindings.

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
