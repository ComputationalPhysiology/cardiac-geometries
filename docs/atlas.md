---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Cardiac atlases

`cardiac-geometries` also have support for downloading and converting geometries coming from atlases. Currently we only support data coming from the Bai et. al atlas which is found at https://zenodo.org/records/4506463

You can use the command
```{code-cell} shell
!cardiac-geometries atlas-bai --help
```
to work with the atlas data. For example say you want to download data for instance number 4, then you can do

```{code-cell} shell
!cardiac-geometries atlas-bai 4 -o atlas-data
```
which will output the data in the folder `atlas-data`
```{code-cell} shell
!ls -R atlas-data
```
Now the data inside `atlas-data/instance_004/original` can be read by `dolfin`.
