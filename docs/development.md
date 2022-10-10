# Development installation

Developers should install some extra dependencies
```
python3 -m pip install install -e ".[dev,docs,test]"
```
as well as the pre-commit hook
```
pre-commit install
```
Alternatively, you can use the `Makefile` and hit `make dev`.
