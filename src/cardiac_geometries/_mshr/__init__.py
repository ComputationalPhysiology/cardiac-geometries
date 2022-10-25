import warnings

from ._biv import create_biv_mesh
from ._lv import create_lv_mesh

warnings.warn(
    "mshr is not supported anymore and will be removed in a future release",
    category=DeprecationWarning,
)

__all__ = ["create_biv_mesh", "create_lv_mesh"]
