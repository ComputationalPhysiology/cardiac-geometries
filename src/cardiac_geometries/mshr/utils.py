from typing import Dict


def default_markers() -> Dict[str, int]:
    """
    Default markers for the mesh boundaries
    """
    return dict(base=10, rv=20, lv=30, epi=40)
