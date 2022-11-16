from typing import Dict
from typing import Tuple


def default_markers() -> Dict[str, Tuple[int, int]]:
    """
    Default markers for the mesh boundaries
    """
    return dict(base=(10, 2), rv=(20, 2), lv=(30, 2), epi=(40, 2))
