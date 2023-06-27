import pytest
from cardiac_geometries import has_mshr
from cardiac_geometries import mshr


require_mshr = pytest.mark.skipif(
    not has_mshr(),
    reason="mshr is required to run the test",
)


@require_mshr
@pytest.mark.mshr
def test_create_lv_mesh():
    geo = mshr.create_lv_mesh()
    assert geo is not None


@require_mshr
@pytest.mark.mshr
def test_create_biv_mesh():
    geo = mshr.create_biv_mesh()
    assert geo is not None
