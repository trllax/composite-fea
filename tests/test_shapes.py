from pathlib import Path
from compfea.shapes import circular_centerline, parse_nodes_elements
import math

def test_circular_centerline_90():
    xyz = circular_centerline(length_mm=100.0, theta_deg=90.0, s0=0.0, long_axis="x")
    assert abs(xyz[-1, 0] - 100.0 / (math.pi / 2)) < 1e-6  # r = L/θ, tip span = r
    assert abs(xyz[-1, 2] - 100.0 / (math.pi / 2)) < 1e-6  # z = r at 90

def test_parse_saved_deck():
    p = Path("results/fin_ubend_90_saved/deck.inp")
    if not p.is_file():
        return
    nodes, els = parse_nodes_elements(p)
    assert len(nodes) > 100 and len(els) > 50
