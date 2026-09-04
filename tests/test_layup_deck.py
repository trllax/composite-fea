"""Unit tests for layup and deck assembly (no ccx)."""

from __future__ import annotations

import math

import pytest

from compfea.deck import StaticStep, assemble, tip_u_clamp_body
from compfea.layup import (
    EngineeringConstants,
    Layup,
    Ply,
    ZoneLayup,
    cross_ply_symmetric,
    orientation_name,
)
from compfea.metrics import f_spring


def test_orientation_names():
    assert orientation_name(0) == "ori_p0"
    assert orientation_name(90) == "ori_p90"
    assert orientation_name(-45) == "ori_m45"


def test_cross_ply_emits_orientation_names_not_angles():
    layup = cross_ply_symmetric(0.1)
    text = layup.to_inp()
    assert "*SHELL SECTION, COMPOSITE, ELSET=blade" in text
    assert "ori_p0" in text and "ori_p90" in text
    section = text.split("COMPOSITE")[1]
    assert ", 0," not in section
    assert text.count("0.1, , cfrp,") == 4
    assert "ENGINEERING CONSTANTS" in text


def test_complex_zone_layup():
    mat = EngineeringConstants(
        135000, 9000, 9000, 0.3, 0.3, 0.45, 4500, 4500, 3000, name="t700"
    )
    layup = Layup(
        materials=(mat,),
        zones=(
            ZoneLayup(
                "zone_root",
                (
                    Ply(0.1, 0, "t700"),
                    Ply(0.1, 45, "t700"),
                    Ply(0.1, -45, "t700"),
                    Ply(0.1, 90, "t700"),
                    Ply(0.1, 90, "t700"),
                    Ply(0.1, -45, "t700"),
                    Ply(0.1, 45, "t700"),
                    Ply(0.1, 0, "t700"),
                ),
            ),
            ZoneLayup(
                "zone_tip",
                (Ply(0.1, 0, "t700"), Ply(0.1, 90, "t700"), Ply(0.1, 0, "t700")),
            ),
        ),
        zero_along="y",
    )
    text = layup.to_inp()
    assert "zone_root" in text and "zone_tip" in text
    assert "ori_p45" in text and "ori_m45" in text
    assert text.count("*SHELL SECTION, COMPOSITE") == 2


def test_deck_assemble_orders_blocks():
    layup = cross_ply_symmetric(0.1)
    mesh = (
        "*NODE\n1, 0, 0, 0\n"
        "*ELEMENT, TYPE=S8R, ELSET=blade\n"
        "1, 1, 1, 1, 1, 1, 1, 1, 1"
    )
    step = StaticStep(tip_u_clamp_body({10: (0.0, -1.0, 2.0)}), inc=100)
    deck = assemble(
        mesh_inp=mesh,
        layup=layup,
        initial_bc="*BOUNDARY\nfixed_end, 1, 6",
        steps=[step],
        heading="unit test deck",
    )
    assert deck.index("*NODE") < deck.index("*MATERIAL")
    assert deck.index("*MATERIAL") < deck.index("*SHELL SECTION")
    assert deck.index("*SHELL SECTION") < deck.index("*STEP")
    assert "*EL PRINT, ELSET=blade, TOTALS=ONLY" in deck
    assert "ELSE" in deck
    assert "10, 3, 3, 2.0000000000" in deck


def test_f_spring_doubles_with_theta():
    c = 100.0
    arm = 100.0
    th = math.pi / 2
    u90 = 0.5 * c * th**2
    u180 = 0.5 * c * (2 * th) ** 2
    _m90, f90 = f_spring(u90, th, arm)
    _m180, f180 = f_spring(u180, 2 * th, arm)
    assert f180 / f90 == pytest.approx(2.0)
