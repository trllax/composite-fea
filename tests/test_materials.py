"""Materials library: units, admissibility, and the round trip to the constants."""

from __future__ import annotations

from pathlib import Path

import pytest

from compfea.layup import (
    ANSYS_EPOXY_CARBON_WOVEN_230_WET,
    UD_CFRP_GENERIC,
    EngineeringConstants,
    woven_from_ud,
)
from compfea.materials import (
    DENSITY_TO_TONNE_PER_MM3,
    STRESS_TO_MPA,
    MaterialError,
    check_admissible,
    constants_for,
    load_materials,
    parse_units_header,
)

ROOT = Path(__file__).resolve().parents[1]
GENERIC = ROOT / "materials" / "generic.csv"

HEADER = "name,kind,e1,e2,e3,nu12,nu13,nu23,g12,g13,g23,density,xt,xc,yt,yc,s12,source"
# One physically sane UD row, in repo units, as a formatting template.
MPA_ROW = "cf,ud,135000,9000,9000,0.3,0.3,0.45,4500,4500,3000,1.6e-9,,,,,,test"
# The same lamina as ANSYS would export it: Pa and kg/m3.
PA_ROW = "cf,ud,1.35e11,9e9,9e9,0.3,0.3,0.45,4.5e9,4.5e9,3e9,1600,,,,,,test"


def write(tmp_path: Path, *lines: str) -> Path:
    path = tmp_path / "mat.csv"
    path.write_text("\n".join(lines) + "\n")
    return path


# --------------------------------------------------------------------------
# units


def test_a_file_with_no_units_line_is_refused(tmp_path):
    """The single most expensive thing to get wrong, so it has no default."""
    path = write(tmp_path, HEADER, MPA_ROW)
    with pytest.raises(MaterialError, match="no '# units:' line"):
        load_materials(path)


def test_a_partial_units_line_is_refused(tmp_path):
    path = write(tmp_path, "# units: stress=MPa", HEADER, MPA_ROW)
    with pytest.raises(MaterialError, match="missing \\['density'\\]"):
        load_materials(path)


def test_an_unknown_unit_is_refused(tmp_path):
    path = write(tmp_path, "# units: stress=bar density=kg/m3", HEADER, MPA_ROW)
    with pytest.raises(MaterialError, match="unknown stress unit"):
        load_materials(path)


def test_pa_and_mpa_files_describe_the_same_lamina(tmp_path):
    """The conversion is the feature: ANSYS data must paste in unconverted."""
    a = load_materials(
        write(tmp_path, "# units: stress=Pa density=kg/m3", HEADER, PA_ROW)
    )
    tmp2 = tmp_path / "b"
    tmp2.mkdir()
    b = load_materials(
        write(tmp2, "# units: stress=MPa density=tonne/mm3", HEADER, MPA_ROW)
    )
    for field in ("e1", "e2", "e3", "g12", "g13", "g23", "density"):
        assert getattr(a["cf"].constants, field) == pytest.approx(
            getattr(b["cf"].constants, field), rel=1e-12
        )
    assert a["cf"].constants.e1 == pytest.approx(135000.0)
    assert a["cf"].constants.density == pytest.approx(1.6e-9)


def test_unit_scales_are_self_consistent():
    assert STRESS_TO_MPA["mpa"] == 1.0
    assert STRESS_TO_MPA["pa"] * 1e11 == pytest.approx(1e5)
    assert STRESS_TO_MPA["gpa"] * 121.0 == pytest.approx(121000.0)
    assert DENSITY_TO_TONNE_PER_MM3["kg/m3"] * 1600.0 == pytest.approx(1.6e-9)
    assert DENSITY_TO_TONNE_PER_MM3["g/cm3"] * 1.6 == pytest.approx(1.6e-9)


def test_parse_units_header_stops_at_the_first_data_line():
    scales = parse_units_header(["# units: stress=GPa density=g/cm3", HEADER, MPA_ROW])
    assert scales == (1e3, 1e-9)


def test_a_file_left_in_pa_but_declared_mpa_is_caught(tmp_path):
    """The mistake the range guard exists for: right numbers, wrong declaration."""
    path = write(tmp_path, "# units: stress=MPa density=tonne/mm3", HEADER, PA_ROW)
    with pytest.raises(MaterialError, match="outside the plausible band"):
        load_materials(path)


def test_a_file_left_in_mpa_but_declared_pa_is_caught(tmp_path):
    path = write(tmp_path, "# units: stress=Pa density=kg/m3", HEADER, MPA_ROW)
    with pytest.raises(MaterialError, match="outside the plausible band"):
        load_materials(path)


# --------------------------------------------------------------------------
# admissibility


def ec(**over) -> EngineeringConstants:
    base = dict(
        e1=135000.0, e2=9000.0, e3=9000.0, nu12=0.30, nu13=0.30, nu23=0.45,
        g12=4500.0, g13=4500.0, g23=3000.0, density=1.6e-9, name="probe",
    )
    base.update(over)
    return EngineeringConstants(**base)


def test_the_shipped_card_is_admissible():
    assert check_admissible(UD_CFRP_GENERIC) == []
    assert check_admissible(woven_from_ud(UD_CFRP_GENERIC)) == []


def test_a_poisson_ratio_past_the_stability_limit_is_refused():
    """|nu12| < sqrt(E1/E2); ccx solves an inadmissible card without complaint.

    nu13 = nu23 = 0 so this pairwise bound and the determinant condition
    coincide -- otherwise the determinant is strictly the tighter of the two and
    would fire first, leaving the pairwise check untested.
    """
    limit = (135000.0 / 9000.0) ** 0.5   # 3.873
    assert check_admissible(ec(nu12=limit * 0.99, nu13=0.0, nu23=0.0)) == []
    with pytest.raises(MaterialError, match="must be <"):
        check_admissible(ec(nu12=limit * 1.01, nu13=0.0, nu23=0.0))


def test_a_negative_or_zero_modulus_is_refused():
    with pytest.raises(MaterialError, match="must be finite and > 0"):
        check_admissible(ec(e2=0.0))
    with pytest.raises(MaterialError, match="must be finite and > 0"):
        check_admissible(ec(g12=-4500.0))


def test_the_determinant_condition_is_checked_independently():
    """A card can pass every pairwise |nu_ij| test and still be inadmissible."""
    bad = ec(e1=10000.0, e2=10000.0, e3=10000.0, nu12=0.7, nu13=0.7, nu23=0.7)
    for pair in ("nu12", "nu13", "nu23"):
        assert abs(getattr(bad, pair)) < 1.0    # every pairwise test passes
    with pytest.raises(MaterialError, match="not positive definite"):
        check_admissible(bad)


def test_a_minor_poisson_ratio_in_the_major_slot_warns_but_loads():
    """Legal, plausible, and almost certainly a transcription error."""
    warns = check_admissible(ec(nu12=0.02))
    assert len(warns) == 1
    assert "minor ratio" in warns[0]


def test_a_ud_card_labelled_woven_warns(tmp_path):
    path = write(
        tmp_path, "# units: stress=MPa density=tonne/mm3", HEADER,
        MPA_ROW.replace("cf,ud,", "cf,woven,"),
    )
    load_materials(path)
    assert any("looks like" in w for w in load_materials.warnings)
    with pytest.raises(MaterialError, match="kind=woven"):
        load_materials(path, on_warning="raise")


def test_density_outside_the_band_is_refused():
    with pytest.raises(MaterialError, match="density"):
        check_admissible(ec(density=1600.0))       # kg/m3 left unconverted


# --------------------------------------------------------------------------
# schema


def test_a_duplicate_name_is_refused(tmp_path):
    path = write(
        tmp_path, "# units: stress=MPa density=tonne/mm3", HEADER, MPA_ROW, MPA_ROW
    )
    with pytest.raises(MaterialError, match="duplicate material name"):
        load_materials(path)


def test_a_missing_source_is_refused(tmp_path):
    path = write(
        tmp_path, "# units: stress=MPa density=tonne/mm3", HEADER,
        MPA_ROW[: MPA_ROW.rindex(",") + 1],
    )
    with pytest.raises(MaterialError, match="source is required"):
        load_materials(path)


def test_a_missing_column_is_refused(tmp_path):
    path = write(
        tmp_path, "# units: stress=MPa density=tonne/mm3",
        HEADER.replace(",g23", ""), MPA_ROW,
    )
    with pytest.raises(MaterialError, match="missing required columns"):
        load_materials(path)


def test_an_unknown_kind_is_refused(tmp_path):
    path = write(
        tmp_path, "# units: stress=MPa density=tonne/mm3", HEADER,
        MPA_ROW.replace("cf,ud,", "cf,prepreg,"),
    )
    with pytest.raises(MaterialError, match="kind must be one of"):
        load_materials(path)


def test_allowables_are_optional_and_scaled(tmp_path):
    row = MPA_ROW.replace(",,,,,,test", ",2231,1082,29,100,60,test")
    path = write(tmp_path, "# units: stress=MPa density=tonne/mm3", HEADER, row)
    rec = load_materials(path)["cf"]
    assert rec.allowables() == {
        "xt": 2231.0, "xc": 1082.0, "yt": 29.0, "yc": 100.0, "s12": 60.0
    }
    pa = MPA_ROW.replace(",,,,,,test", ",2.231e9,1.082e9,2.9e7,1e8,6e7,test")
    tmp2 = tmp_path / "b"
    tmp2.mkdir()
    scaled = load_materials(
        write(tmp2, "# units: stress=Pa density=kg/m3",
              HEADER, pa.replace("135000,9000,9000", "1.35e11,9e9,9e9")
                        .replace("4500,4500,3000", "4.5e9,4.5e9,3e9")
                        .replace(",1.6e-9,", ",1600,"))
    )["cf"]
    assert scaled.xt == pytest.approx(2231.0)


def test_a_blank_allowable_is_none_not_zero(tmp_path):
    path = write(tmp_path, "# units: stress=MPa density=tonne/mm3", HEADER, MPA_ROW)
    assert load_materials(path)["cf"].allowables() == {
        "xt": None, "xc": None, "yt": None, "yc": None, "s12": None
    }


def test_a_negative_allowable_is_refused(tmp_path):
    """Compressive strengths are magnitudes here, not signed stresses."""
    row = MPA_ROW.replace(",,,,,,test", ",2231,-1082,29,100,60,test")
    path = write(tmp_path, "# units: stress=MPa density=tonne/mm3", HEADER, row)
    with pytest.raises(MaterialError, match="positive magnitude"):
        load_materials(path)


# --------------------------------------------------------------------------
# the shipped library


def test_generic_csv_round_trips_the_module_constants():
    """The file replaces two module constants; it must not drift from them."""
    lib = load_materials(GENERIC, on_warning="raise")
    assert lib["cfrp"].constants == UD_CFRP_GENERIC
    assert lib["cfrp_woven"].constants == woven_from_ud(UD_CFRP_GENERIC)
    assert lib["cfrp"].kind == "ud" and lib["cfrp_woven"].kind == "woven"
    for record in lib.values():
        assert record.source, "every row needs provenance"
        assert record.allowables() == dict.fromkeys(
            ("xt", "xc", "yt", "yc", "s12")
        ), "textbook stiffnesses ship with no invented strengths"


def test_the_ansys_library_loads_and_converts():
    """ansys_composites.csv carries the ANSYS Epoxy Carbon library, in Pa.

    The rows are pasted verbatim from ematerials_all.xml, so this exercises the
    Pa -> MPa / kg-m^3 -> tonne-mm^3 conversion and check_admissible on real
    ANSYS data. Any admissibility warning is an error here (on_warning="raise").
    """
    path = ROOT / "materials" / "ansys_composites.csv"
    lib = load_materials(path, on_warning="raise")
    assert set(lib) == {
        "ansys_epoxy_carbon_ud_230_wet",
        "ansys_epoxy_carbon_ud_395_prepreg",
        "ansys_epoxy_carbon_woven_230_wet",
        "ansys_epoxy_carbon_woven_395_prepreg",
        "im7_ud_276",
        "im7_woven_276",
    }
    for name, record in lib.items():
        assert record.source, "every row needs provenance"
        # allowables came across from the ANSYS Stress Limits, as positive MPa
        assert all(v is not None and v > 0 for v in record.allowables().values())
        is_ansys_stock = name.startswith("ansys_")
        assert ("NOT an ANSYS stock card" in record.source) != is_ansys_stock

    # Round trip against the hand-entered module constant. This is blind to a
    # 13<->23 swap (that card has nu13 == nu23 and g13 == g23), so it does not
    # pin the axis mapping on its own.
    assert lib["ansys_epoxy_carbon_woven_230_wet"].constants == (
        ANSYS_EPOXY_CARBON_WOVEN_230_WET
    )

    # The UD 230 Wet row is where the mapping is actually visible: nu_XZ != nu_YZ
    # and G_XZ != G_YZ in the XML, so getting nu13/nu23 and g13/g23 right is a
    # real constraint. XML: nu_XZ=0.27, nu_YZ=0.42; G_XZ=5.0e9, G_YZ=3.08e9.
    ud = lib["ansys_epoxy_carbon_ud_230_wet"].constants
    assert (ud.nu13, ud.nu23) == (0.27, 0.42)
    assert (ud.g13, ud.g23) == (5000.0, 3080.0)


def test_constants_for_is_deterministic_and_only_what_is_used():
    lib = load_materials(GENERIC)
    got = constants_for(lib, ["cfrp_woven", "cfrp", "cfrp"])
    assert [c.name for c in got] == ["cfrp", "cfrp_woven"]
    assert constants_for(lib, ["cfrp"]) == (UD_CFRP_GENERIC,)
    with pytest.raises(MaterialError, match="not in the library"):
        constants_for(lib, ["basalt_ud"])
