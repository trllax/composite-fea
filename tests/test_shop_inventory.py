"""Shop inventory: default lamina fill; fibre moduli stay out of e1."""

from __future__ import annotations

from pathlib import Path

import pytest

from compfea.layup import ANSYS_EPOXY_CARBON_UD_230, ANSYS_EPOXY_CARBON_WOVEN_230_WET
from compfea.shop_inventory import (
    ShopInventoryError,
    default_lamina,
    load_shop_inventory,
    materials_library,
    write_materials_csv,
)

ROOT = Path(__file__).resolve().parents[1]
SHOP = ROOT / "materials" / "shop_inventory.csv"


def test_shipped_shop_matches_im_cards_and_keeps_fibre_out_of_e1():
    skus = load_shop_inventory(SHOP)
    required = {
        "hexcel_uni_379",
        "hexcel_uni_231",
        "twill_3k_198",
        "hexcel_im2_uni_193",
        "hexcel_himax_biax_100",
        "hexcel_im7_twill_205",
    }
    assert required <= set(skus)
    # IM / T700 SKUs carry explicit lamina cards (not blank-E defaults).
    assert not skus["hexcel_im2_uni_193"].used_default_lamina
    assert not skus["hexcel_im7_twill_205"].used_default_lamina
    assert not skus["hexcel_himax_biax_100"].used_default_lamina
    # thickness_mm is the as-laid cured value = gsm/1000 * CURED_PLY_FACTOR (1.2).
    im2 = skus["hexcel_im2_uni_193"]
    assert im2.architecture == "ud" and im2.thickness_mm == pytest.approx(0.2316)
    assert im2.record.constants.e1 == pytest.approx(153000.0)
    biax = skus["hexcel_himax_biax_100"]
    assert biax.architecture == "biax" and biax.record.kind == "woven"
    assert biax.thickness_mm == pytest.approx(0.12)
    assert biax.record.constants.e1 == pytest.approx(59160.0)
    im7 = skus["hexcel_im7_twill_205"]
    assert im7.architecture == "0_90" and im7.thickness_mm == pytest.approx(0.246)
    assert im7.record.constants.e1 == pytest.approx(86000.0)
    # hexcel_uni_379 / _231 are an IM fibre roll (~43 Msi) -> the vetted IM
    # lamina constants (reused from hexcel_im2_uni_193), not the SM default card.
    ud = skus["hexcel_uni_379"]
    assert ud.record.kind == "ud"
    assert ud.thickness_mm == pytest.approx(0.4548)
    assert not ud.used_default_lamina
    assert ud.record.constants.e1 == pytest.approx(153000.0)
    assert ud.record.constants.e1 != pytest.approx(ANSYS_EPOXY_CARBON_UD_230.e1)
    twill = skus["twill_3k_198"]
    assert twill.record.kind == "woven"
    assert twill.fibre_e_msi == pytest.approx(33.5)
    assert twill.fibre_xt_ksi == pytest.approx(653.0)
    assert twill.record.constants.e1 == pytest.approx(
        ANSYS_EPOXY_CARBON_WOVEN_230_WET.e1
    )
    assert twill.record.constants.e1 != pytest.approx(33.5)
    assert twill.record.constants.e1 != pytest.approx(33.5e3)


def test_oz_only_converts_to_gsm(tmp_path):
    path = tmp_path / "shop.csv"
    path.write_text(
        "\n".join(
            [
                "# units: stress=GPa density=g/cm3",
                "name,product,architecture,tow,areal_weight_gsm,areal_weight_oz_yd2,"
                "thickness_mm,width_in,fibre_e_msi,fibre_xt_ksi,e1,e2,e3,nu12,nu13,"
                "nu23,g12,g13,g23,density,qty_m2,source",
                "u,Hex,ud,,,11,,,,,,,,,,,,,,",
            ]
        )
        + "\n"
    )
    sku = load_shop_inventory(path)["u"]
    assert sku.areal_weight_gsm == pytest.approx(11 * 33.9057)
    # no explicit thickness -> derived gsm/1000 * CURED_PLY_FACTOR
    assert sku.thickness_mm == pytest.approx(sku.areal_weight_gsm / 1000 * 1.2)


def test_partial_elastic_row_is_refused(tmp_path):
    path = tmp_path / "shop.csv"
    path.write_text(
        "\n".join(
            [
                "# units: stress=GPa density=g/cm3",
                "name,architecture,areal_weight_gsm,e1,e2,e3,nu12,nu13,nu23,"
                "g12,g13,g23,density,source",
                "bad,ud,200,135,,,,,,,,,,",
            ]
        )
        + "\n"
    )
    with pytest.raises(ShopInventoryError, match="blank"):
        load_shop_inventory(path)


def test_write_materials_csv_round_trip(tmp_path):
    skus = load_shop_inventory(SHOP)
    out = write_materials_csv(skus, tmp_path / "mat.csv")
    from compfea.materials import load_materials

    lib = load_materials(out)
    assert set(lib) == set(materials_library(skus))
    assert lib["hexcel_uni_379"].constants.e1 == pytest.approx(153000.0)


def test_default_lamina_names():
    ud = default_lamina("ud", "x")
    assert ud.name == "x"
    assert ud.e1 == pytest.approx(ANSYS_EPOXY_CARBON_UD_230.e1)
    w = default_lamina("0_90", "y")
    assert w.e1 == pytest.approx(w.e2)
    assert w.e1 == pytest.approx(ANSYS_EPOXY_CARBON_WOVEN_230_WET.e1)
