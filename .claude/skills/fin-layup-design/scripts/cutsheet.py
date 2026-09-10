"""Cut sheet from a ply-book CSV.

  python cutsheet.py cases/fin_test_3/ax/ax2r.csv [more.csv ...]

Per-ply cut list + per-fabric fabric area, nesting allowance, roll yield and dry
mass. Zone planform areas are for FIN_TEST_3.step at the 16 mm mesh (from
geometry.zone_report); regenerate ZONE_AREA_MM2 for a different STEP with

  compfea-build --step X.step --plybook any.csv --materials m.csv \
      --long-axis y --size-mm 16 --out /tmp/z    # prints the zone_report table

Importable: cut_rows(path) -> (ply_rows, fabric_summary, dry_mass_g).
"""
from __future__ import annotations
import csv, sys, pathlib

ZONE_AREA_MM2 = {
    "FULL": 137512.8, "z_3_4ths": 96498.7, "MID": 60705.4,
    "QUARTER": 29984.4, "TIP": 19268.7, "HEAL": 13885.3,
}
# gsm and roll width per fabric in materials/shop_inventory.csv (incl. the older
# IM2/IM7 rows the baseline ply book uses). Extend as stock changes; an unknown
# fabric raises rather than guessing. Roll width None -> yield not estimated.
GSM = {"twill_3k_198": 198.0, "hexcel_uni_231": 231.0, "hexcel_uni_379": 379.0,
       "hexcel_himax_biax_100": 100.0, "hexcel_im7_twill_205": 205.0,
       "hexcel_im2_uni_193": 193.0}
ROLL_W_IN = {"twill_3k_198": 33.5, "hexcel_uni_231": 24.0, "hexcel_uni_379": 24.0,
             "hexcel_himax_biax_100": None, "hexcel_im7_twill_205": None,
             "hexcel_im2_uni_193": None}
IN2MM = 25.4
NEST = 1.35   # offcut / nesting factor for a ply cut to a tapered planform
# wet hand layup runs ~40-50% resin by weight, so cured part mass is roughly
# (dry fabric mass) / 0.55. This is a range, not the x1.2 thickness ratio.
CURED_MASS_OVER_DRY = 1.8


def _rows(path):
    txt = pathlib.Path(path).read_text().splitlines()
    return list(csv.DictReader(l for l in txt
                               if l.strip() and not l.lstrip().startswith("#")))


def cut_rows(path):
    plies, per_mat, mass = [], {}, 0.0
    for r in _rows(path):
        z, m = r["zone"], r["material"]
        if m not in GSM:
            raise SystemExit(f"{path}: unknown fabric {m!r} - add it to cutsheet.GSM")
        area = ZONE_AREA_MM2[z]
        ang = int(float(r["angle_deg"]))
        plies.append(dict(ply=r["ply"], zone=z, material=m, angle=ang,
                          kind=r["kind"], area_mm2=area))
        d = per_mat.setdefault(m, dict(plies=0, area_mm2=0.0, angles=set()))
        d["plies"] += 1
        d["area_mm2"] += area
        d["angles"].add(ang)
        mass += area * GSM[m] * 1e-6
    for m, d in per_mat.items():
        d["fabric_m2"] = d["area_mm2"] / 1e6
        d["nest_m2"] = d["fabric_m2"] * NEST
        w = ROLL_W_IN[m]
        d["roll_m"] = d["nest_m2"] / (w * IN2MM / 1000) if w else None
        d["roll_w_in"] = w
    return plies, per_mat, mass


def print_sheet(path):
    plies, per_mat, mass = cut_rows(path)
    print(f"\n=== {pathlib.Path(path).stem} ===")
    print(f"{'ply':>3} {'zone':<9} {'fabric':<22} {'deg':>4} {'cut cm2':>9}")
    for p in plies:
        print(f"{p['ply']:>3} {p['zone']:<9} {p['material']:<22} "
              f"{p['angle']:>4} {p['area_mm2']/100:>9.0f}")
    print(f"\n{'fabric':<22} {'gsm':>4} {'ply':>4} {'deg':>8} {'m2':>7} "
          f"{'+nest':>7} {'roll m':>13}")
    for m, d in per_mat.items():
        ang = "/".join(str(a) for a in sorted(d["angles"]))
        roll = f"{d['roll_m']:.2f} @ {d['roll_w_in']}in" if d["roll_m"] else "-"
        print(f"{m:<22} {GSM[m]:>4.0f} {d['plies']:>4} {ang:>8} "
              f"{d['fabric_m2']:>7.3f} {d['nest_m2']:>7.3f} {roll:>13}")
    print(f"\ndry fabric mass: {mass:.0f} g  ->  ~{mass * CURED_MASS_OVER_DRY:.0f} g "
          f"cured (wet layup, ~45% resin by wt; not the thickness ratio)")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    for p in sys.argv[1:]:
        print_sheet(p)
