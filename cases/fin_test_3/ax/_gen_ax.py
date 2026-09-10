"""Emit asymmetric single-ply-taper ply books (the gap the symmetric grid can't
reach: a through-section between 3 and 4 UD-equivalent plies), +-45-skin flavor.

Zones (FIN_TEST_3): FULL / z_3_4ths / MID / QUARTER / HEAL / TIP. Ply 1 = -z.
Cured thicknesses come from materials/shop_inventory.csv:
  twill_3k_198 0.2376 (woven), hexcel_uni_231 0.2772 (IM UD),
  hexcel_himax_biax_100 0.12 (woven stand-in).
"""
import csv, pathlib

TH = {"twill_3k_198": 0.2376, "hexcel_uni_231": 0.2772, "hexcel_himax_biax_100": 0.12}
KIND = {"twill_3k_198": "woven", "hexcel_uni_231": "ud", "hexcel_himax_biax_100": "woven"}

OUT = pathlib.Path(r"C:\Users\trl242\projects\composite-fea\cases\fin_test_3\ax")
OUT.mkdir(exist_ok=True)


def book(name, plies, note):
    """plies: list of (zone, material, angle_deg) in -z -> +z stacking order."""
    p = OUT / f"{name}.csv"
    with p.open("w", newline="") as f:
        f.write(f"# {note}\n")
        w = csv.writer(f)
        w.writerow(["ply", "zone", "material", "angle_deg", "thickness_mm", "kind"])
        for i, (z, m, a) in enumerate(plies, 1):
            w.writerow([i, z, m, f"{a:.6f}", f"{TH[m]:.4f}", KIND[m]])
    # quick thickness-by-zone echo
    tz = {}
    for z, m, _ in plies:
        tz[z] = tz.get(z, 0) + TH[m]
    print(f"{name:10s} nplies={len(plies):2d}  " +
          "  ".join(f"{z}={tz[z]:.3f}" for z in tz))
    return p


T45 = ("twill_3k_198", 45.0)
T0 = ("twill_3k_198", 0.0)
U0 = ("hexcel_uni_231", 0.0)
B45 = ("hexcel_himax_biax_100", 45.0)

# All: +-45 twill skins (ply 1 and last), UD spanwise core, single-ply taper.
# QUARTER carries a UD + a +-45 twill (heel bending + torsion). TIP a light +-45.
COMMON_TAIL = [("QUARTER", *U0), ("QUARTER", *T45), ("TIP", *T45)]

designs = {
    # open-blade UD count / MID step / (opt) extra twill@0 in open blade
    "ax1": ([("FULL", *T45), ("FULL", *U0), ("FULL", *T0), ("MID", *U0)]
            + COMMON_TAIL + [("FULL", *T45)],
            "AX1 +-45 twill skins; open blade tw45/UD/tw0/tw45 (1 UD + 1 tw0); MID 1 UD"),
    "ax2": ([("FULL", *T45), ("FULL", *U0), ("FULL", *U0), ("MID", *U0)]
            + COMMON_TAIL + [("FULL", *T45)],
            "AX2 open blade tw45/UD/UD/tw45 (2 UD); MID 1 UD  [~17 N symmetric point, sanity]"),
    "ax3": ([("FULL", *T45), ("FULL", *U0), ("FULL", *U0), ("MID", *T45)]
            + COMMON_TAIL + [("FULL", *T45)],
            "AX3 open blade 2 UD; MID step is a single tw@45 (lighter than UD)"),
    "ax4": ([("FULL", *T45), ("FULL", *U0), ("FULL", *T0), ("FULL", *T0), ("MID", *U0)]
            + COMMON_TAIL + [("FULL", *T45)],
            "AX4 open blade tw45/UD/tw0/tw0/tw45 (1 UD + 2 tw0); MID 1 UD"),
    "ax5": ([("FULL", *T45), ("FULL", *U0), ("FULL", *T0), ("MID", *U0), ("MID", *T45)]
            + COMMON_TAIL + [("FULL", *T45)],
            "AX5 = AX1 open blade; MID step UD + tw@45 (heavier MID, more torsion at the fold)"),
    # biax analog of ax1 (skins from the purchased biax, thinner +-45 shell)
    "ax1b": ([("FULL", *B45), ("FULL", *U0), ("FULL", *T0), ("MID", *U0)]
             + [("QUARTER", *U0), ("QUARTER", *B45), ("TIP", *B45), ("FULL", *B45)],
             "AX1B biax +-45 skins; open blade bx45/UD/tw0/bx45; MID 1 UD"),
}

for n, (plies, note) in designs.items():
    book(n, plies, note)
