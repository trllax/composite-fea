"""B-minimising rebuilds: interleave the +-45 plies and keep the 0deg material
count EVEN and single-material within a zone, so each zone stack is a palindrome
in BOTH angle and material -> B ~ 0.

Key lesson from this pass: a `twill_3k_198 @ 0` ply mixed with `hexcel_uni_231 @
0` in the same stack is NOT material-symmetric -> irreducible B. The clean
families are (a) all-UD at 0deg (ax3 lineage) or (b) all-twill (ax0).
"""
import csv, pathlib

TH = {"twill_3k_198": 0.2376, "hexcel_uni_231": 0.2772}
K = {"twill_3k_198": "woven", "hexcel_uni_231": "ud"}
OUT = pathlib.Path(__file__).parent
T45 = ("twill_3k_198", 45.0); T0 = ("twill_3k_198", 0.0); U0 = ("hexcel_uni_231", 0.0)

BOOKS = {
    # super-soft ~6-8 N, B~0: ALL twill_3k_198 (one fabric). +-45 skins, a 0deg
    # twill pair for span stiffness, +-45 single-ply MID/QUARTER steps.
    # root {FULL:tw45,tw0,tw0,tw45 ; MID:tw45 ; QUARTER:tw0} -> palindrome.
    "ax0": ([("FULL", *T45), ("QUARTER", *T0), ("FULL", *T0), ("MID", *T45),
             ("FULL", *T0), ("QUARTER", *T0), ("TIP", *T45), ("FULL", *T45)],
            "AX0 super-soft, single-fabric (twill_3k_198 only). +-45 skins, 0deg twill "
            "spanwise pair, +-45 MID/QUARTER steps. Palindromic stacks -> B~0."),
    # ax3 interleaved: open blade all-UD -> B~0. FULL{tw45,UD,UD,tw45} MID{tw45}
    # QUARTER{UD,tw45} TIP{tw45}. Root (QUARTER+MID+FULL) is the 7-ply palindrome
    # [45/0/45/0/45/0/45]; MID-only element is [45/0/45/0/45]. Both sym -> B = 0.
    "ax3r": ([("FULL", *T45), ("FULL", *U0), ("MID", *T45), ("QUARTER", *U0),
              ("QUARTER", *T45), ("FULL", *U0), ("TIP", *T45), ("FULL", *T45)],
             "AX3R = ax3 (all-UD at 0deg), +-45 interleaved so the root stack is the "
             "palindrome [45/0/45/0/45/0/45] -> B = 0 in every load-bearing region."),
    # ax5 interleaved -- keeps the softening tw0 ply, so expect a small residual B
    # in the open blade (tw0 vs UD material mismatch). Included for the record.
    "ax5r": ([("FULL", *T45), ("QUARTER", *U0), ("MID", *T45), ("FULL", *T0),
              ("MID", *U0), ("FULL", *U0), ("QUARTER", *T45), ("TIP", *T45), ("FULL", *T45)],
             "AX5R = ax5 interleaved. Residual B from the tw0/UD 0deg-material mix; "
             "see ax3r for the B~0 alternative at a similar flex."),
}

for name, (plies, note) in BOOKS.items():
    p = OUT / f"{name}.csv"
    with p.open("w", newline="") as f:
        f.write(f"# {note}\n")
        w = csv.writer(f); w.writerow(["ply", "zone", "material", "angle_deg", "thickness_mm", "kind"])
        for i, (z, m, a) in enumerate(plies, 1):
            w.writerow([i, z, m, f"{a:.6f}", f"{TH[m]:.4f}", K[m]])
    print(name, len(plies), "plies")
