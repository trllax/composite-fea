#!/usr/bin/env python3
"""FIN_TEST_3 tip-weight: a dead weight hung from a blade held along its own axis.

Bench apparatus (see ``cases/fin_test_3/README.md``): the blade is clamped at the
HEAL and held so its axis is vertical, tip up. A weight hangs from the tip and
pulls straight down -- *along* the undeformed axis, in compression -- until the
blade lies over to a 90 degree tip tangent. The operator nudges the tip so it
buckles the intended way. The weight that holds 90 degrees is the bench number.

This drives the tip clip patch axially under displacement control, leaves it free
to rotate and swing, and reads the axial reaction at the HEAL as the equivalent
weight ``W``. Tip tangent ``theta`` is read from the deformed shape. Output is
``W_vs_theta.csv`` + an SVG and ``W`` interpolated at ``theta = 90 deg``.

This is a real solve; launch it detached (``setsid nohup ... &``) and poll the
run dir, per CLAUDE.md.
"""

from __future__ import annotations

import argparse
import dataclasses
import math
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from compfea.build import build
from compfea.deck import StaticStep, assemble, axial_drive_body
from compfea.geometry import Mesh
from compfea.layup import coverages_from_mesh, mesh_elsets_for_stacks
from compfea.materials import load_materials
from compfea.plybook import load_plybook, resolve
from compfea.run import parse_dat_energy, parse_dat_totals, solve
from compfea.step_mesh import mesh_step
from compfea.ubend import tip_length_mm

ROOT = Path(__file__).resolve().parents[2]
STEP = ROOT / "FIN_TEST_3.step"
LONG_AXIS = "y"
# The fin diverges even at 0.25 on the circular-arc path (cases/step_fin); the
# buckling ramp is gentler but still wants a tight cap. Re-calibrate per part.
FIN_STATIC_LINE = "0.002, 1.0, 1.E-10, 0.05"


def mesh_and_layup(*, plybook: Path, materials: Path, size_mm: float, out: Path):
    """STEP + ply book + materials -> (Mesh, Layup). Same path as run_ubend.py."""
    build(
        step=STEP,
        plybook=plybook,
        materials=materials,
        long_axis=LONG_AXIS,
        out=out,
        size_mm=size_mm,
        clamp_coverage="HEAL",
    )
    library = load_materials(materials)
    rows = load_plybook(plybook)
    raw = mesh_step(
        STEP, size_mm=size_mm, clamp_coverage="HEAL", long_axis=LONG_AXIS
    )
    coverages = coverages_from_mesh(raw.elsets)
    layup, stacks = resolve(
        rows, library, coverages, long_axis=LONG_AXIS, elsets=raw.elsets
    )
    mesh = Mesh(
        nodes=raw.nodes,
        elements=raw.elements,
        nsets=raw.nsets,
        elsets=mesh_elsets_for_stacks(stacks, all_elements=tuple(raw.elements)),
        heading=raw.heading,
    )
    return mesh, layup


def _axis_index(long_axis: str) -> int:
    if long_axis not in ("x", "y"):
        raise ValueError(f"long_axis must be 'x' or 'y', not {long_axis!r}")
    return 0 if long_axis == "x" else 1


def _station(mesh: Mesh, nset: str, axis: int) -> float:
    ids = mesh.nsets[nset]
    return sum(mesh.nodes[n][axis] for n in ids) / len(ids)


def tip_clip(
    mesh: Mesh, *, tip_nset: str = "far_face", long_axis: str = "y", n: int = 3
):
    """The ``n`` tip-edge nodes nearest the chord midpoint -- where a clip grips.

    Driving a small clip patch rather than the whole tip chord is deliberate: the
    real fixture is a string tied near mid-chord, and pinning the whole edge in
    three translations would over-constrain the section.
    """
    across = 1 - _axis_index(long_axis)
    ids = list(mesh.nsets[tip_nset])
    lo = min(mesh.nodes[i][across] for i in ids)
    hi = max(mesh.nodes[i][across] for i in ids)
    mid = 0.5 * (lo + hi)
    ids.sort(key=lambda i: abs(mesh.nodes[i][across] - mid))
    return tuple(sorted(ids[: max(1, n)]))


def _clamp_face(mesh: Mesh, axis: int) -> tuple[float, float]:
    """(s0, sign): the clamp edge nearest the tip and the tip-ward direction."""
    root = [mesh.nodes[n][axis] for n in mesh.nsets["fixed_end"]]
    tip_mid = _station(mesh, "far_face", axis)
    if tip_mid >= 0.5 * (min(root) + max(root)):
        return max(root), 1.0
    return min(root), -1.0


def _bow(mesh: Mesh, tip_offset_mm: float, *, long_axis: str, free_mm: float) -> Mesh:
    """Add a small first-mode z lift, zero at the clamp, ``tip_offset_mm`` at the tip.

    Only needed if the solve will not start: the unsymmetric shop ply book
    already couples the axial drive into bending, which usually removes the
    bifurcation without any help. ``w(s) = A (1 - cos(pi s / 2 L_free))`` with s
    measured tip-ward from the clamp face and clipped to [0, L_free], so heel
    nodes and anything past the clamp are untouched.

    This is a ``z = f(s)`` lift, not the arc-length-preserving map ``geometry.py``
    uses for real camber: it stretches the span by ~integral(w'^2)/2, which grows
    as A^2. At the sub-millimetre A used to seed a buckle that is ~1e-4 mm and
    irrelevant; do not reach for a large A here.
    """
    if tip_offset_mm == 0.0:
        return mesh
    axis = _axis_index(long_axis)
    s0, sign = _clamp_face(mesh, axis)
    nodes = {}
    for nid, xyz in mesh.nodes.items():
        s = min(max(sign * (xyz[axis] - s0), 0.0), free_mm)
        dz = tip_offset_mm * (1.0 - math.cos(math.pi * s / (2.0 * free_mm)))
        nodes[nid] = (xyz[0], xyz[1], xyz[2] + dz)
    return dataclasses.replace(mesh, nodes=nodes)


def tip_band(
    mesh: Mesh,
    *,
    tip_nset: str = "far_face",
    long_axis: str = "y",
    near_mm: float = 6.0,
    far_mm: float = 45.0,
    k: int = 5,
) -> tuple[int, ...]:
    """Up to ``k`` nodes ``near_mm``..``far_mm`` inboard of the tip, near mid-chord.

    The second station for the tip tangent: close enough that the clip->band
    segment is a real tangent (checked at runtime -- its length over the
    undeformed length must stay ~1 for an inextensible blade), far enough that
    the deflection difference carries signal.
    """
    axis = _axis_index(long_axis)
    across = 1 - axis
    _, tipward = _clamp_face(mesh, axis)
    tip_s = _station(mesh, "far_face", axis)
    edge = [mesh.nodes[i][across] for i in mesh.nsets[tip_nset]]
    mid = 0.5 * (min(edge) + max(edge))
    half = 0.35 * (max(edge) - min(edge))
    band = [
        n
        for n, xyz in mesh.nodes.items()
        if near_mm <= -tipward * (xyz[axis] - tip_s) <= far_mm
        and abs(xyz[across] - mid) <= half
    ]
    if len(band) < 2:
        raise ValueError(
            f"tip_band found {len(band)} nodes in {near_mm}..{far_mm} mm inboard "
            "of the tip; widen the window or refine the mesh"
        )
    band.sort(key=lambda n: abs(mesh.nodes[n][across] - mid))
    return tuple(sorted(band[:k]))


_DISP_HEADER = re.compile(
    r"^\s*displacements \(vx,vy,vz\) for set (?P<nset>\S+) and time\s+(?P<time>\S+)\s*$"
)


def read_set_disp(dat_path: Path, nset: str) -> dict[float, tuple[float, float, float]]:
    """Mean (ux, uy, uz) over ``nset`` per printed time, from *NODE PRINT / U."""
    want = nset.upper()
    out: dict[float, list[tuple[float, float, float]]] = {}
    current: float | None = None
    for line in dat_path.read_text().splitlines():
        header = _DISP_HEADER.match(line)
        if header:
            current = float(header["time"]) if header["nset"] == want else None
            if current is not None:
                out.setdefault(current, [])
            continue
        if current is None or not line.strip():
            continue
        fields = line.split()
        if len(fields) == 4 and fields[0].isdigit():
            out[current].append(tuple(float(v) for v in fields[1:]))
        else:
            current = None
    return {
        t: tuple(sum(v[i] for v in rows) / len(rows) for i in range(3))
        for t, rows in out.items()
        if rows
    }


def station_tangent_deg(
    near_undef: tuple[float, float, float],
    near_def: tuple[float, float, float],
    far_undef: tuple[float, float, float],
    far_def: tuple[float, float, float],
) -> tuple[float, float]:
    """Tip rotation from two stations, plus a quality ratio.

    ``near`` is the tip clip, ``far`` the inboard band. Returns
    ``(theta_deg, ratio)`` where ``theta`` is the angle between the deformed
    clip->band segment and the *undeformed* one -- so it is measured off the
    blade's own initial axis (any bow tilt or chord misalignment of the two
    stations cancels), single-valued past 90 deg via ``atan2``, and carries the
    transverse component through the 3-D cross product. ``ratio`` is the deformed
    segment length over the undeformed one: it must stay near 1 for an
    inextensible blade, and a value well off 1 means the segment spans curvature
    and ``theta`` is a chord angle, not a tangent.
    """
    du = np.asarray(near_undef, float) - np.asarray(far_undef, float)
    dd = np.asarray(near_def, float) - np.asarray(far_def, float)
    theta = math.degrees(
        math.atan2(float(np.linalg.norm(np.cross(du, dd))), float(np.dot(du, dd)))
    )
    ratio = float(np.linalg.norm(dd) / np.linalg.norm(du))
    return theta, ratio


def _centroid(mesh: Mesh, nset: str) -> np.ndarray:
    ids = mesh.nsets[nset]
    return np.mean([mesh.nodes[n] for n in ids], axis=0)


def collect(
    run_dir: Path, mesh: Mesh, *, long_axis: str, time_min: float = 0.0
) -> pd.DataFrame:
    """Per-increment W, axial shortening, energy and tip angle, from the .dat only.

    ``W`` is the **clip** RF along the drive axis -- the tip-fixture force, i.e.
    the equivalent hung weight, which does not include the blade's own body load
    (that goes to the heel). Without gravity the two are equal and opposite.

    Tip angle comes from two *NODE PRINT / U stations (deck nodes), not the .frd:
    ccx's expanded solid mesh fans one planform station into many wherever the
    shell normal is tilted, so a midsurface fit there is unreliable on a
    non-flat blade. ``time_min`` skips a gravity-settle pre-step.
    """
    axis = _axis_index(long_axis)
    comp = ("fx", "fy", "fz")[axis]
    dat = run_dir / "job.dat"

    totals = parse_dat_totals(dat)
    clip_rf = totals[totals["nset"] == "clip"].set_index("time")
    energy = parse_dat_energy(dat)
    blade = energy[energy["elset"] == "blade"].set_index("time")
    clip_u = read_set_disp(dat, "clip")
    band_u = read_set_disp(dat, "tip_band")

    clip0 = _centroid(mesh, "clip")
    band0 = _centroid(mesh, "tip_band")

    # All four print at the default (every-increment) frequency, so their times
    # line up exactly; require that rather than nearest-match so a reading is
    # never paired with another increment's.
    rows = []
    for t, cu in sorted(clip_u.items()):
        if t <= time_min + 1e-9:
            continue
        if t not in band_u or t not in clip_rf.index or t not in blade.index:
            continue
        theta, ratio = station_tangent_deg(
            tuple(clip0),
            tuple(clip0 + np.asarray(cu)),
            tuple(band0),
            tuple(band0 + np.asarray(band_u[t])),
        )
        rows.append(
            {
                "time": t,
                "u_axial_mm": abs(cu[axis]),
                "theta_deg": theta,
                "tangent_ratio": ratio,
                "W_N": abs(float(clip_rf.loc[t, comp])),
                "energy_Nmm": float(blade.loc[t, "energy"]),
            }
        )
    if not rows:
        raise RuntimeError("no increments with clip, band, RF and energy all printed")
    return pd.DataFrame(rows).sort_values("theta_deg").reset_index(drop=True)


def report(
    df: pd.DataFrame,
    run_dir: Path,
    *,
    bench_n: float | None,
    bench_deg: float | None = None,
) -> None:
    df.to_csv(run_dir / "W_vs_theta.csv", index=False)

    fig, ax = plt.subplots(figsize=(5.5, 4.0))
    ax.plot(df["theta_deg"], df["W_N"], "o-", lw=1.4, ms=4, label="model")
    ax.axvline(90.0, color="0.6", lw=1.0, ls="--")
    if bench_n is not None and bench_deg is not None:
        ax.plot([bench_deg], [bench_n], "D", color="C3", ms=8, label="bench")
    elif bench_n is not None:
        ax.axhline(bench_n, color="C3", lw=1.0, ls=":", label=f"bench ~{bench_n:g} N")
    ax.legend()
    ax.set_xlabel("tip tangent theta (deg)")
    ax.set_ylabel("clip fixture force W (N)")
    ax.set_title("FIN_TEST_3 tip-weight: W(theta)")
    fig.tight_layout()
    fig.savefig(run_dir / "W_vs_theta.svg", format="svg")
    plt.close(fig)

    if bench_n is not None and bench_deg is not None:
        lo, hi = df["theta_deg"].min(), df["theta_deg"].max()
        if lo < bench_deg < hi:
            w_at = float(np.interp(bench_deg, df["theta_deg"], df["W_N"]))
            print(
                f"W(theta={bench_deg:g} deg) = {w_at:.3f} N   "
                f"bench {bench_n:g} N   ratio {w_at / bench_n:.2f}"
            )

    worst_ratio = float((df["tangent_ratio"] - 1.0).abs().max())
    if worst_ratio > 0.03:
        print(
            f"WARNING: tip-tangent segment length is off by up to "
            f"{worst_ratio * 100:.1f}% from inextensible -- theta is a chord, "
            "not a tangent; tighten tip_band's window"
        )

    mono = df.sort_values("theta_deg")["W_N"].to_numpy()
    if np.any(np.diff(mono) < -1e-6 * np.abs(mono[:-1]).clip(min=1.0)):
        print(
            "WARNING: W(theta) is not monotonic -- a limit point breaks both the "
            "interpolation below and the 'holds stably' premise"
        )

    if (df["theta_deg"] < 90.0).any() and (df["theta_deg"] > 90.0).any():
        w90 = float(np.interp(90.0, df["theta_deg"], df["W_N"]))
        line = f"W(theta=90 deg) = {w90:.3f} N"
        if bench_n and bench_deg is None:
            line += f"   bench ~{bench_n:g} N   ratio {w90 / bench_n:.2f}"
        print(line)
    else:
        print(
            f"theta range {df['theta_deg'].min():.1f}..{df['theta_deg'].max():.1f} "
            "deg does not bracket 90; raise --uy-frac"
        )

    # Castigliano cross-check: W should equal dU/d(axial shortening).
    d = df.sort_values("u_axial_mm")
    u = d["u_axial_mm"].to_numpy()
    e = d["energy_Nmm"].to_numpy()
    if len(u) >= 3:
        dude = np.gradient(e, u)
        w = d["W_N"].to_numpy()
        dev = np.abs(dude[1:-1] - w[1:-1]) / w[1:-1]
        print(f"dU/d(delta) vs RF: max deviation {float(np.nanmax(dev)) * 100:.1f}%")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--plybook", type=Path, default=ROOT / "cases/fin_test_3/plybook_shop.csv"
    )
    p.add_argument(
        "--materials", type=Path, default=ROOT / "materials/from_shop.csv"
    )
    p.add_argument("--size-mm", type=float, default=16.0)
    p.add_argument(
        "--uy-frac",
        type=float,
        default=0.45,
        help="axial drive as a fraction of the free tip length "
        "(~0.39 L reaches theta = 90 deg for the shop plybook; 0.45 brackets it)",
    )
    p.add_argument("--clip-n", type=int, default=3)
    p.add_argument(
        "--bow-tip-mm",
        type=float,
        default=0.0,
        help="first-mode z bow amplitude (mm) if the solve will not start",
    )
    p.add_argument("--static-line", default=FIN_STATIC_LINE)
    p.add_argument(
        "--gravity",
        type=float,
        default=0.0,
        help="blade self-weight: g in m/s^2 (e.g. 9.81). Adds a settle step "
        "before the drive; W is still the clip fixture force only",
    )
    p.add_argument(
        "--bench-n",
        type=float,
        default=None,
        help="measured bench weight, N, for the plot/ratio only",
    )
    p.add_argument(
        "--bench-deg",
        type=float,
        default=None,
        help="tip angle the bench weight was recorded at; plots it as a point "
        "and prints the ratio there",
    )
    p.add_argument("--threads", type=int, default=5)
    p.add_argument("--timeout-s", type=float, default=36000.0)
    p.add_argument(
        "--run-dir",
        type=Path,
        default=ROOT / "results" / "fin_test_3_tipweight",
    )
    args = p.parse_args(argv)

    if not STEP.is_file():
        print(f"missing {STEP}", file=sys.stderr)
        return 2

    args.run_dir.mkdir(parents=True, exist_ok=True)
    mesh, layup = mesh_and_layup(
        plybook=args.plybook,
        materials=args.materials,
        size_mm=args.size_mm,
        out=args.run_dir / "build",
    )
    length = tip_length_mm(mesh, long_axis=LONG_AXIS)
    mesh = _bow(mesh, args.bow_tip_mm, long_axis=LONG_AXIS, free_mm=length)
    clip = tip_clip(mesh, long_axis=LONG_AXIS, n=args.clip_n)
    band = tip_band(mesh, long_axis=LONG_AXIS)
    mesh = dataclasses.replace(
        mesh, nsets={**mesh.nsets, "clip": clip, "tip_band": band}
    )

    axis = _axis_index(LONG_AXIS)
    _, tipward = _clamp_face(mesh, axis)
    # Compression drives the clip toward the clamp face -- against tip-ward. The
    # hung weight and gravity pull the same way, so gravity points there too.
    value = -tipward * args.uy_frac * length
    grav_dir = [0.0, 0.0, 0.0]
    grav_dir[axis] = -tipward

    steps: list[StaticStep] = []
    if args.gravity > 0.0:
        g_mm = args.gravity * 1000.0
        steps.append(
            StaticStep(
                axial_drive_body(
                    "clip",
                    axis + 1,
                    0.0,
                    drive=False,
                    gravity=(g_mm, grav_dir),
                    read_nset="fixed_end",
                    tangent_nsets=("tip_band",),
                ),
                inc=40000,
                static_line=args.static_line,
            )
        )
    steps.append(
        StaticStep(
            # A gravity settle step, if any, already ramped GRAV to full; ccx
            # carries it forward held, so the drive step does not restate it.
            axial_drive_body(
                "clip",
                axis + 1,
                value,
                read_nset="fixed_end",
                tangent_nsets=("tip_band",),
            ),
            inc=40000,
            static_line=args.static_line,
        )
    )
    final_time = float(len(steps))
    deck = assemble(
        mesh_inp=mesh.to_inp(),
        layup=layup,
        initial_bc="*BOUNDARY\nfixed_end, 1, 6",
        steps=steps,
        heading=(
            f"FIN_TEST_3 tip-weight plybook={args.plybook.name} "
            f"L={length:.3f} mm axial drive {value:.2f} mm ({args.uy_frac:g} L), "
            f"bow {args.bow_tip_mm:g} mm, gravity {args.gravity:g} m/s^2"
        ),
    )
    deck_path = args.run_dir / "deck.inp"
    deck_path.write_text(deck)
    print(
        f"deck {deck_path} elems={len(mesh.elements)} L={length:.2f} mm "
        f"clip={clip} drive_dof={axis + 1} value={value:.2f} mm "
        f"steps={len(steps)} gravity={args.gravity:g} "
        f"mats={[m.name for m in layup.materials]}"
    )

    result = solve(
        deck_path,
        args.run_dir / "ccx",
        job_name="job",
        timeout_s=args.timeout_s,
        final_time=final_time,
        threads=args.threads,
    )
    print(f"OK wall={result.wall_time_s:.1f}s increments={result.increments}")

    if args.gravity > 0.0:
        settle = parse_dat_totals(args.run_dir / "ccx" / "job.dat")
        settle = settle[
            (settle["nset"] == "fixed_end")
            & ((settle["time"] - 1.0).abs() <= 1e-6)
        ]
        if not settle.empty:
            comp = ("fx", "fy", "fz")[axis]
            print(f"blade self-weight (heel RF at settle) = "
                  f"{abs(float(settle[comp].iloc[0])):.3f} N")

    df = collect(
        args.run_dir / "ccx",
        mesh,
        long_axis=LONG_AXIS,
        time_min=1.0 if args.gravity > 0.0 else 0.0,
    )
    report(df, args.run_dir, bench_n=args.bench_n, bench_deg=args.bench_deg)
    print(f"wrote {args.run_dir / 'W_vs_theta.csv'} and .svg  ({len(df)} points)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
