#!/usr/bin/env python3
"""20 mm-chord strip: rigid tip fixture via prescribed U on a tip band.

CalculiX forbids *RIGID BODY on shell nodes, and twin+RIGID left the shells
unloaded (ELSE=0). Instead enforce the fixture kinematically:

  - identify the last ``fixture_mm`` of span as the clamp band
  - prescribe translations on every band node so the band moves as a rigid
    body whose tip midchord follows the circular-arc tip path and whose
    orientation is tip-tangent θ (rotation about the chord axis)
  - report energy F = M/arm with M = 2U/θ, plus RF TOTALS on the band
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

from compfea.deck import StaticStep, assemble, tip_u_clamp_body
from compfea.geometry import Mesh, Outline, mesh_outline
from compfea.layup import ANSYS_EPOXY_CARBON_UD_230, Layup, Ply
from compfea.run import parse_dat_energy, solve
from compfea.ubend import force_at_theta, step_index_for, tip_length_mm

ROOT = Path(__file__).resolve().parents[2]


def _build_strip(
    *,
    chord: float,
    span: float,
    size: float,
    fixture_mm: float,
) -> tuple[Mesh, float]:
    outline = Outline.rectangle(chord=chord, span=span)
    raw = mesh_outline(
        outline,
        size=size,
        heading=f"rigid tip fixture strip chord={chord} span={span}",
    )
    tip_y = max(p[1] for p in raw.nodes.values())
    shell_fix = sorted(
        n for n, (x, y, z) in raw.nodes.items() if tip_y - y <= fixture_mm + 1e-9
    )
    nsets = dict(raw.nsets)
    nsets["tip_fixture"] = tuple(shell_fix)
    mesh = Mesh(
        nodes=raw.nodes,
        elements=raw.elements,
        nsets=nsets,
        elsets=raw.elsets,
        heading=raw.heading,
    )
    arm = tip_length_mm(mesh, tip_nset="far_face", root_nset="fixed_end", long_axis="y")
    return mesh, arm


def fixture_displacements(
    mesh: Mesh,
    theta_rad: float,
    *,
    fixture_nset: str = "tip_fixture",
    tip_nset: str = "far_face",
    long_axis: str = "y",
) -> dict[int, tuple[float, float, float]]:
    """Rigid-band U: tip midchord on circular arc, band rotated by tip tangent θ."""
    if theta_rad <= 0:
        raise ValueError("theta must be > 0")
    if long_axis != "y":
        raise ValueError("fixture_displacements currently supports long_axis='y' only")
    length = tip_length_mm(mesh, tip_nset=tip_nset, long_axis=long_axis)
    r = length / theta_rad
    s_t = r * math.sin(theta_rad)
    z_lift = r * (1.0 - math.cos(theta_rad))

    root_axis = [mesh.nodes[n][1] for n in mesh.nsets["fixed_end"]]
    tip_axis = [mesh.nodes[n][1] for n in mesh.nsets[tip_nset]]
    tip_y = 0.5 * (min(tip_axis) + max(tip_axis))
    root_hi, root_lo = max(root_axis), min(root_axis)
    if tip_y >= root_hi - 1e-9:
        s0 = root_hi
        sign = 1.0
    elif tip_y <= root_lo + 1e-9:
        s0 = root_lo
        sign = -1.0
    else:
        raise ValueError("tip overlaps clamp")
    target_s = s0 + sign * s_t

    tip_pts = [mesh.nodes[n] for n in mesh.nsets[tip_nset]]
    cx = sum(p[0] for p in tip_pts) / len(tip_pts)
    cz = sum(p[2] for p in tip_pts) / len(tip_pts)
    c, s = math.cos(theta_rad), math.sin(theta_rad)

    out: dict[int, tuple[float, float, float]] = {}
    for nid in mesh.nsets[fixture_nset]:
        x0, y0, z0 = mesh.nodes[nid]
        ox = x0 - cx
        oy = y0 - tip_y  # negative inboard when tip at high y
        oz = z0 - cz
        # rotate offset about +x by θ (span y → z)
        new_oy = oy * c - oz * s
        new_oz = oy * s + oz * c
        # tip mid goes to (cx, target_s, cz + z_lift); keep chord offset ox
        tx = cx + ox
        ty = target_s + new_oy
        tz = cz + z_lift + new_oz
        out[nid] = (tx - x0, ty - y0, tz - z0)
    return out


def _layup() -> Layup:
    ud = ANSYS_EPOXY_CARBON_UD_230
    t = 0.25
    return Layup.uniform(
        [
            Ply(t, 0.0, material=ud),
            Ply(t, 90.0, material=ud),
            Ply(t, 90.0, material=ud),
            Ply(t, 0.0, material=ud),
        ],
        long_axis="y",
        material=ud,
    )


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--chord", type=float, default=20.0)
    p.add_argument("--span", type=float, default=200.0)
    p.add_argument("--fixture-mm", type=float, default=20.0)
    p.add_argument("--size", type=float, default=8.0)
    p.add_argument("--end-deg", type=float, default=90.0)
    p.add_argument("--start-deg", type=float, default=5.0)
    p.add_argument("--step-deg", type=float, default=5.0)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument(
        "--run-dir",
        type=Path,
        default=ROOT / "results" / "tip_rigid_fixture_90",
    )
    p.add_argument("--static-line", default="0.01, 1.0, 1.E-10, 0.1")
    args = p.parse_args(argv)

    mesh, arm = _build_strip(
        chord=args.chord,
        span=args.span,
        size=args.size,
        fixture_mm=args.fixture_mm,
    )
    n = int(round((args.end_deg - args.start_deg) / args.step_deg))
    angles = [args.start_deg + i * args.step_deg for i in range(n + 1)]
    angles[-1] = args.end_deg

    steps: list[StaticStep] = []
    for deg in angles:
        th = math.radians(deg)
        node_u = fixture_displacements(mesh, th)
        body = tip_u_clamp_body(
            node_u,
            tip_nset="tip_fixture",
            energy_elset="blade",
            node_file=(deg == args.end_deg),
        )
        steps.append(StaticStep(body=body, static_line=args.static_line, inc=4000))

    deck = assemble(
        mesh_inp=mesh.to_inp(),
        layup=_layup(),
        initial_bc="*BOUNDARY\nfixed_end, 1, 6",
        steps=steps,
        heading=(
            f"prescribed-U rigid tip fixture chord={args.chord} "
            f"fixture={args.fixture_mm}mm arm={arm:.3f} "
            f"{angles[0]:g}->{angles[-1]:g}deg"
        ),
    )
    args.run_dir.mkdir(parents=True, exist_ok=True)
    deck_path = args.run_dir / "deck.inp"
    deck_path.write_text(deck)
    print(
        f"deck {deck_path} shell_elems={len(mesh.elements)} "
        f"fixture_nodes={len(mesh.nsets['tip_fixture'])} "
        f"arm={arm:.2f} mm steps={len(angles)}"
    )

    result = solve(
        deck_path,
        args.run_dir / "ccx",
        job_name="job",
        timeout_s=7200.0,
        final_time=float(len(angles)),
        threads=args.threads,
    )
    dat_path = args.run_dir / "ccx" / "job.dat"
    energy = parse_dat_energy(dat_path)
    print(f"OK wall={result.wall_time_s:.1f}s increments={result.increments}")

    for deg in (angles[0], angles[-1]) if len(angles) > 1 else angles:
        idx = step_index_for(angles, deg)
        u, m_e, f_e = force_at_theta(energy, theta_deg=deg, step_index=idx, arm_mm=arm)
        print(f"  energy θ={deg:g}°  U={u:.6g}  M={m_e:.6g}  F={f_e:.6g} N")

    # last RF totals block for tip_fixture
    dat = dat_path.read_text()
    print("  --- last tip_fixture RF totals ---")
    blocks = [m.start() for m in re.finditer(r"total force \(fx,fy,fz\) for set TIP_FIXTURE", dat)]
    if blocks:
        start = blocks[-1]
        print(dat[start : start + 200].strip())
    else:
        print("  (no TIP_FIXTURE total force block found)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
