"""Undeformed mesh + tip-drive / circular-arc pose plots for U-bend runs.

True FEA deformed shapes need DISP in the ``.frd`` (``*NODE FILE`` on selected
steps). Until that exists, these plots show the mid-surface mesh and the
circular-arc centerline / tip-edge targets the U-clamp path drives toward.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np


_NODE_RE = re.compile(
    r"^\s*(\d+)\s*,\s*([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s*,"
    r"\s*([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s*,"
    r"\s*([+-]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s*$"
)
_EL_RE = re.compile(r"^\s*(\d+)\s*,\s*(.+)$")


def parse_nodes_elements(inp_path: str | Path) -> tuple[dict[int, tuple[float, float, float]], list[list[int]]]:
    """Parse *NODE / *ELEMENT blocks from a deck or job.inp (S8R corners ok)."""
    lines = Path(inp_path).read_text().splitlines()
    nodes: dict[int, tuple[float, float, float]] = {}
    elements: list[list[int]] = []
    mode = None
    for line in lines:
        u = line.strip().upper()
        if u.startswith("*NODE"):
            mode = "node"
            continue
        if u.startswith("*ELEMENT"):
            mode = "elem"
            continue
        if u.startswith("*"):
            mode = None
            continue
        if mode == "node":
            m = _NODE_RE.match(line)
            if m:
                nodes[int(m.group(1))] = (
                    float(m.group(2)),
                    float(m.group(3)),
                    float(m.group(4)),
                )
        elif mode == "elem":
            m = _EL_RE.match(line)
            if m:
                ids = [int(x) for x in m.group(2).split(",") if x.strip()]
                if len(ids) >= 4:
                    elements.append(ids[:4])  # S8R corners
    if not nodes or not elements:
        raise ValueError(f"no mesh in {inp_path}")
    return nodes, elements


def circular_centerline(
    *,
    length_mm: float,
    theta_deg: float,
    s0: float,
    long_axis: str = "x",
    n: int = 64,
) -> np.ndarray:
    """Ideal circular-arc centerline in (span, z) for tip-tangent θ.

    Returns Nx3 xyz with transverse coord 0.
    """
    th = math.radians(theta_deg)
    if th <= 0:
        raise ValueError("theta must be > 0")
    r = length_mm / th
    s = np.linspace(0.0, length_mm, n)
    phi = s / r
    span = s0 + r * np.sin(phi)
    z = r * (1.0 - np.cos(phi))
    xyz = np.zeros((n, 3))
    if long_axis == "x":
        xyz[:, 0] = span
        xyz[:, 2] = z
    else:
        xyz[:, 1] = span
        xyz[:, 2] = z
    return xyz


def plot_planform(
    nodes: dict[int, tuple[float, float, float]],
    elements: list[list[int]],
    out_path: str | Path,
    *,
    title: str = "Undeformed mid-surface",
    tip_ids: Sequence[int] | None = None,
    heal_ids: Sequence[int] | None = None,
) -> Path:
    """Top view (x-y) of shell mesh; optional tip/HEAL node markers."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 3.2))
    for el in elements:
        xs = [nodes[n][0] for n in el] + [nodes[el[0]][0]]
        ys = [nodes[n][1] for n in el] + [nodes[el[0]][1]]
        ax.plot(xs, ys, color="#4a6fa5", lw=0.4, alpha=0.7)
    if heal_ids:
        hx = [nodes[n][0] for n in heal_ids if n in nodes]
        hy = [nodes[n][1] for n in heal_ids if n in nodes]
        ax.scatter(hx, hy, s=8, c="#27ae60", label="HEAL clamp", zorder=3)
    if tip_ids:
        tx = [nodes[n][0] for n in tip_ids if n in nodes]
        ty = [nodes[n][1] for n in tip_ids if n in nodes]
        ax.scatter(tx, ty, s=12, c="#c0392b", label="tip drive", zorder=3)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title(title)
    if heal_ids or tip_ids:
        ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    return out_path


def plot_bend_side(
    nodes: dict[int, tuple[float, float, float]],
    out_path: str | Path,
    *,
    length_mm: float,
    s0: float,
    long_axis: str = "x",
    angles_deg: Sequence[float] = (45.0, 90.0, 135.0, 180.0),
    title: str = "Circular-arc tip-drive poses",
) -> Path:
    """Side view (span-z): undeformed mid-chord line + circular targets."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    axis = 0 if long_axis == "x" else 1
    # mid-chord undeformed: unique span stations, mean transverse, mean z
    spans = {}
    for x, y, z in nodes.values():
        s = (x, y)[axis]
        key = round(s, 3)
        spans.setdefault(key, []).append((x, y, z))
    keys = sorted(spans)
    und_s = np.array(keys)
    und_z = np.array([np.mean([p[2] for p in spans[k]]) for k in keys])

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(und_s, und_z, color="#2c3e50", lw=1.6, label="undeformed")
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(angles_deg)))
    for deg, c in zip(angles_deg, colors):
        if deg <= 0:
            continue
        xyz = circular_centerline(
            length_mm=length_mm, theta_deg=float(deg), s0=s0, long_axis=long_axis
        )
        ax.plot(xyz[:, axis], xyz[:, 2], color=c, lw=1.4, label=f"θ={deg:g}° target")
        ax.scatter([xyz[-1, axis]], [xyz[-1, 2]], color=c, s=28, zorder=3)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("span (mm)")
    ax.set_ylabel("z (mm)")
    ax.set_title(title)
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    return out_path


def parse_nset(inp_path: str | Path, name: str) -> list[int]:
    """Collect node ids from ``*NSET, NSET=name`` (handles generate)."""
    lines = Path(inp_path).read_text().splitlines()
    want = name.lower()
    ids: list[int] = []
    mode = False
    gen = False
    for line in lines:
        u = line.strip()
        if u.upper().startswith("*NSET"):
            mode = False
            gen = "GENERATE" in u.upper()
            # NSET=foo
            m = re.search(r"NSET\s*=\s*([^\s,]+)", u, re.I)
            mode = bool(m and m.group(1).lower() == want)
            continue
        if u.startswith("*"):
            mode = False
            continue
        if not mode:
            continue
        nums = [int(x) for x in u.replace(",", " ").split() if x.strip().lstrip("-").isdigit()]
        if gen and len(nums) >= 3:
            a, b, step = nums[0], nums[1], nums[2]
            ids.extend(range(a, b + 1, step))
        else:
            ids.extend(nums)
    return ids
