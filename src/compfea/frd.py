"""Streaming reader for CalculiX ``.frd`` results -- displacements only.

``run.py`` parses the ``.dat`` and that stays the rule for reactions and ELSE
energy. Nodal displacements have no ``.dat`` route on a solve already on disk,
so a deformed shape comes from here instead. Four things about the format are
load-bearing and none of them announce themselves:

- **Records are fixed width, not whitespace separated.** A negative component
  abuts the field before it::

      -1       272 0.00000E+00 0.00000E+00-3.73579E-06

  Splitting on whitespace finds four fields where there are five and reads
  ``0.00000E+00-3.73579E-06`` as one number, so every component after the first
  negative one lands on the wrong axis. Silently, with a plausible plot.
- **The field widths depend on the format code**, the trailing integer on the
  ``2C`` and ``100CL`` header lines. 0 is short (5-wide ids), 1 is long
  (10-wide), 2 is binary. Hardcoding the long widths is the same class of bug
  as ``.split()``, so the code branches on it and refuses binary.
- **``-1`` is overloaded.** The element block uses it too, with a different
  layout (``-1       318    4    0    1``). A reader keying on the prefix alone
  reads element ids and node counts as coordinates. Block state is tracked.
- **Node ids are CalculiX's expanded solid mesh, not the deck's shell nodes.**
  ccx expands a composite shell into stacked solids, one layer per ply, so a
  261-station quad8 mesh comes back as 2484 nodes with ids matching nothing in
  the ``.inp``. Never join these onto deck nodes. The ``2C`` block in the same
  file is the only consistent frame for them.

Files run 15-130 MB. ``index_disp_blocks`` makes one header-only pass and
records byte offsets; ``read_disp`` seeks to one block and reads only it.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

# End-of-step times land on integers; the tolerance only has to separate
# adjacent increments.
FRD_TIME_TOL = 1e-6

# format code -> (node id slice, three value slices)
_LAYOUT: dict[int, tuple[slice, tuple[slice, slice, slice]]] = {
    0: (slice(3, 8), (slice(8, 20), slice(20, 32), slice(32, 44))),
    1: (slice(3, 13), (slice(13, 25), slice(25, 37), slice(37, 49))),
}

Vec3 = tuple[float, float, float]


@dataclass(frozen=True)
class DispBlock:
    """One DISP block's location, so it can be read without a full rescan."""

    step: int
    time: float
    n_nodes: int
    offset: int
    fmt: int


def _layout(fmt: int) -> tuple[slice, tuple[slice, slice, slice]]:
    if fmt == 2:
        raise ValueError(
            "binary .frd (format 2) is not supported; re-run ccx without it"
        )
    try:
        return _LAYOUT[fmt]
    except KeyError:
        raise ValueError(f"unknown .frd format code {fmt}") from None


def parse_record(line: str, fmt: int = 1) -> tuple[int, Vec3]:
    """One fixed-width ``-1`` record: node id and three components.

    Never uses ``.split()`` -- see the module docstring.
    """
    node_slice, value_slices = _layout(fmt)
    return (
        int(line[node_slice]),
        (
            float(line[value_slices[0]]),
            float(line[value_slices[1]]),
            float(line[value_slices[2]]),
        ),
    )


def _trailing_int(line: str, default: int = 1) -> int:
    parts = line.split()
    try:
        return int(parts[-1])
    except (IndexError, ValueError):
        return default


def _header_time(line: str) -> float:
    """Step time from a ``100CL`` header (fixed 12-wide field)."""
    try:
        return float(line[12:24])
    except ValueError:
        return float(line.split()[2])


def index_disp_blocks(frd_path: str | Path) -> list[DispBlock]:
    """Locate every DISP block. One header-only pass; records byte offsets.

    Peak memory is one line plus a small dataclass per block, so a 125 MB file
    with a thousand blocks costs tens of kilobytes.
    """
    path = Path(frd_path)
    blocks: list[DispBlock] = []
    step = 0
    time: float | None = None
    n_nodes = 0
    fmt = 1
    with path.open() as handle:
        while True:
            line = handle.readline()
            if not line:
                break
            if line[4:6] == "2C":
                fmt = _trailing_int(line)
                continue
            if line.startswith("    1PSTEP"):
                parts = line.split()
                if len(parts) >= 4:
                    step = int(parts[3])
                continue
            if line.startswith("  100CL"):
                time = _header_time(line)
                parts = line.split()
                if len(parts) >= 4:
                    n_nodes = int(parts[3])
                fmt = _trailing_int(line, fmt)
                continue
            if line.startswith(" -4"):
                parts = line.split()
                if len(parts) > 1 and parts[1] == "DISP" and time is not None:
                    # records begin after the -5 component definition lines
                    data_offset = handle.tell()
                    while True:
                        probe = handle.tell()
                        nxt = handle.readline()
                        if not nxt or not nxt.startswith(" -5"):
                            data_offset = probe
                            handle.seek(probe)
                            break
                    blocks.append(
                        DispBlock(
                            step=step,
                            time=time,
                            n_nodes=n_nodes,
                            offset=data_offset,
                            fmt=fmt,
                        )
                    )
                continue
    return blocks


def read_nodes(frd_path: str | Path) -> dict[int, Vec3]:
    """Undeformed coordinates from the ``2C`` block, keyed by frd node id.

    These are the expanded solid nodes, not the deck's shell nodes.
    """
    path = Path(frd_path)
    nodes: dict[int, Vec3] = {}
    fmt = 1
    inside = False
    with path.open() as handle:
        for line in handle:
            if line[4:6] == "2C":
                fmt = _trailing_int(line)
                inside = True
                continue
            if not inside:
                continue
            if line.startswith(" -3"):
                break
            if line.startswith(" -1"):
                node, vec = parse_record(line, fmt)
                nodes[node] = vec
    if not nodes:
        raise ValueError(f"no nodal coordinates (2C block) in {path}")
    return nodes


def read_disp(frd_path: str | Path, block: DispBlock) -> dict[int, Vec3]:
    """Read exactly one indexed DISP block by seeking to its offset."""
    path = Path(frd_path)
    out: dict[int, Vec3] = {}
    with path.open() as handle:
        handle.seek(block.offset)
        while True:
            line = handle.readline()
            if not line or line.startswith(" -3"):
                break
            if line.startswith(" -1"):
                node, vec = parse_record(line, block.fmt)
                out[node] = vec
    return out


def disp_at_step(
    frd_path: str | Path,
    step: int,
    *,
    blocks: Sequence[DispBlock] | None = None,
) -> tuple[DispBlock, dict[int, Vec3]]:
    """The end-of-step DISP for one ``*STEP``.

    A step is solved over several increments and each writes a block, so the
    end-of-step one is picked by time. A step with no DISP raises rather than
    returning a neighbour: ``*NODE FILE`` is requested on selected steps only
    and ccx then carries it forward, so the steps before the first request have
    no shape at all, and labelling one angle's shape with another angle's
    number is exactly the silent-wrong-answer this repo is built to avoid.
    """
    index = list(blocks) if blocks is not None else index_disp_blocks(frd_path)
    if not index:
        raise LookupError(
            f"{frd_path} carries no DISP blocks; the run predates *NODE FILE. "
            "Re-solve with file_deg set to the angles you want shapes for."
        )
    candidates = [b for b in index if b.step == int(step)]
    if not candidates:
        have = sorted({b.step for b in index})
        raise LookupError(
            f"no DISP at step {step} in {frd_path}; steps with DISP: "
            f"{', '.join(str(s) for s in have)}"
        )
    block = min(candidates, key=lambda b: abs(b.time - float(b.step)))
    # Each *STATIC step runs a unit period, so the end of step N is tot-time N.
    # Without this check a step that never reached its end -- a solve killed or
    # timed out mid-step -- hands back its last converged increment, and the
    # caller plots that pose under the angle it asked for. That is the
    # "returns the last available increment" fallback CLAUDE.md forbids, and it
    # fails silently because a 176-degree curve looks entirely reasonable on a
    # plot labelled 180.
    if abs(block.time - float(block.step)) > FRD_TIME_TOL:
        raise LookupError(
            f"step {step} in {frd_path} never reached its end: the last DISP "
            f"is at tot-time {block.time:g}, not {float(step):g}. This solve "
            "stopped part way through the step; its shape is not the pose you "
            "asked for."
        )
    return block, read_disp(frd_path, block)


def deformed_midsurface(
    nodes: Mapping[int, Vec3],
    disp: Mapping[int, Vec3],
    *,
    long_axis: str = "y",
    chord_mm: float | None = None,
    round_to: int = 4,
) -> pd.DataFrame:
    """Deformed mid-surface along one chord station, root first.

    The expanded mesh stacks a node column through the thickness at each
    planform station -- one layer per ply, so a 4-ply laminate gives 9 levels.
    The mid-surface point is the **midpoint of the extreme pair**, the two ends
    of the shell director, not the mean of the whole column. The two agree only
    when the levels are symmetric *and* evenly weighted, which holds for the
    equal-thickness stacks measured so far but is not a property of the format:
    unequal ply thicknesses bias a multiplicity-weighted mean off the
    mid-surface while the director midpoint stays on it.

    ``chord_mm`` picks the transverse station; the default is the one nearest
    the middle of the part, which is the section a side view wants.

    Columns: ``span0_mm`` (undeformed span, the material ordering),
    ``span_mm``, ``z_mm``, ``transverse_mm`` (all deformed).
    """
    if long_axis not in ("x", "y"):
        raise ValueError(f"long_axis must be 'x' or 'y', not {long_axis!r}")
    axis = 0 if long_axis == "x" else 1
    across = 1 - axis

    missing = len(nodes.keys() - disp.keys())
    if missing:
        raise ValueError(
            f"{missing} of {len(nodes)} nodes have no displacement; this DISP "
            "block does not match this coordinate block"
        )

    # station -> list of (undeformed z, deformed xyz)
    stations: dict[tuple[float, float], list[tuple[float, Vec3]]] = defaultdict(list)
    for node, xyz in nodes.items():
        dx, dy, dz = disp[node]
        key = (round(xyz[axis], round_to), round(xyz[across], round_to))
        stations[key].append(
            (xyz[2], (xyz[0] + dx, xyz[1] + dy, xyz[2] + dz))
        )

    transverse = sorted({key[1] for key in stations})
    target = (
        0.5 * (transverse[0] + transverse[-1]) if chord_mm is None else float(chord_mm)
    )
    pick = min(transverse, key=lambda v: abs(v - target))

    rows = []
    for (span0, across0), column in stations.items():
        if across0 != pick:
            continue
        z_lo = min(z for z, _ in column)
        z_hi = max(z for z, _ in column)
        # Average within each extreme level, so a repeated extreme is handled.
        lo = [p for z, p in column if z == z_lo]
        hi = [p for z, p in column if z == z_hi]
        face = []
        for group in (lo, hi):
            n = len(group)
            face.append(
                (
                    sum(p[axis] for p in group) / n,
                    sum(p[2] for p in group) / n,
                    sum(p[across] for p in group) / n,
                )
            )
        rows.append(
            {
                "span0_mm": span0,
                "span_mm": 0.5 * (face[0][0] + face[1][0]),
                "z_mm": 0.5 * (face[0][1] + face[1][1]),
                "transverse_mm": 0.5 * (face[0][2] + face[1][2]),
            }
        )
    if not rows:
        raise ValueError(f"no stations at transverse coordinate {pick:g}")
    return pd.DataFrame(rows).sort_values("span0_mm").reset_index(drop=True)
