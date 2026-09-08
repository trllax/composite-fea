"""Assemble a CalculiX .inp from mesh text, layup, BCs, and steps."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from compfea.layup import Layup


@dataclass(frozen=True)
class StaticStep:
    """One *STEP, NLGEOM static increment block."""

    body: str
    """Keyword lines inside the step after *STATIC line (BCs, loads, prints)."""

    inc: int = 8000
    static_line: str = "0.005, 1.0, 1.E-8, 0.02"
    nlgeom: bool = True

    def to_inp(self) -> str:
        step = "*STEP, NLGEOM" if self.nlgeom else "*STEP"
        if self.inc:
            step += f", INC={self.inc}"
        return "\n".join(
            [
                step,
                "*STATIC",
                self.static_line,
                self.body.rstrip(),
                "*END STEP",
            ]
        )


def assemble(
    *,
    mesh_inp: str,
    layup: Layup,
    initial_bc: str = "",
    steps: list[StaticStep] | tuple[StaticStep, ...] = (),
    heading: str = "",
) -> str:
    """Build a full deck.

    ``mesh_inp`` should contain *NODE / *ELEMENT / *NSET / *ELSET only (no
    materials or steps). ``initial_bc`` is optional pre-step *BOUNDARY (e.g.
    fixed foot). ``steps`` are appended in order.
    """
    parts: list[str] = []
    if heading:
        for line in heading.strip().splitlines():
            parts.append(line if line.startswith("**") else f"** {line}")
    parts.append(mesh_inp.strip())
    parts.append(layup.to_inp())
    if initial_bc.strip():
        parts.append(initial_bc.strip())
    for step in steps:
        parts.append(step.to_inp())
    return "\n\n".join(parts) + "\n"


#: Any value above a step's increment count leaves only the end-of-step block,
#: which is the one anyone reads. Measured: FREQUENCY=100 over 71 increments cut
#: a 4.4 MB .dat to 62 kB. Unthrottled S output is ~70x larger and every block
#: but the last is discarded in post.
STRESS_FREQUENCY = 1000000


def tip_u_clamp_body(
    node_u: dict[int, tuple[float, float, float]],
    *,
    energy_elset: str = "blade",
    tip_nset: str = "far_face",
    node_file: bool = False,
    stress: bool = False,
) -> str:
    """Step body: prescribe tip node translations; print RF + ELSE.

    ``node_file=True`` adds ``*NODE FILE`` / U so the ``.frd`` gets DISP at
    the end of this step (omit on intermediate steps to keep FRDs small).

    ``stress=True`` adds ``*EL PRINT ... GLOBAL=YES`` / S. ``GLOBAL=YES`` is
    deliberate: a bare S card reports each ply in its own ``*ORIENTATION`` frame,
    so plies at different angles arrive in different bases, while the global card
    puts every ply in one basis for the cost of a ply-angle rotation that post
    has to do anyway. Neither frame co-rotates -- see ``compfea.stress``.
    """
    lines = ["*BOUNDARY"]
    for n, (ux, uy, uz) in sorted(node_u.items()):
        lines.append(f"{n}, 1, 1, {ux:.10f}")
        lines.append(f"{n}, 2, 2, {uy:.10f}")
        lines.append(f"{n}, 3, 3, {uz:.10f}")
    lines += [
        f"*NODE PRINT, NSET={tip_nset}, TOTALS=YES",
        "RF",
        f"*EL PRINT, ELSET={energy_elset}, TOTALS=ONLY",
        "ELSE",
    ]
    if stress:
        lines += [
            f"*EL PRINT, ELSET={energy_elset}, GLOBAL=YES, "
            f"FREQUENCY={STRESS_FREQUENCY}",
            "S",
        ]
    if node_file:
        lines += ["*NODE FILE", "U"]
    return "\n".join(lines)


def axial_drive_body(
    driven_nset: str,
    dof: int,
    value: float,
    *,
    read_nset: str = "fixed_end",
    energy_elset: str = "blade",
    tangent_nsets: Sequence[str] = (),
    node_file: bool = False,
    node_file_frequency: int = 0,
) -> str:
    """Step body: ramp one translation DOF on ``driven_nset`` along the blade axis.

    Models a dead weight hung from a slender blade whose axis is held along the
    weight's line of action: the weight compresses the blade past its buckling
    load and it lies over. ``dof`` (1/2/3) is the axis of the *undeformed* blade
    and ``value`` is the axial displacement it is ramped to -- toward the clamp,
    so the caller signs it. Every other DOF on ``driven_nset`` is left free, so
    the tip rotates and swings out as it buckles: the free transverse translation
    is what keeps this a clamped-*free* column rather than a clamped-pinned one,
    which would be several times stiffer.

    Nothing here breaks the buckling symmetry -- that is the mesh's job. Give
    ``mesh_outline`` a small ``camber`` (the gate uses ``CircularCamber``) or bow
    the STEP mesh's nodes (the case runner's ``_bow``) so the column starts a
    known amount off-axis and the response is a smooth imperfect-column path with
    no bifurcation to diverge on.

    ``read_nset`` is the fully clamped end; its RF component along ``dof`` is the
    equivalent hung weight. ``driven_nset`` RF and per-node U are printed too as
    cross-checks, blade ``ELSE`` for a Castigliano check of the reaction, and
    ``tangent_nsets`` get a per-node U print each so a tip tangent angle can be
    built from two stations in the .dat -- deck nodes, not the .frd's expanded
    solid mesh.
    """
    if dof not in (1, 2, 3):
        raise ValueError(f"dof must be 1, 2 or 3, not {dof}")
    if node_file_frequency < 0:
        raise ValueError("node_file_frequency must be >= 0")

    lines = [
        "*BOUNDARY",
        f"{driven_nset}, {dof}, {dof}, {value:.10f}",
        f"*NODE PRINT, NSET={read_nset}, TOTALS=YES",
        "RF",
        f"*NODE PRINT, NSET={driven_nset}, TOTALS=YES",
        "RF",
    ]
    for nset in (driven_nset, *tangent_nsets):
        lines += [f"*NODE PRINT, NSET={nset}", "U"]
    lines += [
        f"*EL PRINT, ELSET={energy_elset}, TOTALS=ONLY",
        "ELSE",
    ]
    if node_file:
        if node_file_frequency:
            lines += [f"*NODE FILE, FREQUENCY={node_file_frequency}", "U"]
        else:
            lines += ["*NODE FILE", "U"]
    return "\n".join(lines)
