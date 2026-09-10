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
    drive: bool = True,
    gravity: tuple[float, Sequence[float]] | None = None,
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

    ``read_nset`` is the fully clamped end; ``driven_nset`` RF along ``dof`` is
    the equivalent hung weight (the fixture force, not counting body load), and
    ``read_nset`` RF is the whole reaction. blade ``ELSE`` is for a Castigliano
    check, and ``tangent_nsets`` get a per-node U print each so a tip tangent
    angle can be built from two stations in the .dat -- deck nodes, not the
    .frd's expanded solid mesh.

    ``drive=False`` drops the ``*BOUNDARY`` line: a pre-step that applies only
    ``gravity`` and lets the blade settle before the drive starts. ``gravity`` is
    ``(g, (nx, ny, nz))`` -- ccx ``*DLOAD, GRAV`` on ``energy_elset``, ``g`` in
    the deck's length units per s^2 (9810 for mm), direction a unit vector. It is
    held while later steps ramp, so a self-weight run is two steps: settle under
    gravity, then drive.
    """
    if dof not in (1, 2, 3):
        raise ValueError(f"dof must be 1, 2 or 3, not {dof}")
    if node_file_frequency < 0:
        raise ValueError("node_file_frequency must be >= 0")
    if not drive and gravity is None:
        raise ValueError("drive=False needs a gravity load, else the step is empty")

    lines: list[str] = []
    if drive:
        lines += ["*BOUNDARY", f"{driven_nset}, {dof}, {dof}, {value:.10f}"]
    if gravity is not None:
        g, (nx, ny, nz) = gravity
        lines += [
            "*DLOAD",
            f"{energy_elset}, GRAV, {g:.6g}, {nx:.6g}, {ny:.6g}, {nz:.6g}",
        ]
    lines += [
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


def twist_couple_body(
    *,
    le_nset: str,
    te_nset: str,
    force: float,
    dof: int = 2,
    op_new: bool = False,
    read_nsets: Sequence[str] = ("fixed_end",),
    energy_elset: str = "blade",
) -> str:
    """Step body: a small force couple that twists the tip about its chord line.

    Appended to a solve that has already driven the blade over to a ~90 degree
    tip tangent (``axial_drive_body``). ``force`` is applied along ``dof`` (the
    blade long axis) to every node of ``le_nset`` and ``-force`` to every node of
    ``te_nset`` -- equal and opposite about the chord midline, i.e. a couple that
    rotates the tip section about its own tangent (global z at 90 degrees).

    A ``*CLOAD``, not a ``*BOUNDARY``: after the blade has rolled with an
    unsymmetric layup the tip-edge nodes sit at unknown, non-uniform axial
    positions, so a prescribed *total* displacement there cannot be formed
    without first solving. A force is incremental by construction. Run this as a
    ``+force`` / ``-force`` pair: the second step passes ``op_new=True`` so its
    ``*CLOAD`` replaces the first's rather than adding to it, and the swing
    differences out any bend-twist pre-load.

    ``force`` must be small enough that the resulting twist stays linear (a
    degree or so); the caller reads the achieved angle from the ``U`` prints and
    ``k_twist = couple / angle``.

    The held axial drive (``clip``) and root clamp (``fixed_end``) from earlier
    steps carry forward -- ccx keeps a prior ``*BOUNDARY`` unless it is restated.

    Prints ``RF`` TOTALS on ``le_nset`` / ``te_nset`` (near zero if their drive
    axis is free, as it must be for the ``*CLOAD`` to act -- a large value means
    the DOF is constrained and the couple did nothing) and on each ``read_nsets``
    (a load-path check), ``U`` on ``le_nset`` / ``te_nset`` (the achieved twist,
    and a net-translation contamination check), and ``energy_elset`` ``ELSE``.
    """
    if dof not in (1, 2, 3):
        raise ValueError(f"dof must be 1, 2 or 3, not {dof}")
    if force == 0.0:
        raise ValueError("twist_couple_body needs a non-zero force")

    cload = "*CLOAD, OP=NEW" if op_new else "*CLOAD"
    lines = [
        cload,
        f"{le_nset}, {dof}, {force:.10f}",
        f"{te_nset}, {dof}, {-force:.10f}",
    ]
    for nset in (le_nset, te_nset, *read_nsets):
        lines += [f"*NODE PRINT, NSET={nset}, TOTALS=YES", "RF"]
    for nset in (le_nset, te_nset):
        lines += [f"*NODE PRINT, NSET={nset}", "U"]
    lines += [
        f"*EL PRINT, ELSET={energy_elset}, TOTALS=ONLY",
        "ELSE",
    ]
    return "\n".join(lines)
