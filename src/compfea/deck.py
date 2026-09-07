"""Assemble a CalculiX .inp from mesh text, layup, BCs, and steps."""

from __future__ import annotations

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
