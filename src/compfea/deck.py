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


def tip_u_clamp_body(
    node_u: dict[int, tuple[float, float, float]],
    *,
    energy_elset: str = "blade",
    tip_nset: str = "far_face",
    node_file: bool = False,
) -> str:
    """Step body: prescribe tip node translations; print RF + ELSE.

    ``node_file=True`` adds ``*NODE FILE`` / U so the ``.frd`` gets DISP at
    the end of this step (omit on intermediate steps to keep FRDs small).
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
    if node_file:
        lines += ["*NODE FILE", "U"]
    return "\n".join(lines)
