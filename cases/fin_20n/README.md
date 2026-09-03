# fin_20n

Freediving fin — a flexible water-propulsion blade — pinned against a physical
measurement. This is the repo's only case validated against hardware rather
than a closed-form answer, which is why its tolerance is fixed.

It is also the worked example of a part-specific load case. Everything
fin-shaped lives here, not in `src/compfea/` and not in `CLAUDE.md`.

## Load case

Displacement control (see `CLAUDE.md`). The part-specific bindings are:

| generic name | this part |
| --- | --- |
| driven reference node | `tip_ref`, blade tip |
| driven DOF | 6 (rotation), to 1.5708 rad |
| reaction node set | `foot_pocket` |

Reaction force is read from the `foot_pocket` totals in the `.dat` file.

## Pinned value

TODO: record the measured reaction force, the tolerance, the measurement
conditions, and who took the measurement. A pinned number with no provenance
is not a validation case.
