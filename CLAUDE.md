# composite-fea

Parametric large-deflection FEA of laminated composite parts. Sweep layups,
get force-displacement, pick candidates for layup.

Stack: gmsh (mesh) -> CalculiX ccx (solve) -> pandas (post). Everything is CLI,
text in and text out. No GUI in the loop.

Parts are thin laminated shells loaded well past small-deflection theory. The
worked example is a freediving fin (`cases/fin_20n`), but nothing in
`src/compfea/` is fin-specific and nothing new should be.

## Environment

macOS, Apple Silicon. Python and the solver both come from conda-forge via
mamba. Do not add uv, pip-tools, or poetry to this repo.

```
mamba env create -f environment.yml
mamba activate composite-fea
pip install -e . --no-deps        # editable install of src/compfea only
ccx -v                            # must work before anything else
```

Runtime deps are declared in both places on purpose: `environment.yml` is what
installs them, `pyproject.toml` records what the code imports. Add to both.

## Units — read this first

CalculiX has no unit system. Everything is **mm, N, tonne, MPa, s**.
Density is therefore in tonne/mm^3 (CFRP ~1.6e-9, not 1600).

A wrong unit here produces a plausible-looking number, not an error. If you
add or change a material card, state the units you used in the commit message.

## Layout

```
src/compfea/
  geometry.py     gmsh API -> planform outline -> quad8 mesh + zone ELSETs -> .inp
  layup.py        design vector -> *SHELL SECTION, COMPOSITE blocks
  deck.py         assemble complete .inp from templates/base.inp
  run.py          subprocess ccx, validate convergence, parse .dat
  sweep.py        parameter grid -> process pool -> results.parquet
  post.py         pandas/seaborn -> SVG
templates/base.inp
cases/
  smoke_cantilever/    32 elements, ~1 s, hand CLPT + closed-form elastica
  fin_20n/             freediving fin, pinned to a physical measurement
results/               gitignored
tests/
```

## Modeling conventions

- Shells are S8R. CalculiX internally expands shell elements into stacked
  solid elements, so element counts and runtimes are larger than the shell
  mesh suggests. This is expected.
- Laminates are zone-based: each `ELSET` gets its own
  `*SHELL SECTION, COMPOSITE`. There is no ply-based draping. Ply drops are
  approximated by giving inboard zones longer stacks than outboard zones.
- Fiber reference direction is set with `*ORIENTATION`. A 0-degree ply runs
  along the part's long axis, and **which global axis that is belongs to the
  case, not to this file** -- `cases/cantilever_ansys` and everything derived
  from it use +y. `layup.py` takes `long_axis` with no default for that reason;
  a default is silently wrong for half the repo. State it in the case README.
  Positive angle is counter-clockwise about +z, the standard convention. Writing
  the mirrored form instead agrees at 0, +/-45 and 90 degrees -- so a cross-ply
  or a balanced +/-45 stack cannot detect it -- and flips the sign of bend-twist
  coupling on anything unbalanced.
- Ply angles are `*ORIENTATION` names, not degrees. ccx reads field 4 of a
  composite ply line as the name of an `*ORIENTATION` card, so a layup emits
  one card per distinct angle, each rotated about the shell normal. Writing
  the angle there fails with "nonexistent orientation". Two distinct angles must
  never round onto one name: ccx does not reject duplicate NAMEs, it keeps one,
  and both plies end up at an angle nobody asked for.

  ```
  *ORIENTATION, NAME=ori_p45, SYSTEM=RECTANGULAR
  0.70710678, 0.70710678, 0., -0.70710678, 0.70710678, 0.
  *SHELL SECTION, COMPOSITE, ELSET=zone_a
  0.25, , cfrp, ori_p45
  ```
- Material is `*ELASTIC, TYPE=ENGINEERING CONSTANTS` (9 constants).
- Meshes come from `geometry.py` (gmsh). Four things there are load-bearing and
  none of them announce themselves:
  - **quad8 is gmsh element type 16**, not 10. Type 10 is the 9-node quad, and it
    is what `getElementType("Quadrangle", 2)` returns. `Mesh.SecondOrderIncomplete
    = 1` is what makes the mesh serendipity and therefore S8R.
  - **The curve-loop direction sets the element normal, and the normal decides
    which ply is at -z.** A clockwise outline gives every element a -z normal and
    inverts every unsymmetric stack in the model, silently. `geometry.py` forces
    the outline counter-clockwise and then checks every element.
  - gmsh's quad8 node order already matches CalculiX's, so connectivity passes
    through unpermuted. Asserted in `tests/test_geometry.py`, not assumed.
  - Recombination is not guaranteed to give all quads. A mixed mesh is refused,
    because only quad8 elements are read back: the triangles would be dropped,
    leaving a hole in the part and nodes attached to nothing, and ccx solves a
    deck with a hole in it without complaint.
  - Spanwise position is measured from the root edge, and `geometry.py` refuses
    an outline whose root is not the inboard end. Zone fractions and the camber
    map both depend on it, and both go wrong quietly if it is not so.
- Initial curvature is applied **after** meshing, as an arc-length-preserving map
  of the flat blade. Lifting with `z = f(y)` instead stretches the span and moves
  the stiffness as the cube of the length. Only single (developable) curvature is
  supported: a flat laminate cannot take double curvature without in-plane strain,
  and that residual stress is modelled nowhere here.
- ccx has no trailing comments. `**` must start its own line; put it after data
  on a card and the card fails to parse.
- The "ccx's S8R composite shell is ~0.25% softer than hand CLPT" note that used
  to live here was **wrong, and wrong in the expensive direction**: it explained
  away a reference error as a solver property. There are two ways to reduce a
  laminate to a beam `EI`, and on a strip with free long edges the right one is
  `b (D11 - D12^2/D22)`, not `D11 * b`. They differ by 0.247%. Against the free-
  edge value ccx lands within 0.01% on a converged mesh -- there is no
  formulation offset to budget for. See `cases/smoke_cantilever/README.md` for
  the mesh study.
- Residual disagreement with a closed form is usually **discretization**, and it
  is worth proving which before writing a tolerance around it: a formulation
  offset does not move under refinement and a mesh error does. Refine and look.

## The load case

**Displacement control, always.** Prescribe the driven DOF at a reference node
with `*BOUNDARY` and read the reaction out of a node set at the fixed end. Do
not sweep load and solve for deflection — it will not survive the softening
regime.

```
*STEP, NLGEOM, INC=2000
*STATIC
0.02, 1.0, 1.E-5, 0.05
*BOUNDARY
load_ref, 6, 6, 1.5708
*NODE PRINT, NSET=fixed_end, TOTALS=YES
RF
*END STEP
```

`NLGEOM` is mandatory. Deflections are large. Per-part node set names and the
driven DOF live in that part's case directory, not here.

## Solver rules

- `OMP_NUM_THREADS=1` per solve. Parallelism lives at the sweep level: N
  independent single-core solves, not one solve on N cores. Sparse direct
  solvers on this model size scale badly past ~4 threads.
- **Do not use PaStiX on shell models.** There are documented cases of it
  returning wrong answers and diverging on thin-shell bending where SPOOLES
  and Pardiso are correct. If you have a reason to try it, `cases/fin_20n`
  must pass first.
- The conda-forge `calculix=2.23` build links **SPOOLES only**. `SOLVER=PARDISO`
  answers "the PARDISO library is not linked", and so do PASTIX, TAUCS and SGI.
  The only alternatives here are `ITERATIVE SCALING` and `ITERATIVE CHOLESKY`,
  which measured ~140x slower on a shell bending model. Note `SOLVER=` belongs
  on `*STATIC`; put it on `*STEP` and it is silently ignored.
- `sweep.py` takes an explicit `--jobs`, defaulting to `cpu_count() - 2`.
  Never saturate the machine.

## What counts as a result

A run is valid only if **both** hold:

1. `ccx` exits zero, and
2. the last increment in the `.sta` file reached total time 1.0.

Neither condition is redundant. `ccx` skips keywords it does not recognise, so
a deck that ends up with no `*STEP` is reported as "Job finished" with exit 0
and an empty `.sta` — a clean exit status on its own proves nothing.

A partially converged solve produces a reaction force that looks fine and is
meaningless. `run.py` raises on failure — it does not return a number. Do not
add a fallback path that returns the last available increment.

Parse the `.dat` file for reactions. Do not parse `.frd`.

Reaction **moments** need a different route: `*NODE PRINT` has no `RM` label,
`*RIGID BODY` is rejected on shell nodes, and `*SECTION PRINT` silently returns
zeros for `*SHELL SECTION, COMPOSITE` (it works for a plain shell section, which
is what makes the zeros so easy to trust). Integrate ply stresses instead and
cross-check against energy -- see `cases/cantilever_ansys`.

## Regression cases — do not weaken these

- `cases/smoke_cantilever` must pass after any change to `deck.py`,
  `layup.py`, or `geometry.py`. It runs in about a second. It checks the laminate
  against hand CLPT, the large-deflection path against the closed-form elastica,
  and -- through the tip draw-in of a reversed unsymmetric stack -- that the first
  ply line really is the -z ply. Its README records what each check is blind to;
  read that before assuming a green run covers something.
- `cases/fin_20n` is pinned against a physical measurement and has a fixed
  tolerance.

IMPORTANT: if a change breaks either case, fix the change. Do not widen the
tolerance, do not update the expected value, do not mark it xfail. These two
cases are the only thing standing between a refactor that scrambles ply
ordering and fifty commits of plausible garbage.

## Long-running work

Never run a solve or a sweep in the foreground of a tool call — the Bash tool
times out well before these finish.

`sweep.py` detaches itself (`setsid nohup`), writes `results/<run_id>/status.json`
and `results/<run_id>/sweep.log`, and returns immediately with the run id.
Poll it with short commands (`jq . status.json`, `tail -n 40 sweep.log`).

Cache on a hash of the design vector so reruns are free.

## Code conventions

snake_case everywhere, including dataframe columns. Figures out as SVG.
Type hints on public functions. pytest for tests.

## How to work in this repo

Follow this loop. It exists because a converged-looking wrong answer is the
default failure mode here, and it is expensive to notice late.

1. **Plan before touching anything.** For any change bigger than a one-sentence
   diff, use plan mode and produce the plan before editing. Say which files
   change and what the check will be.
2. **Get intent before assuming physics.** Ask when the objective, material
   data, boundary conditions, or the failure criterion are ambiguous — use
   `AskUserQuestion` rather than picking a plausible value. A guessed modulus
   or an over-constrained BC produces a clean number that is wrong. Code style
   and file layout are not worth asking about; physics is.
3. **Execute to a stated, testable goal.** Every change lands with a check that
   returns pass/fail — a pytest case, a deck diff against a hand-written `.inp`,
   or a regression case. State the goal before starting and show the actual
   command output as evidence, not an assertion that it worked.
4. **Contest the result adversarially.** Before calling work done, spawn a
   subagent with fresh context to review the diff against the stated goal. It
   sees the diff and the goal, not the reasoning that produced them. Ask it for
   gaps in correctness, units, ply ordering, and convergence handling — not
   style. Fix real findings; ignore the ones that are just an agent finding
   something to say.
5. **Commit** once the checks pass and the review is clean.

Spawn subagents for exploration and verification whenever the alternative is
reading a pile of files into this conversation — mesh dumps, `.frd` contents,
sweep logs, solver source. Investigation belongs in a separate context window;
only the conclusion belongs here.

Match the subagent's model to the task, with the `model` argument on the Agent
tool. Cheap work does not need a frontier model, and physics review does not
survive a cheap one.

- **haiku** — mechanical retrieval with an unambiguous answer: grep for a
  keyword, tail a `sweep.log`, count elements in a mesh, list the node sets in
  an `.inp`, confirm a file exists. No judgement, no synthesis.
- **sonnet** — the default for exploration and summarizing: trace how a value
  flows through `src/compfea/`, triage a failed sweep across many logs, work
  out which increment a solve stalled at. Multi-file reading where the answer
  is still a matter of fact.
- **opus** — anything where being wrong is expensive and looks fine: the
  step-4 adversarial diff review, unit checks, ply ordering, convergence
  handling, `.inp` card semantics, or interpreting a result that seems
  anomalously stiff or soft. Never downgrade step 4 to save time.

When in doubt between two tiers, take the higher one for anything that touches
the deck or the physics, and the lower one for anything that only reads files.

Build order when a module does not exist yet: `run.py` -> `layup.py` and
`deck.py` -> `geometry.py` -> `sweep.py` and `post.py`. Stop after each for
review. Do not scaffold all five at once.

## Out of scope for you

Judging whether a converged result is physically correct. You will optimize
toward whatever objective `sweep.py` reports, including a layup that is stiff
because a boundary condition is over-constrained. Flag anything that looks
anomalously stiff or soft rather than reporting it as a win.
