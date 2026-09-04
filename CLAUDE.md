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
  geometry.py     gmsh API -> planform surface -> quad8 mesh -> .inp
  layup.py        design vector -> *SHELL SECTION, COMPOSITE blocks
  deck.py         assemble complete .inp from templates/base.inp
  run.py          subprocess ccx, validate convergence, parse .dat
  sweep.py        parameter grid -> process pool -> results.parquet
  post.py         pandas/seaborn -> SVG
templates/base.inp
cases/
  smoke_cantilever/    8 elements, <5 s, closed-form large-deflection answer
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
- Fiber reference direction is set with `*ORIENTATION`, +x along the part's
  long axis.
- Ply angles are `*ORIENTATION` names, not degrees. ccx reads field 4 of a
  composite ply line as the name of an `*ORIENTATION` card, so a layup emits
  one card per distinct angle, each rotated about the shell normal. Writing
  the angle there fails with "nonexistent orientation".

  ```
  *ORIENTATION, NAME=ori_p45, SYSTEM=RECTANGULAR
  0.70710678, 0.70710678, 0., -0.70710678, 0.70710678, 0.
  *SHELL SECTION, COMPOSITE, ELSET=zone_a
  0.25, , cfrp, ori_p45
  ```
- Material is `*ELASTIC, TYPE=ENGINEERING CONSTANTS` (9 constants).
- ccx has no trailing comments. `**` must start its own line; put it after data
  on a card and the card fails to parse.
- ccx's S8R composite shell is systematically ~0.25% softer than hand CLPT,
  measured in the linear range, so it is a formulation offset and not a
  convergence error. Anything pinned tighter than ~0.5% sits inside it.

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
  `layup.py`, or `geometry.py`. It runs in under 5 seconds.
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
