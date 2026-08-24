# Installation

Full installation recipe for `xenium-preprocess-pipeline`.

> **BREAKING CHANGE for existing installs.** The `scripts/submit_*.sbatch` stubs no longer
> activate the Python env by `$HOME/.local/bin` + an env NAME, and no longer hardcode the R
> module. They now REQUIRE a gitignored `scripts/env.local.conf` (resolved absolute paths) and
> refuse to run without one. **Migration: run `scripts/write-env-config.sh` once after
> upgrading** (step 6 below) — that's the whole fix.
>
> **Why:** `$HOME` and a bare env NAME can resolve differently inside a Slurm job than in the
> shell that submitted it — that's what actually broke this pipeline (twice). `sbatch
> --export=ALL` forwards env VAR *values* faithfully, but does nothing to stop `$HOME` itself
> from resolving to a different instance once the job's own shell starts — so anything that
> re-derives a path from `$HOME`, or an env name, at job time is still exposed even with
> `--export=ALL`. A file of already-resolved absolute paths on shared storage is immune to that
> class of failure; that's the only thing that changed, and it's why there's no `$HOME`-derived
> fallback — a silent fallback would just reintroduce the same bug one layer down. See
> `scripts/lib/env_config.sh` for the full rationale and the exact failure mode this replaced.

## Prerequisites

- Linux with a Slurm scheduler (`sbatch`, `--dependency=afterok:` support).
- [Micromamba](https://mamba.readthedocs.io/) (or `conda`). The sbatch
  layer (`scripts/submit_*.sbatch`) activates the Python env by its
  ABSOLUTE PREFIX, resolved once by step 6 below into
  `scripts/env.local.conf` — it does not assume micromamba lives at any
  particular path (e.g. `$HOME/.local/bin`) or that a bare env NAME
  resolves the same way inside a Slurm job as it did in your login shell.
  **`micromamba activate` (step 2 below) requires the shell hook to
  already be sourced** — a fresh install of the raw `micromamba` binary
  does not wire this up by itself. If `micromamba activate xenium` fails
  with `critical libmamba Shell not initialized` / `'micromamba' is
  running as a subprocess and can't modify the parent shell`, run this
  once (add it to your shell rc to persist across sessions):
  ```bash
  eval "$(micromamba shell hook --shell bash)"   # or --shell zsh
  ```
- [`uv`](https://docs.astral.sh/uv/) for the editable Python installs.
- R 4.4+ with Seurat, Matrix, spacexr, SPLIT, SpatialExperiment.
- Outbound network access to `github.com` / `api.github.com` for step 4 —
  `spacexr` and `SPLIT` install via `remotes::install_github()`, not CRAN.
  GitHub credentials are **not** required (both source repos are public),
  but they raise GitHub's anonymous API rate limit; see step 4 § GitHub
  access for why this matters on shared infrastructure.

## 1. Clone the repository

```bash
git clone https://github.com/settylab/xenium-preprocess-pipeline
cd xenium-preprocess-pipeline
```

## 2. Create the Python environment

```bash
scripts/create-env.sh -n xenium -f environments/xenium.yml
micromamba activate xenium
```

`scripts/create-env.sh` is a thin wrapper around `micromamba create`.
If you don't set `MAMBA_ROOT_PREFIX`, it's a plain passthrough. If you
**do** set `MAMBA_ROOT_PREFIX` to keep this install fully isolated
(e.g. off `$HOME`, for a scratch/test install, or to keep multiple
installs from sharing state) — for example:

```bash
export MAMBA_ROOT_PREFIX=/abs/path/to/isolated/root
export XDG_CACHE_HOME="$MAMBA_ROOT_PREFIX/xdg-cache"
export XDG_CONFIG_HOME="$MAMBA_ROOT_PREFIX/xdg-config"
scripts/create-env.sh -n xenium -f environments/xenium.yml
```

— the wrapper also exports `CONDA_PKGS_DIRS="$MAMBA_ROOT_PREFIX/pkgs"`
and asserts afterward that `~/.mamba/pkgs` was not touched. Without
this, micromamba's `pkgs_dirs` silently resolves to
`[$MAMBA_ROOT_PREFIX/pkgs, ~/.mamba/pkgs]` — an undocumented second
entry — so an "isolated" install can still quietly write package-cache
state to `~/.mamba/pkgs`. Nothing errors when that happens; the
wrapper's post-create check exists because this failure is otherwise
invisible. A bare `micromamba create` still works exactly as before if
you don't need isolation.

`scripts/uv-pip-install.sh` (step 3) covers uv's own package cache,
but plenty of libraries importable from this env's packages fall back
to the [XDG Base
Directory](https://specifications.freedesktop.org/basedir-spec/latest/)
spec (`$XDG_CACHE_HOME`, default `~/.cache`; `$XDG_CONFIG_HOME`,
default `~/.config`) for their OWN caches, with no repo-specific
wrapper to intercept them — e.g. `matplotlib` writes a font-list cache
(`fontlist-*.json`) and a config dir the first time it's imported,
confirmed to fire during step 5's `pytest` run. Exporting
`XDG_CACHE_HOME`/`XDG_CONFIG_HOME` once, alongside `MAMBA_ROOT_PREFIX`
above, redirects this whole class of dependency at the source instead
of chasing each library that writes to `$HOME` one at a time. If you
don't set `MAMBA_ROOT_PREFIX` (no isolation requested), leave these
unset too — everything falls back to the normal `$HOME` locations.

## 3. Install the three packages

Each package installs its own CLI entry point:

```bash
scripts/uv-pip-install.sh -e packages/xenium-preprocess
scripts/uv-pip-install.sh -e packages/ref-build
scripts/uv-pip-install.sh -e packages/rctd-split
```

`scripts/uv-pip-install.sh` is a thin wrapper around `uv pip install`.
Like `scripts/create-env.sh` (step 2), it's a plain passthrough unless
`MAMBA_ROOT_PREFIX` is set — but if it IS set, plain `uv pip install`
still writes its content-addressed cache to `~/.cache/uv` regardless
(uv doesn't read micromamba's config), so an "isolated" install can
quietly leave files behind in `$HOME` anyway. The wrapper overrides
`UV_CACHE_DIR` into `$MAMBA_ROOT_PREFIX/uv-cache` and asserts
afterward that `~/.cache/uv` was not touched — the same isolation
guarantee `create-env.sh` gives the conda side. A bare `uv pip install
-e ...` still works exactly as before if you don't need isolation.

Verify:

```bash
xenium-preprocess --help
ref-build         --help
rctd-split        --help
```

## 4. R side

### On a cluster with Lmod

**Run `ml` inside a subshell `( ... )`, not directly in your login shell.**
Lmod's `ml` sets `PATH` *and* `PYTHONPATH` for the rest of the shell
session it runs in — `PYTHONPATH` unconditionally, to point at every
Lmod-managed `python3.x/site-packages` dir the module ships (Graphviz,
SciPy-bundle, Python-bundle-PyPI, …). `micromamba activate` never
touches `PYTHONPATH`, so once it's set it stays set: even after
re-running `micromamba activate xenium`, the venv's own `pytest`
binary still imports the wrong (Lmod) `_pytest` package via the
leftover `PYTHONPATH` and crashes — and because PyYAML genuinely *is*
installed correctly in the venv, the resulting error
(`ModuleNotFoundError: No module named 'yaml'`, or an `ImportError`
importing `_pytest.config`) points nowhere near the real cause. A
subshell confines `ml`'s env changes to itself — nothing about the
outer shell's `PATH`/`PYTHONPATH` changes once the `)` closes, so step
5 sees the untouched, correctly-activated `xenium` env:

```bash
(
    ml fhR/4.4.1-foss-2023b
    scripts/install-r-packages.sh --r-lib-dir ~/R/x86_64-pc-linux-gnu-library/4.4
)
```

If you ever DO run `ml` directly (outside a subshell) and see import
errors afterward, the fix is `micromamba activate xenium && unset
PYTHONPATH` — `scripts/check-python-env.sh` (step 5) catches and
explains this case instead of leaving you to chase a misleading
`ModuleNotFoundError`.

`remotes::install_github()` checks the installed SHA across **all**
`.libPaths()` entries, not just the first — so simply prepending a
fresh library directory (an earlier version of this doc did exactly
that) does not protect the install from being silently **skipped** if
a contaminated default library (`$R_LIBS_USER`, typically
`~/R/x86_64-pc-linux-gnu-library/4.4`) already has a matching SHA. On
any shared login node where multiple users' R sessions write to a
common `$R_LIBS_USER` default, this is the common case, not the
exception — and nothing errors when it happens; you're left believing
you have an isolated install when you're actually running on
contaminated shared state. `scripts/install-r-packages.sh` sets
`R_LIBS_USER` explicitly to the target dir, strips any existing
`$HOME/R/*` entry from `.libPaths()`, passes explicit `lib=` +
`force=TRUE` to `remotes::install_github()` so the install cannot be
silently skipped, and asserts afterward that spacexr + SPLIT actually
landed in the target dir **and load successfully**.

`remotes::install_github()` also defaults `upgrade="ask"` under a
non-interactive `Rscript`, which `remotes` resolves to `upgrade="always"`
— silently upgrading every dependency it can and rebuilding it from
source, discarding fhR's pre-built versions in the process.
`scripts/install-r-packages.sh` pins `upgrade="never"` at the call site
and via `R_REMOTES_UPGRADE`, so the install only ever adds the handful
of packages fhR doesn't already ship. If fhR ships a dependency too old
for spacexr/SPLIT, the install now fails loudly instead of silently
diverging from the validated fhR stack.

### GitHub access

`remotes::install_github()` fetches both packages from GitHub's REST API
(`api.github.com`), not CRAN. Verified empirically (both source repos —
[`dmcable/spacexr`](https://github.com/dmcable/spacexr) and
[`bdsc-tds/SPLIT`](https://github.com/bdsc-tds/SPLIT) — confirmed public
via the API, `"private": false`): **GitHub credentials are not a hard
requirement.** With no `GITHUB_PAT`/`GITHUB_TOKEN` set and no git
credential store present, `scripts/install-r-packages.sh` downloads both
packages successfully over GitHub's unauthenticated API.

The real constraint is GitHub's **anonymous rate limit** — 60 requests
per source IP per hour, versus 5,000/hour authenticated — and each
`install_github()` call spends 2-3 of those. On a shared cluster where
many users' traffic egresses through the same IP, that ceiling is easy
to exhaust incidentally (any concurrent unauthenticated GitHub API
traffic on the same egress IP counts against it, not just this script),
and `remotes` then fails partway through the install with an HTTP 403
rate-limit error.

This is also why the failure is easy to miss in testing: `remotes`
silently picks up `GITHUB_PAT`/`GITHUB_TOKEN`, or credentials already
sitting in the git credential store (e.g. from a prior `gh auth login`)
— printing `Using GitHub PAT from the git credential store.` and nothing
else — so on any account with ambient credentials this step always
works, with no indication that a credential-free account would be
running against a much smaller, shared quota instead. **Recommended:**
set `GITHUB_PAT` or `GITHUB_TOKEN` before running this step, or run `gh
auth login` once (also picked up automatically), to raise the ceiling to
5,000 requests/hour.

### Off-cluster

Uncomment the R block in `environments/xenium.yml` and re-create the
env, or install R 4.4+ separately and install Seurat + Matrix + readr
+ SpatialExperiment via `install.packages()` /
`BiocManager::install()`. Then run:

```bash
scripts/install-r-packages.sh --r-lib-dir /abs/path/to/your/R/library
```

## 5. Verify the install

The `pytest` suites need the `[test]` extra — it's gated behind
`[project.optional-dependencies]` in each package's `pyproject.toml` and
is NOT part of `environments/xenium.yml`, so install it explicitly
before running bare `pytest` (omitting this step is the most common
"why does `pytest` fail with `ModuleNotFoundError: No module named
'pytest'`" report):

```bash
scripts/uv-pip-install.sh -e "packages/xenium-preprocess[test]"
scripts/uv-pip-install.sh -e "packages/ref-build[test]"
scripts/uv-pip-install.sh -e "packages/rctd-split[test]"

# Sanity check: confirms python/pytest actually resolve inside the
# `xenium` env before running anything. Catches the ml/PYTHONPATH
# shadowing described in step 4 loudly and accurately, instead of
# `pytest` failing later with a misleading ModuleNotFoundError.
scripts/check-python-env.sh --env-name xenium

# Python side
pytest packages/xenium-preprocess/tests
pytest packages/ref-build/tests
pytest packages/rctd-split/tests

# Driver dry-run
./scripts/submit_workflow.sh --dry-run \
    --sample-id     TEST \
    --output-root   /tmp/workflow_runs \
    --run-id        v0 \
    --flex-h5ad     /tmp/nowhere.h5ad \
    --celltype-marker-json /tmp/nowhere.json \
    --proseg-dir    /tmp/nowhere \
    --xenium-ranger-dir /tmp/nowhere
```

The `--dry-run` prints the sbatch invocations without submitting.

## 6. Record the resolved environment

The final install step: resolve the micromamba binary, the env
**prefix** (not name), the R library dir from step 4, and the R module
into absolute values, and write them to the gitignored
`scripts/env.local.conf`. Every `scripts/submit_*.sbatch` stub sources
this instead of re-deriving `$HOME`-relative paths inside the Slurm
job — `$HOME` can resolve to a different location inside a job than in
your login shell, which is exactly what broke this pipeline twice
before (a missing R library, and micromamba's own root-prefix
resolution) — see `scripts/lib/env_config.sh` for the full rationale.

```bash
./scripts/write-env-config.sh \
    --env-name  xenium \
    --r-lib-dir ~/R/x86_64-pc-linux-gnu-library/4.4
```

`--env-name xenium` reads back the prefix that `scripts/create-env.sh`
(step 2) already recorded to `scripts/.env-prefix-xenium` — it does not
query `micromamba env list`, whose output accumulates one row per
`xenium` env ever created on this account, across every
`MAMBA_ROOT_PREFIX` ever used (not just this install), and so cannot be
`awk`-matched by name safely once you have more than one. If you
skipped `scripts/create-env.sh` and created the env some other way (a
bare `micromamba create` / `conda create`), there is no receipt to read
back — pass `--env-prefix /abs/path/to/the/env` directly instead.

Then verify it end-to-end (also runs automatically before every
`submit_workflow.sh` dispatch — `--skip-preflight` to opt out):

```bash
./scripts/env-preflight.sh
```
