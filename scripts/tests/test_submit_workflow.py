"""Integration tests for submit_workflow.sh.

Runs the driver against mock `sbatch` + mock package CLIs (xenium-preprocess,
ref-build, rctd-split) placed on PATH. The mock sbatch executes each
sbatch script synchronously, so the full chain runs end-to-end in the
test process — without a live cluster and without any of the three
packages installed.

Coverage (per dispatch plan §4 + comment 5251080220 §2 overwrite audit):

  * `--sample-id` + `--flex-h5ad` required; missing either exits non-zero.
  * The three jobs submit in order with `--parsable`.
  * Step 3's `--dependency` is `afterok:JOB1`; step 4's is `afterok:JOB3`.
  * RUN_ID propagation:
      - fresh run (no --run-id, no $RUN_ID)   → JOB1's fake id is used.
      - --run-id flag override                → CLI value wins.
      - RUN_ID env var                        → env value wins.
      - --run-id > $RUN_ID env                → CLI value wins.
  * `--force` guard: refuses if the run folder exists (exit 3), unless
    --force is passed (proceeds with a stderr warning).
  * End-to-end target layout: after a full chain run, the run folder
    contains every file listed in the "Final layout" table of comment
    5251080220 §1.
  * `--dry-run` submits nothing (log stays empty).
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest


DRIVER_DIR = Path(__file__).resolve().parent.parent
DRIVER = DRIVER_DIR / "submit_workflow.sh"
MOCKS_DIR = Path(__file__).resolve().parent / "mocks"


@pytest.fixture
def env(tmp_path, monkeypatch):
    """Fresh temp env for each test: an OUTPUT_ROOT under tmp_path, a
    prepped PATH with mocks in front, and log/state files for the mock
    sbatch.

    Also stubs `micromamba` + `ml` on the fixture's PATH front so the
    sbatch script's env-activation block (which shells out to micromamba
    activate + Lmod's `ml`) is a no-op during tests. Otherwise `set -e` +
    a real micromamba+broken condarc / a missing Lmod would tear the
    sbatch script down before it ever reaches the mock package CLI.
    """
    output_root = tmp_path / "runs"
    output_root.mkdir()

    log = tmp_path / "sbatch.log"
    log.touch()
    state = tmp_path / "sbatch.state"
    state.write_text("1000")

    # Stub micromamba + ml so the sbatch script's env-activation block is
    # a no-op. `micromamba shell hook --shell bash` must print nothing (so
    # the `eval "$(...)"` does nothing) and `micromamba activate` +
    # subsequent commands must exit 0. `ml` must accept any arg (`ml fhR/…`)
    # and exit 0.
    stub_bin = tmp_path / "stub_bin"
    stub_bin.mkdir()
    micromamba_stub = stub_bin / "micromamba"
    micromamba_stub.write_text(
        "#!/usr/bin/env bash\n"
        "# Test stub: `shell hook` emits nothing (empty eval), everything\n"
        "# else exits 0 without side effects.\n"
        "exit 0\n"
    )
    micromamba_stub.chmod(0o755)
    ml_stub = stub_bin / "ml"
    ml_stub.write_text("#!/usr/bin/env bash\nexit 0\n")
    ml_stub.chmod(0o755)

    original_path = os.environ.get("PATH", "")
    # Order: MOCKS (xenium-preprocess / ref-build / rctd-split / sbatch)
    # before stub_bin (micromamba, ml) before the inherited PATH (so the
    # real micromamba on the host doesn't shadow the stub).
    new_path = f"{MOCKS_DIR}:{stub_bin}:{original_path}"

    e = {
        "PATH": new_path,
        "MOCK_SBATCH_LOG": str(log),
        "MOCK_SBATCH_STATE": str(state),
        "OUTPUT_ROOT": str(output_root),
        # HOME points at tmp_path — not the parent's $HOME — because the
        # sbatch script does `export PATH="$HOME/.local/bin:$PATH"` up
        # front, which would re-prepend the parent's real micromamba
        # (from ~/.local/bin) and shadow the stub above.
        "HOME": str(tmp_path),
    }
    return {
        "env": e,
        "output_root": output_root,
        "log": log,
        "state": state,
        "tmp_path": tmp_path,
    }


def _run_driver(env, *args, extra_env=None):
    """Invoke the driver with the given args. Returns
    CompletedProcess (never raises on non-zero exit)."""
    e = dict(env["env"])
    if extra_env:
        e.update(extra_env)
    return subprocess.run(
        [str(DRIVER), *args],
        env=e,
        capture_output=True,
        text=True,
    )


def _parse_log(log_path):
    """Parse the mock sbatch log into a list of records (one per submit).
    Each record: {jobid, dependency, script, export, exit}."""
    text = log_path.read_text()
    if not text.strip():
        return []
    records = []
    cur = None
    for line in text.splitlines():
        if line == "---":
            if cur:
                records.append(cur)
            cur = {}
            continue
        if cur is None:
            continue
        key, _, val = line.partition(":")
        cur[key.strip()] = val.strip()
    if cur:
        records.append(cur)
    return records


def _flex_h5ad(tmp_path):
    p = tmp_path / "flex.h5ad"
    p.write_bytes(b"")
    return str(p)


def _marker_json(tmp_path):
    p = tmp_path / "markers.json"
    p.write_text("{}")
    return str(p)


# --------------------------------------------------------------------------
# Arg parsing
# --------------------------------------------------------------------------

def test_missing_sample_id_exits_nonzero(env, tmp_path):
    r = _run_driver(env, "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path))
    assert r.returncode != 0
    assert "--sample-id is required" in r.stderr


def test_missing_flex_h5ad_exits_nonzero(env):
    r = _run_driver(env, "--sample-id", "MH10")
    assert r.returncode != 0
    assert "--flex-h5ad is required" in r.stderr


def test_missing_celltype_marker_json_exits_nonzero(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path))
    assert r.returncode != 0
    assert "--celltype-marker-json is required" in r.stderr


def test_help_flag_exits_zero(env):
    r = _run_driver(env, "--help")
    assert r.returncode == 0
    assert "Usage: submit_workflow.sh" in r.stdout


def test_unknown_arg_exits_nonzero(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--bogus")
    assert r.returncode != 0
    assert "unknown arg" in r.stderr


# --------------------------------------------------------------------------
# Dependency chain + parsable submission
# --------------------------------------------------------------------------

def test_three_jobs_submitted_in_order(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path))
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"

    records = _parse_log(env["log"])
    assert len(records) == 3, f"expected 3 submits, got {len(records)}: {records}"

    j1, j3, j4 = records
    assert j1["script"].endswith("submit_step1.sbatch")
    assert j3["script"].endswith("submit_step3.sbatch")
    assert j4["script"].endswith("submit_step4.sbatch")


def test_step3_afterok_step1(env, tmp_path):
    _run_driver(env,
                "--sample-id", "MH10",
                "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path))
    records = _parse_log(env["log"])
    assert len(records) == 3
    j1_id = records[0]["jobid"]
    assert records[0]["dependency"] == ""
    assert records[1]["dependency"] == f"afterok:{j1_id}"


def test_step4_afterok_step3(env, tmp_path):
    _run_driver(env,
                "--sample-id", "MH10",
                "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path))
    records = _parse_log(env["log"])
    j3_id = records[1]["jobid"]
    assert records[2]["dependency"] == f"afterok:{j3_id}"


# --------------------------------------------------------------------------
# RUN_ID propagation — the load-bearing precedence rule.
# --------------------------------------------------------------------------

def _run_id_from(record):
    """Extract the RUN_ID value from an --export payload string."""
    payload = record.get("export", "")
    m = re.search(r"(?:^|,)RUN_ID=([^,]+)", payload)
    return m.group(1) if m else None


def test_run_id_defaults_to_job1_slurm_id(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path))
    assert r.returncode == 0
    records = _parse_log(env["log"])
    j1_id, _, j4 = records
    # Step 1's --export doesn't set RUN_ID (defaults inside the sbatch to
    # $SLURM_JOB_ID); step 3+4 explicitly carry RUN_ID=<JOB1 id>.
    assert _run_id_from(records[0]) is None
    assert _run_id_from(records[1]) == j1_id["jobid"]
    assert _run_id_from(records[2]) == j1_id["jobid"]


def test_cli_run_id_wins(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "cli_run")
    assert r.returncode == 0
    records = _parse_log(env["log"])
    # All three jobs receive RUN_ID=cli_run (JOB1 gets it explicitly too).
    for rec in records:
        assert _run_id_from(rec) == "cli_run", rec


def test_env_run_id_used_when_no_cli(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    extra_env={"RUN_ID": "env_run"})
    assert r.returncode == 0
    records = _parse_log(env["log"])
    for rec in records:
        assert _run_id_from(rec) == "env_run", rec


def test_cli_run_id_beats_env_run_id(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "cli_wins",
                    extra_env={"RUN_ID": "env_loses"})
    assert r.returncode == 0
    records = _parse_log(env["log"])
    for rec in records:
        assert _run_id_from(rec) == "cli_wins", rec


def _export_field(record, key):
    """Extract the value of `key` from an --export payload string."""
    payload = record.get("export", "")
    m = re.search(rf"(?:^|,){re.escape(key)}=([^,]+)", payload)
    return m.group(1) if m else None


def test_celltype_marker_json_threaded_to_all_jobs(env, tmp_path):
    marker = _marker_json(tmp_path)
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path),
                    "--celltype-marker-json", marker)
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    records = _parse_log(env["log"])
    assert len(records) == 3
    for rec in records:
        assert _export_field(rec, "CELLTYPE_MARKER_JSON") == marker, rec


# --------------------------------------------------------------------------
# --force guard
# --------------------------------------------------------------------------

def test_force_guard_refuses_when_run_folder_exists(env, tmp_path):
    (env["output_root"] / "MH10" / "MH10_myrun").mkdir(parents=True)
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "myrun")
    assert r.returncode == 3
    assert "already exists" in r.stderr
    # No sbatch submissions on refusal.
    assert _parse_log(env["log"]) == []


def test_force_guard_proceeds_with_force(env, tmp_path):
    # --force semantics (post-comment-5260289249): DELETE the existing run
    # folder, then submit the full chain fresh. Meant for a from-scratch
    # re-run under an already-used run-id.
    run_dir = env["output_root"] / "MH10" / "MH10_myrun"
    run_dir.mkdir(parents=True)
    # Leave a marker file — --force must remove it (destructive semantics).
    (run_dir / "stale.txt").write_text("stale")
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "myrun",
                    "--force")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    assert "removing existing run folder" in r.stderr
    assert not (run_dir / "stale.txt").exists(), "--force did not delete pre-existing content"
    records = _parse_log(env["log"])
    assert len(records) == 3


def test_force_guard_silent_on_fresh_folder(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "fresh_run")
    assert r.returncode == 0
    # No warning about existing folder.
    assert "already exists" not in r.stderr
    assert "re-running under --force" not in r.stderr


def test_force_guard_only_fires_when_run_id_bound_up_front(env, tmp_path):
    # Even with a stale folder for SOME run-id, if the driver picks its
    # own RUN_ID (= JOB1's SLURM_JOB_ID), the guard is not applicable —
    # the fresh id names a fresh folder by construction.
    (env["output_root"] / "MH10" / "MH10_1001").mkdir(parents=True)
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path))
    # Passes: driver never checked for the folder because RUN_ID wasn't
    # bound before submission.
    assert r.returncode == 0, r.stderr


# --------------------------------------------------------------------------
# End-to-end target layout
# --------------------------------------------------------------------------

def test_end_to_end_layout(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "e2e_test")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"

    run_dir = env["output_root"] / "MH10" / "MH10_e2e_test"

    # Bytes-locked layout from comment 5251080220 §1.
    expected = [
        "spatial_adata/MH10_xenium_ranger.h5ad",
        "spatial_adata/MH10_proseg_raw.h5ad",
        "spatial_adata/MH10_proseg_purified.h5ad",
        "rctd/MH10_reference_post_rules.h5ad",
        "rctd/MH10_test_object.rds",
        "rctd/MH10_reference.rds",
        "rctd/MH10_rctd_results.rds",
        "config.yaml",
        "logs/xenium-preprocess.log",
        "logs/ref-build.log",
        "logs/rctd-split.log",
    ]
    missing = [rel for rel in expected if not (run_dir / rel).exists()]
    assert not missing, f"missing paths under {run_dir}: {missing}"


def test_resolved_config_merged_across_steps(env, tmp_path):
    _run_driver(env,
                "--sample-id", "MH10",
                "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                "--run-id", "e2e_merge")
    cfg = (env["output_root"] / "MH10" / "MH10_e2e_merge"
           / "config.yaml").read_text()
    # All three top-level keys survive (invariant enforced by the shared
    # _merge_config helper in each pipeline — the mocks here emulate it).
    assert "step1:" in cfg
    assert "step3:" in cfg
    assert "step4:" in cfg
    # The flex path was recorded verbatim.
    assert _flex_h5ad(tmp_path) in cfg or "flex_h5ad_path" in cfg


# --------------------------------------------------------------------------
# --dry-run
# --------------------------------------------------------------------------

def test_dry_run_submits_nothing(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--dry-run")
    assert r.returncode == 0
    # Mock sbatch is bypassed under --dry-run, so the log stays empty.
    assert _parse_log(env["log"]) == []
    # But the driver still prints the chain summary.
    assert "Workflow chain submitted" in r.stdout


# --------------------------------------------------------------------------
# --donor-h5ad / --fallback-donor-h5ad — Tracy's per-run inputs
# (settylab/TracyY123-nexus#26 comment 5257935882). Optional, repeatable;
# threaded to step 3 as numbered env vars, expanded in submit_step3.sbatch
# into --donor-h5ad / --fallback-donor-h5ad flags on `ref-build run`.
# --------------------------------------------------------------------------

def _step3_log(env, sample, run_id):
    """Return the mock ref-build's ref-build.log contents (donors + fallbacks
    recorded, one per line)."""
    p = env["output_root"] / sample / f"{sample}_{run_id}" / "logs" / "ref-build.log"
    return p.read_text() if p.exists() else ""


def test_donor_h5ad_empty_by_default(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "no_donors")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    records = _parse_log(env["log"])
    # COUNT=0 threaded to every job (so downstream can inspect the origin).
    for rec in records:
        assert _export_field(rec, "DONOR_H5AD_COUNT") == "0", rec
        assert _export_field(rec, "FALLBACK_H5AD_COUNT") == "0", rec
    # And ref-build received no donor / fallback flags.
    log = _step3_log(env, "MH10", "no_donors")
    assert "donors_count=0" in log
    assert "fallback_donors_count=0" in log


def test_single_donor_h5ad_threaded_to_ref_build(env, tmp_path):
    donor = str(tmp_path / "donor.h5ad")
    Path(donor).write_bytes(b"")
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--donor-h5ad", donor,
                    "--run-id", "one_donor")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    records = _parse_log(env["log"])
    for rec in records:
        assert _export_field(rec, "DONOR_H5AD_COUNT") == "1", rec
        assert _export_field(rec, "DONOR_H5AD_1") == donor, rec
    log = _step3_log(env, "MH10", "one_donor")
    assert "donors_count=1" in log
    assert f"donor={donor}" in log


def test_multiple_donor_h5ads_threaded_repeatable(env, tmp_path):
    donors = []
    for i in range(3):
        p = tmp_path / f"donor_{i}.h5ad"
        p.write_bytes(b"")
        donors.append(str(p))
    extra_args = []
    for d in donors:
        extra_args += ["--donor-h5ad", d]
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    *extra_args,
                    "--run-id", "multi_donor")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    records = _parse_log(env["log"])
    for rec in records:
        assert _export_field(rec, "DONOR_H5AD_COUNT") == "3", rec
        for i, d in enumerate(donors, 1):
            assert _export_field(rec, f"DONOR_H5AD_{i}") == d, (i, rec)
    log = _step3_log(env, "MH10", "multi_donor")
    assert "donors_count=3" in log
    for d in donors:
        assert f"donor={d}" in log


def test_fallback_donor_h5ad_threaded(env, tmp_path):
    fb = str(tmp_path / "fallback.h5ad")
    Path(fb).write_bytes(b"")
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--fallback-donor-h5ad", fb,
                    "--run-id", "one_fb")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    records = _parse_log(env["log"])
    for rec in records:
        assert _export_field(rec, "FALLBACK_H5AD_COUNT") == "1", rec
        assert _export_field(rec, "FALLBACK_H5AD_1") == fb, rec
    log = _step3_log(env, "MH10", "one_fb")
    assert "fallback_donors_count=1" in log
    assert f"fallback_donor={fb}" in log


def test_donor_and_fallback_donor_h5ads_together(env, tmp_path):
    """The MH3-style invocation: 2 donors + 1 fallback."""
    donor_a = str(tmp_path / "donorA.h5ad"); Path(donor_a).write_bytes(b"")
    donor_b = str(tmp_path / "donorB.h5ad"); Path(donor_b).write_bytes(b"")
    fb      = str(tmp_path / "fallback.h5ad"); Path(fb).write_bytes(b"")
    r = _run_driver(env,
                    "--sample-id", "MH3",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--donor-h5ad", donor_a,
                    "--donor-h5ad", donor_b,
                    "--fallback-donor-h5ad", fb,
                    "--run-id", "mh3_e2e")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    log = _step3_log(env, "MH3", "mh3_e2e")
    assert "donors_count=2" in log
    assert f"donor={donor_a}" in log
    assert f"donor={donor_b}" in log
    assert "fallback_donors_count=1" in log
    assert f"fallback_donor={fb}" in log


# --------------------------------------------------------------------------
# --celltype-col-for-ref-build — Tracy's step-3 celltype-column override
# (settylab/TracyY123-nexus#26 comment 5258483286). Optional; threaded via
# CELLTYPE_COL_FOR_REF_BUILD env var; expanded in submit_step3.sbatch into
# `ref-build run --celltype-col <name>`. Unset ⇒ ref-build's config default.
# --------------------------------------------------------------------------

def test_celltype_col_for_ref_build_unset_by_default(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "no_ct_col")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    records = _parse_log(env["log"])
    # Not threaded at all when the flag is omitted (step 3 falls through
    # to ref-build's config default).
    for rec in records:
        assert _export_field(rec, "CELLTYPE_COL_FOR_REF_BUILD") is None, rec
    log = _step3_log(env, "MH10", "no_ct_col")
    # Mock ref-build records an empty celltype_col when --celltype-col
    # wasn't passed by submit_step3.sbatch.
    assert "celltype_col=" in log
    assert "celltype_col=refined" not in log


def test_celltype_col_for_ref_build_threaded_to_step3(env, tmp_path):
    col = "refined_celltype_update_lymphocyes_new"
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--celltype-col-for-ref-build", col,
                    "--run-id", "ct_col_set")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    records = _parse_log(env["log"])
    # Threaded on every job's --export payload (common exports).
    for rec in records:
        assert _export_field(rec, "CELLTYPE_COL_FOR_REF_BUILD") == col, rec
    # And expanded into `ref-build run --celltype-col <col>` at step 3.
    log = _step3_log(env, "MH10", "ct_col_set")
    assert f"celltype_col={col}" in log


# --------------------------------------------------------------------------
# Slurm log routing (settylab/TracyY123-nexus#26 comments 5259180881 +
# 5274187257). Route per-step stdout/stderr into the run-scoped logs/ folder
# so the run directory is self-contained, and name the log by the pipeline
# package invoked (xenium-preprocess / ref-build / rctd-split) rather than
# the internal "stepN" label. Verified via the mock sbatch's `output:` field.
# --------------------------------------------------------------------------

# Log-file suffix per step, in submission order (step 1, 3, 4).
_STAGE_LOG_SUFFIXES = ("xenium-preprocess", "ref-build", "rctd-split")


def test_slurm_output_routed_to_run_logs_with_run_id(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path),
                    "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "log_relocate")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    records = _parse_log(env["log"])
    assert len(records) == 3

    run_logs = env["output_root"] / "MH10" / "MH10_log_relocate" / "logs"
    for suffix, rec in zip(_STAGE_LOG_SUFFIXES, records):
        assert rec["output"] == str(run_logs / f"slurm-%j-{suffix}.log"), rec
    # And the logs directory was created before submission.
    assert run_logs.is_dir()


def test_slurm_output_uses_jobid_when_run_id_auto(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path),
                    "--celltype-marker-json", _marker_json(tmp_path))
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    records = _parse_log(env["log"])
    j1_id = records[0]["jobid"]

    # JOB1's --output uses sbatch's %j placeholder in BOTH the run-folder
    # slot (SAMPLE_<jobid>) and the filename slot — expanded on the compute
    # node when the log file is opened.
    assert records[0]["output"] == str(
        env["output_root"] / "MH10" / "MH10_%j" / "logs" / "slurm-%j-xenium-preprocess.log"
    ), records[0]

    # JOB3 and JOB4 (RUN_ID bound to JOB1's parsable id) get the explicit path.
    run_logs = env["output_root"] / "MH10" / f"MH10_{j1_id}" / "logs"
    assert records[1]["output"] == str(run_logs / "slurm-%j-ref-build.log"), records[1]
    assert records[2]["output"] == str(run_logs / "slurm-%j-rctd-split.log"), records[2]
    # The parent <sample> dir was created up front; the run-scoped logs/
    # dir was created immediately after JOB1's id came back (before JOB3/4
    # submit).
    assert run_logs.is_dir()


def test_workflow_submit_log_mirrors_summary(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path),
                    "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "log_mirror")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"

    wf_log = (env["output_root"] / "MH10" / "MH10_log_mirror"
              / "logs" / "workflow-submit.log")
    assert wf_log.exists(), f"missing workflow-submit.log at {wf_log}"
    text = wf_log.read_text()
    assert "Workflow chain submitted" in text
    assert "sample_id  = MH10" in text
    assert "run_id     = log_mirror" in text
    # Job ids show up on the step lines.
    records = _parse_log(env["log"])
    for step, rec in zip((1, 3, 4), records):
        assert rec["jobid"] in text, (step, rec)


def test_dry_run_writes_no_workflow_submit_log(env, tmp_path):
    _run_driver(env,
                "--sample-id", "MH10",
                "--flex-h5ad", _flex_h5ad(tmp_path),
                "--celltype-marker-json", _marker_json(tmp_path),
                "--run-id", "dr",
                "--dry-run")
    # No filesystem side effects under --dry-run (no logs/ dir, no summary log).
    assert not (env["output_root"] / "MH10" / "MH10_dr" / "logs").exists()


# --------------------------------------------------------------------------
# --reuse-run-dir — resume-friendly counterpart to --force. Keeps the run
# folder intact (each step overwrites the files it writes; other files
# preserved). settylab/TracyY123-nexus#26 comment 5260289249.
# --------------------------------------------------------------------------

def _seed_step1_outputs(env, sample, run_id):
    """Pretend step 1 has already run: create the artifacts a resume path
    would need to find in the run folder."""
    run_dir = env["output_root"] / sample / f"{sample}_{run_id}"
    (run_dir / "spatial_adata").mkdir(parents=True, exist_ok=True)
    (run_dir / "rctd").mkdir(parents=True, exist_ok=True)
    (run_dir / "logs").mkdir(parents=True, exist_ok=True)
    (run_dir / "spatial_adata" / f"{sample}_proseg_raw.h5ad").write_bytes(b"step1")
    (run_dir / "spatial_adata" / f"{sample}_xenium_ranger.h5ad").write_bytes(b"step1")
    (run_dir / "rctd" / f"{sample}_test_object.rds").write_bytes(b"step1")
    return run_dir


def _seed_step3_outputs(env, sample, run_id):
    """Seed step-1 + step-3 outputs, as if the chain reached the end of step 3."""
    run_dir = _seed_step1_outputs(env, sample, run_id)
    (run_dir / "rctd" / f"{sample}_reference.rds").write_bytes(b"step3")
    (run_dir / "rctd" / f"{sample}_reference_post_rules.h5ad").write_bytes(b"step3")
    return run_dir


def test_reuse_run_dir_keeps_existing_folder(env, tmp_path):
    run_dir = _seed_step1_outputs(env, "MH10", "reuse")
    marker = run_dir / "spatial_adata" / "MH10_proseg_raw.h5ad"
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "reuse",
                    "--reuse-run-dir")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    assert "reuse-run-dir" in r.stderr or "keeping existing run folder" in r.stderr
    # Pre-existing file survived — --reuse-run-dir does NOT wipe.
    assert marker.exists()


def test_reuse_run_dir_and_force_are_mutually_exclusive(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "conflict",
                    "--reuse-run-dir",
                    "--force")
    assert r.returncode != 0
    assert "mutually exclusive" in r.stderr


def test_reuse_run_dir_on_fresh_folder_is_noop(env, tmp_path):
    # No existing folder — --reuse-run-dir shouldn't error, just proceed.
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "fresh_reuse",
                    "--reuse-run-dir")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    records = _parse_log(env["log"])
    assert len(records) == 3


# --------------------------------------------------------------------------
# --start-step — skip earlier steps and resume mid-chain. Tracy's use case
# from comment 5260289249: step 1 already ran, wants to redo step 3 + 4
# under the same run-id without wiping step 1's outputs.
# --------------------------------------------------------------------------

def test_start_step_default_is_1(env, tmp_path):
    # Explicit --start-step 1 must be equivalent to omitting the flag.
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--start-step", "1",
                    "--run-id", "explicit1")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    records = _parse_log(env["log"])
    assert len(records) == 3
    assert records[0]["script"].endswith("submit_step1.sbatch")


def test_start_step_invalid_value_rejected(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--start-step", "2",
                    "--run-id", "bad")
    assert r.returncode != 0
    assert "one of 1|xenium-preprocess, 3|ref-build, 4|rctd-split" in r.stderr


def test_start_step_3_requires_run_id(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--start-step", "3")
    assert r.returncode != 0
    assert "requires --run-id" in r.stderr


def test_start_step_3_missing_run_folder_fails_loud(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--start-step", "3",
                    "--run-id", "no_such_run")
    assert r.returncode != 0
    assert "does not exist" in r.stderr or "nothing to resume" in r.stderr


def test_start_step_3_missing_prereqs_fails_loud(env, tmp_path):
    # Folder exists but step-1 outputs are absent — fail-loud with the
    # missing paths listed.
    (env["output_root"] / "MH10" / "MH10_partial").mkdir(parents=True)
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--start-step", "3",
                    "--run-id", "partial")
    assert r.returncode != 0
    assert "requires these prior outputs" in r.stderr
    assert "test_object.rds" in r.stderr


def test_start_step_3_submits_step3_and_step4_only(env, tmp_path):
    _seed_step1_outputs(env, "MH10", "resume3")
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--start-step", "3",
                    "--run-id", "resume3")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    records = _parse_log(env["log"])
    assert len(records) == 2, f"expected 2 submits (step 3 + step 4), got {len(records)}"
    assert records[0]["script"].endswith("submit_step3.sbatch")
    assert records[1]["script"].endswith("submit_step4.sbatch")
    # Step 3 has NO dependency (step 1 was skipped).
    assert records[0]["dependency"] == "", records[0]
    # Step 4 chains off step 3's jobid.
    assert records[1]["dependency"] == f"afterok:{records[0]['jobid']}", records[1]


def _write_config(tmp_path, body: str) -> str:
    p = tmp_path / "workflow_config.yaml"
    p.write_text(body)
    return str(p)


def test_workflow_config_seeds_driver_defaults(env, tmp_path):
    # --config <path> populates sample_id / flex_h5ad / celltype-marker-json
    # + a per-step named param so the chain runs without repeating them on
    # the CLI. Assert (1) the run submits, (2) the ref_build named param
    # threads through to the sbatch --export payload.
    flex = _flex_h5ad(tmp_path)
    marker = _marker_json(tmp_path)
    cfg = _write_config(tmp_path, textwrap.dedent(f"""\
        sample_id: MH10
        flex_h5ad: {flex}
        celltype_marker_json: {marker}
        ref_build:
          donor_borrow_cap: 137
    """))
    r = _run_driver(env, "--config", cfg)
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    records = _parse_log(env["log"])
    assert len(records) == 3
    step3 = next(r for r in records if r["script"].endswith("submit_step3.sbatch"))
    assert "STEP3_DONOR_BORROW_CAP=137" in step3["export"]


def test_workflow_config_cli_wins_over_yaml(env, tmp_path):
    # CLI flags override YAML defaults — precedence rule Tracy insisted on.
    flex = _flex_h5ad(tmp_path)
    marker = _marker_json(tmp_path)
    cfg = _write_config(tmp_path, textwrap.dedent(f"""\
        sample_id: WRONG_FROM_YAML
        flex_h5ad: {flex}
        celltype_marker_json: {marker}
        ref_build:
          donor_borrow_cap: 999
    """))
    r = _run_driver(env, "--config", cfg,
                    "--sample-id", "MH10",
                    "--step3-donor-borrow-cap", "42")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    records = _parse_log(env["log"])
    step1 = next(r for r in records if r["script"].endswith("submit_step1.sbatch"))
    # CLI --sample-id wins over YAML sample_id.
    assert "SAMPLE=MH10" in step1["export"]
    step3 = next(r for r in records if r["script"].endswith("submit_step3.sbatch"))
    # CLI --step3-donor-borrow-cap wins over YAML ref_build.donor_borrow_cap.
    assert "STEP3_DONOR_BORROW_CAP=42" in step3["export"]


def test_workflow_config_missing_file_fails_loud(env, tmp_path):
    r = _run_driver(env, "--config", str(tmp_path / "no_such.yaml"),
                    "--sample-id", "MH10")
    assert r.returncode != 0
    assert "not found" in r.stderr


def test_workflow_config_rejects_numeric_step_keys(env, tmp_path):
    # settylab/TracyY123-nexus#26 comment 5321822161: the config schema
    # is semantic-only. Numeric `step1:` / `step3:` / `step4:` at the top
    # level must fail loud so a typo (or a copy-paste from an older
    # sketch) doesn't silently no-op the whole section.
    flex = _flex_h5ad(tmp_path)
    marker = _marker_json(tmp_path)
    for legacy in ("step1", "step3", "step4"):
        cfg = _write_config(tmp_path, textwrap.dedent(f"""\
            sample_id: MH10
            flex_h5ad: {flex}
            celltype_marker_json: {marker}
            {legacy}:
              donor_borrow_cap: 100
        """))
        r = _run_driver(env, "--config", cfg)
        assert r.returncode != 0, (
            f"expected --config to reject `{legacy}:` key; stderr:\n{r.stderr}"
        )
        # Error must name both the offending key and the accepted form.
        assert legacy in r.stderr, r.stderr
        assert "xenium_preprocess" in r.stderr or "ref_build" in r.stderr \
            or "rctd_split" in r.stderr, r.stderr


def test_workflow_config_accepts_all_three_semantic_step_keys(env, tmp_path):
    # Round-trip check: xenium_preprocess / ref_build / rctd_split all
    # thread through to the correct STEP{1,3,4}_* env var on the sbatch
    # --export payload.
    flex = _flex_h5ad(tmp_path)
    marker = _marker_json(tmp_path)
    cfg = _write_config(tmp_path, textwrap.dedent(f"""\
        sample_id: MH10
        flex_h5ad: {flex}
        celltype_marker_json: {marker}
        xenium_preprocess:
          qc_min_counts_cell: 11
        ref_build:
          donor_borrow_cap: 33
        rctd_split:
          umi_min: 55
    """))
    r = _run_driver(env, "--config", cfg)
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    records = _parse_log(env["log"])
    step1 = next(r for r in records if r["script"].endswith("submit_step1.sbatch"))
    step3 = next(r for r in records if r["script"].endswith("submit_step3.sbatch"))
    step4 = next(r for r in records if r["script"].endswith("submit_step4.sbatch"))
    assert "STEP1_QC_MIN_COUNTS_CELL=11" in step1["export"]
    assert "STEP3_DONOR_BORROW_CAP=33"    in step3["export"]
    assert "STEP4_UMI_MIN=55"             in step4["export"]


def test_start_step_semantic_alias_ref_build(env, tmp_path):
    # `ref-build` is the semantic alias for numeric `3`; the two must
    # take identical paths through the driver.
    _seed_step1_outputs(env, "MH10", "resume_semantic")
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--start-step", "ref-build",
                    "--run-id", "resume_semantic")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    records = _parse_log(env["log"])
    assert len(records) == 2, f"expected 2 submits (step 3 + step 4), got {len(records)}"
    assert records[0]["script"].endswith("submit_step3.sbatch")
    assert records[1]["script"].endswith("submit_step4.sbatch")


def test_start_step_semantic_alias_rctd_split(env, tmp_path):
    # `rctd-split` is the semantic alias for numeric `4`.
    _seed_step1_outputs(env, "MH10", "resume4_semantic")
    _seed_step3_outputs(env, "MH10", "resume4_semantic")
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--start-step", "rctd-split",
                    "--run-id", "resume4_semantic")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    records = _parse_log(env["log"])
    assert len(records) == 1, f"expected 1 submit (step 4 only), got {len(records)}"
    assert records[0]["script"].endswith("submit_step4.sbatch")


def test_start_step_semantic_alias_xenium_preprocess(env, tmp_path):
    # `xenium-preprocess` is the semantic alias for numeric `1` (full chain).
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--start-step", "xenium-preprocess",
                    "--run-id", "start1_semantic")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    records = _parse_log(env["log"])
    assert len(records) == 3


def test_start_step_3_implies_reuse_run_dir(env, tmp_path):
    _seed_step1_outputs(env, "MH10", "implicit_reuse")
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--start-step", "3",
                    "--run-id", "implicit_reuse")
    # No --reuse-run-dir on the CLI, yet the driver proceeds — proves
    # --start-step 3 implies --reuse-run-dir.
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"


def test_start_step_3_incompatible_with_force(env, tmp_path):
    _seed_step1_outputs(env, "MH10", "step3_force")
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--start-step", "3",
                    "--run-id", "step3_force",
                    "--force")
    assert r.returncode != 0
    assert "incompatible" in r.stderr


def test_start_step_4_missing_step3_outputs_fails_loud(env, tmp_path):
    # Step-1 outputs present, but step-3 reference.rds absent.
    _seed_step1_outputs(env, "MH10", "step4_missing_ref")
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--start-step", "4",
                    "--run-id", "step4_missing_ref")
    assert r.returncode != 0
    assert "reference.rds" in r.stderr


def test_start_step_4_submits_only_step4(env, tmp_path):
    _seed_step3_outputs(env, "MH10", "resume4")
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--start-step", "4",
                    "--run-id", "resume4")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    records = _parse_log(env["log"])
    assert len(records) == 1, f"expected 1 submit (step 4), got {len(records)}"
    assert records[0]["script"].endswith("submit_step4.sbatch")
    assert records[0]["dependency"] == ""


def test_start_step_summary_marks_skipped_steps(env, tmp_path):
    _seed_step1_outputs(env, "MH10", "sum3")
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path), "--celltype-marker-json", _marker_json(tmp_path),
                    "--start-step", "3",
                    "--run-id", "sum3")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    # Summary shows step 1 as skipped, step 3 / step 4 as submitted jobids.
    assert "step 1     = skipped" in r.stdout
    assert "start_step = 3" in r.stdout


# --------------------------------------------------------------------------
# Per-step parameter exposure — Tracy's request on
# settylab/TracyY123-nexus#26 comment 5277569725.
#
# Every step's config surface must be regulable at the submit_workflow.sh
# level via three layers (in precedence order, last wins):
#   (2) --stepN-config <path>
#   (3) --override stepN.<dotted.key>=<yaml-val>
#   (4) --stepN-<param> <value>            (named flag)
#
# The tests exercise each layer + the base64 encoding roundtrip that
# ferries layers 2+3 through slurm's --export=ALL,K=V payload.
# --------------------------------------------------------------------------

def _step1_log(env, sample, run_id):
    p = env["output_root"] / sample / f"{sample}_{run_id}" / "logs" / "xenium-preprocess.log"
    return p.read_text() if p.exists() else ""


def _step4_log(env, sample, run_id):
    p = env["output_root"] / sample / f"{sample}_{run_id}" / "logs" / "rctd-split.log"
    return p.read_text() if p.exists() else ""


# ---- Layer 4: named flags -----------------------------------------------

def test_step1_named_flags_threaded(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path),
                    "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "s1flags",
                    "--step1-x-source", "maxpost_counts",
                    "--step1-qc-min-counts-cell", "15",
                    "--step1-gex-only", "false",
                    "--step1-force-rerun")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    log = _step1_log(env, "MH10", "s1flags")
    assert "x_source=maxpost_counts" in log
    assert "qc_min_counts_cell=15" in log
    assert "gex_only=false" in log
    assert "force_rerun=1" in log


def test_step3_named_flags_threaded(env, tmp_path):
    target = tmp_path / "target.json"
    target.write_text("{}")
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path),
                    "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "s3flags",
                    "--step3-donor-borrow-cap", "80",
                    "--step3-cell-min-instance", "25",
                    "--step3-min-umi", "15",
                    "--step3-random-seed", "123",
                    "--step3-celltype-target-list", str(target))
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    log = _step3_log(env, "MH10", "s3flags")
    assert "donor_borrow_cap=80" in log
    assert "cell_min_instance=25" in log
    assert "min_umi=15" in log
    assert "random_seed=123" in log
    assert f"celltype_target_list={target}" in log


def test_step4_named_flags_threaded(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path),
                    "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "s4flags",
                    "--step4-umi-min", "20",
                    "--step4-counts-min", "8",
                    "--step4-cell-min-instance", "30",
                    "--step4-doublet-mode", "full",
                    "--step4-postprocess-min-counts", "75",
                    "--step4-keep-intermediate")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    log = _step4_log(env, "MH10", "s4flags")
    assert "umi_min=20" in log
    assert "counts_min=8" in log
    assert "cell_min_instance=30" in log
    assert "doublet_mode=full" in log
    assert "postprocess_min_counts=75" in log
    assert "keep_intermediate=1" in log


def test_named_flags_unset_do_not_thread(env, tmp_path):
    # Baseline invariant: when no named flags are passed, none of the
    # STEPN_<PARAM> env vars leak into the --export payload — each step
    # falls through to its CLI's own default.yaml value.
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path),
                    "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "no_flags")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    records = _parse_log(env["log"])
    for rec in records:
        for var in (
            "STEP1_X_SOURCE", "STEP1_QC_MIN_COUNTS_CELL", "STEP1_GEX_ONLY",
            "STEP1_FORCE_RERUN",
            "STEP3_DONOR_BORROW_CAP", "STEP3_CELL_MIN_INSTANCE",
            "STEP3_MIN_UMI", "STEP3_RANDOM_SEED",
            "STEP3_CELLTYPE_TARGET_LIST",
            "STEP4_UMI_MIN", "STEP4_COUNTS_MIN", "STEP4_CELL_MIN_INSTANCE",
            "STEP4_DOUBLET_MODE", "STEP4_POSTPROCESS_MIN_COUNTS",
            "STEP4_KEEP_INTERMEDIATE",
            "STEP1_OVERRIDES_B64", "STEP3_OVERRIDES_B64",
            "STEP4_OVERRIDES_B64",
        ):
            assert _export_field(rec, var) is None, (var, rec)


# ---- Layer 3: --override (single-key) -----------------------------------

def test_override_single_key_step3(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path),
                    "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "ov3",
                    "--override", "step3.census.random_seed=999")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    # ref-build received --config pointing at a yaml that carries the
    # override under census.random_seed.
    log = _step3_log(env, "MH10", "ov3")
    assert "config_path=" in log
    assert "config_path=\n" not in log       # non-empty path
    # The mock cats the yaml; check the merged content.
    assert "census:" in log
    assert "random_seed: 999" in log


def test_override_nested_list_yaml_parsed(env, tmp_path):
    # YAML-typed value: a flow-style list — the driver's yaml.safe_load
    # of the RHS must produce a Python list, and yaml.safe_dump must
    # emit a valid block list the step CLI's --config loader can read.
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path),
                    "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "ov_list",
                    "--override", "step4.postprocess.leiden.resolutions=[0.5, 0.7, 0.9]")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    log = _step4_log(env, "MH10", "ov_list")
    assert "postprocess:" in log
    assert "leiden:" in log
    assert "resolutions:" in log
    # Emitted as a YAML block list (default_flow_style=False).
    assert "- 0.5" in log
    assert "- 0.7" in log
    assert "- 0.9" in log


def test_override_multiple_keys_merged(env, tmp_path):
    # Two overrides on the same step get merged into a single yaml.
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path),
                    "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "ov_multi",
                    "--override", "step3.census.random_seed=42",
                    "--override", "step3.census.donor_borrow_cap=200")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    log = _step3_log(env, "MH10", "ov_multi")
    assert "random_seed: 42" in log
    assert "donor_borrow_cap: 200" in log


def test_override_wrong_prefix_rejected(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path),
                    "--celltype-marker-json", _marker_json(tmp_path),
                    "--override", "step2.foo=1")
    assert r.returncode != 0
    assert "step1./step3./step4." in r.stderr


def test_override_bare_key_rejected(env, tmp_path):
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path),
                    "--celltype-marker-json", _marker_json(tmp_path),
                    "--override", "foo=1")
    assert r.returncode != 0
    assert "step1./step3./step4." in r.stderr


# ---- Layer 2: --stepN-config <path> -------------------------------------

def test_step_config_yaml_passed_via_config_flag(env, tmp_path):
    cfg = tmp_path / "step3_user.yaml"
    cfg.write_text(textwrap.dedent("""\
        census:
          cell_min_instance: 33
        rctd_reference_build:
          require_int: false
    """))
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path),
                    "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "s3cfg",
                    "--step3-config", str(cfg))
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    log = _step3_log(env, "MH10", "s3cfg")
    assert "cell_min_instance: 33" in log
    assert "require_int: false" in log


def test_step_config_and_override_merged(env, tmp_path):
    # Overrides layer on top of --stepN-config's contents. When both touch
    # the SAME key, the override wins.
    cfg = tmp_path / "step3_base.yaml"
    cfg.write_text(textwrap.dedent("""\
        census:
          cell_min_instance: 33
          donor_borrow_cap: 111
    """))
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path),
                    "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "s3merge",
                    "--step3-config", str(cfg),
                    "--override", "step3.census.donor_borrow_cap=222",
                    "--override", "step3.census.random_seed=7")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    log = _step3_log(env, "MH10", "s3merge")
    # Preserved from --step3-config
    assert "cell_min_instance: 33" in log
    # Overridden by --override
    assert "donor_borrow_cap: 222" in log
    assert "donor_borrow_cap: 111" not in log
    # New key from --override alone
    assert "random_seed: 7" in log


# ---- Precedence: named flag beats --override + --stepN-config -----------

def test_named_flag_overrides_config_yaml(env, tmp_path):
    # --step3-config sets random_seed=1; --step3-random-seed=9 should win.
    cfg = tmp_path / "seed_cfg.yaml"
    cfg.write_text("census:\n  random_seed: 1\n")
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path),
                    "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "s3prec",
                    "--step3-config", str(cfg),
                    "--step3-random-seed", "9")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    log = _step3_log(env, "MH10", "s3prec")
    # random_seed: 1 is in the yaml but ref-build sees --random-seed 9 too;
    # our mock records whichever `--random-seed` value ref-build gets last,
    # which is the CLI-flag one. This mirrors ref-build's own precedence
    # (CLI flag > user YAML > default.yaml).
    assert "random_seed=9" in log


# ---- Dry-run visibility --------------------------------------------------

def test_dry_run_shows_new_step_env_vars(env, tmp_path):
    # Under --dry-run the driver prints the sbatch command with the full
    # --export payload, so operators can eyeball what each step receives
    # (Tracy's task: "update --dry-run to show what each step would
    # receive"). No fs side effects.
    r = _run_driver(env,
                    "--sample-id", "MH10",
                    "--flex-h5ad", _flex_h5ad(tmp_path),
                    "--celltype-marker-json", _marker_json(tmp_path),
                    "--run-id", "dr_show",
                    "--step1-qc-min-counts-cell", "50",
                    "--step3-min-umi", "12",
                    "--step4-doublet-mode", "full",
                    "--override", "step4.postprocess.qc.min_counts=99",
                    "--dry-run")
    assert r.returncode == 0, f"stderr:\n{r.stderr}\nstdout:\n{r.stdout}"
    combined = r.stdout + r.stderr
    assert "STEP1_QC_MIN_COUNTS_CELL=50" in combined
    assert "STEP3_MIN_UMI=12" in combined
    assert "STEP4_DOUBLET_MODE=full" in combined
    assert "STEP4_OVERRIDES_B64=" in combined
    # No submissions.
    assert _parse_log(env["log"]) == []


def test_help_lists_new_flags(env):
    r = _run_driver(env, "--help")
    assert r.returncode == 0
    for flag in (
        "--step1-x-source", "--step1-qc-min-counts-cell", "--step1-gex-only",
        "--step1-force-rerun",
        "--step3-donor-borrow-cap", "--step3-cell-min-instance",
        "--step3-min-umi", "--step3-random-seed",
        "--step3-celltype-target-list",
        "--step4-umi-min", "--step4-counts-min",
        "--step4-cell-min-instance", "--step4-doublet-mode",
        "--step4-postprocess-min-counts", "--step4-keep-intermediate",
        "--stepN-config", "--override",
    ):
        assert flag in r.stdout, flag
