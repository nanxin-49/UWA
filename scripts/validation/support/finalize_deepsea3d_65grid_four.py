"""Close the saved four-point validation with process, SHA, and figure receipts."""
import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/validation/deepsea3d_65grid_four_20261027"
JOBS = ["run_f4750_01", "run_f4750_02", "run_f7250_01", "run_f4125_01",
        "run_f6562.5_01", "run_audit_01", "run_render_01", "run_render_02",
        "run_render_03", "run_figaudit_01"]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    audit = json.loads((OUT / "independent_audit.json").read_text(encoding="utf-8"))
    figures = json.loads((OUT / "figure_object_audit.json").read_text(encoding="utf-8"))
    assert audit["status"] == figures["status"] == "PASS"
    rows = list(csv.DictReader((OUT / "four_point_blind_metrics.csv").open(newline="", encoding="utf-8")))
    assert len(rows) == 4 and all(r["prediction_gates_pass"] == "True" for r in rows)
    costs = {name: int((OUT / name / "process_wall_s.txt").read_text()) for name in JOBS}
    total = 95 + sum(costs.values())
    assert total <= 720
    assert all(costs[name] <= 150 for name in JOBS if name.startswith("run_f"))
    assert all(int((OUT / name / "exit_status.txt").read_text()) == 0 for name in JOBS[1:] if name != "run_render_01")
    assert int((OUT / "run_f4750_01/exit_status.txt").read_text()) == 1
    assert int((OUT / "run_render_01/exit_status.txt").read_text()) == 1
    for line in (OUT / "output.sha256").read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ", 1)
        assert sha(OUT / name) == digest
    for manifest_name in ("figure_output.sha256", "figure_output_v3.sha256"):
        for line in (OUT / manifest_name).read_text(encoding="utf-8").splitlines():
            digest, name = line.split("  ", 1)
            assert sha(OUT / name) == digest
    required = ["prediction_lock.json", "prediction_lock_audit.json", "65grid_four.mat",
                "four_point_blind_metrics.csv", "independent_audit.json", "artifact_assertions.json",
                "figure_object_audit.json", "four_point_reflection_complex_v3.fig",
                "four_point_reflection_complex_v3.png", "prediction_vs_FE_sensitivity_v3.fig",
                "prediction_vs_FE_sensitivity_v3.png"]
    code = ["scripts/validation/support/freeze_deepsea3d_65grid_four.py",
            "scripts/validation/validate_deepsea3d_65grid_four_frequency.m",
            "scripts/validation/run_deepsea3d_65grid_four_frequency.m",
            "scripts/validation/support/audit_deepsea3d_65grid_four_saved.py",
            "scripts/reporting/render_deepsea3d_65grid_four_v3.m",
            "scripts/validation/audit_deepsea3d_65grid_four_figures.m",
            "scripts/validation/support/launch_deepsea3d_65grid_four_remote_v3.sh",
            "scripts/validation/support/launch_deepsea3d_65grid_four_audit.sh",
            "scripts/validation/support/launch_deepsea3d_65grid_four_render_v3.sh",
            "scripts/validation/support/launch_deepsea3d_65grid_four_figaudit.sh"]
    rss = []
    for name in JOBS:
        match = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)",
                          (OUT / name / "stderr.log").read_text(encoding="utf-8", errors="replace"))
        if match:
            rss.append(int(match.group(1)))
    checkpoint = {
        "run_id": "deepsea3d_65grid_four_20261027", "status": "COMPACT_DUAL_65GRID_HOLDOUT_SCREENING_PASS",
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "execution_mode": "Ubuntu SSH through matlab-ubuntu", "matlab_release": "R2023b",
        "completed_steps": ["prediction_lock", "independent_pre_FE_lock_audit", "four_new_FE_references",
                            "four_full_FK_five_configuration_checks", "independent_saved_audit",
                            "two_FIG_PNG_pairs", "FIG_object_audit", "SHA_and_budget_reconciliation"],
        "config": {"frequencies_hz": [4750, 7250, 4125, 6562.5],
                   "groups": ["routine", "routine", "cancellation_stress", "cancellation_stress"],
                   "model": "fixed_eta_five_anchor_Q_piecewise_complex_linear",
                   "prediction_e_R_e_Delta_a_D_strict_lt": 1e-3,
                   "reference_relative_strict_lt": 1e-4,
                   "FullFK_increment_numerical_strict_lt": 1e-3},
        "input_versions_sha256": {name: sha(OUT / name) for name in ["prediction_lock.json", "predictions_pre_FE.csv"]},
        "code_versions_sha256": {name: sha(ROOT / name) for name in code},
        "artifacts_sha256": {name: sha(OUT / name) for name in required},
        "failed_or_unrun_checks": {"physical_reference_failures": 0, "prediction_gate_failures": 0,
                                   "setup_failure_before_FE": 1, "first_render_failure_before_figures": 1,
                                   "continuous_band_and_CIR_not_run": True},
        "resource_usage": {"conservative_version_probe_charge_s": 95,
                           "process_wall_s": costs, "conservative_total_wall_s": total,
                           "budget_s": 720, "max_observed_process_RSS_kB": max(rss)},
        "audit_status_and_path": {"status": "PASS", "saved_audit": str(OUT / "independent_audit.json"),
                                  "figure_audit": str(OUT / "figure_object_audit.json"),
                                  "returned_SHA_count": len((OUT / "output.sha256").read_text().splitlines())},
        "next_allowed_step": "separately authorized time-window and band-extremum validation; no automatic FE expansion or CIR certification",
        "limitations": ["FE sensitivity is empirical, not a rigorous Helmholtz error bound",
                        "95 s MATLAB version probe was observed but its full stdout/stderr was not persisted",
                        "rough group delay cannot be unwrapped reliably from these sparse points"]}
    path = OUT / "checkpoint.json"
    if path.exists():
        snapshot = OUT / "checkpoint_v1.json"
        if not snapshot.exists():
            snapshot.write_bytes(path.read_bytes())
    with path.open("w", encoding="utf-8") as handle:
        json.dump(checkpoint, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    print(json.dumps({"status": checkpoint["status"], "total_wall_s": total,
                      "budget_s": 720, "returned_SHA_count": checkpoint["audit_status_and_path"]["returned_SHA_count"],
                      "max_RSS_kB": max(rss)}))

if __name__ == "__main__":
    main()
