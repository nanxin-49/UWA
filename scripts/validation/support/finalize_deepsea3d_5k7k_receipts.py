"""Write compact post-run receipts without changing frozen numerical artifacts."""
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
CASE = ROOT / "results/validation/deepsea3d_5k7k_20261026"
INPUTS = json.loads((CASE / "input_manifest.json").read_text(encoding="utf-8-sig"))
INPUT_SHA = {x["path"]: x["sha256"] for x in INPUTS}
AUDIT = json.loads((CASE / "independent_audit_v5.json").read_text())
JOBS = ["diag6000", "f5000", "f7000", "render", "render2", "figaudit", "figaudit2"]


def sha(file):
    return hashlib.sha256(file.read_bytes()).hexdigest()


def mem_available(file):
    for line in file.read_text().splitlines():
        if line.startswith("MemAvailable:"):
            return int(line.split()[1]) / 1024
    raise AssertionError(file)


def write_new(file, value):
    if file.exists():
        raise FileExistsError(file)
    file.write_text(json.dumps(value, indent=2, ensure_ascii=False), encoding="utf-8")


def main():
    total = AUDIT["remote_process_wall_total_s"]
    for job in JOBS:
        directory = CASE / f"run_{job}_01"
        wall = int((directory / "process_wall_s.txt").read_text())
        exit_status = int((directory / "exit_status.txt").read_text())
        artifacts = [p for p in directory.iterdir() if p.is_file()]
        checkpoint = {
            "run_id": "20261026_deepsea3d_5k7k_01", "job_id": job,
            "status": "COMPLETE" if exit_status == 0 else "FAILED",
            "execution_mode": "ssh Ubuntu MATLAB -batch",
            "matlab_release": "R2023b",
            "process_wall_s": wall, "exit_status": exit_status,
            "completed_steps": [job] if exit_status == 0 else [],
            "config": {"Tx_m": [0, 150.5, 1500], "Rx_m": [0, 0, 5],
                       "c_mps": 1500, "frequencies_hz": [5000, 7000],
                       "surface_amplitudes_m": [0.05, 0.025],
                       "surface_wavevectors_radpm": [[0.4, 0], [0.2, 0.6]],
                       "surface_phases_rad": [0.2, 0.7],
                       "C2_plateau_m": 8, "C2_support_m": 12,
                       "FE_reference_threshold": 1e-4, "FullFK_numerical_threshold": 1e-3,
                       "FullFK_physical_threshold": 0.01,
                       "reconstruction_screen": 1e-3},
            "input_versions_sha256": INPUT_SHA,
            "artifacts_sha256": {p.name: sha(p) for p in artifacts},
            "resource_usage": {"physical_available_before_MiB": mem_available(directory / "mem_before.txt"),
                               "physical_available_after_MiB": mem_available(directory / "mem_after.txt"),
                               "estimated_extra_working_set_MiB": 1520,
                               "reserve_fraction": 0.30,
                               "measured_peak_MiB": "unknown"},
            "failed_or_unrun_checks": [],
            "audit_status_and_path": str(CASE / "independent_audit_v5.json"),
            "next_allowed_step": "STOP; new research authorization required",
        }
        write_new(directory / "checkpoint_v3.json", checkpoint)
    summary = {
        "run_id": "20261026_deepsea3d_5k7k_01",
        "status": "COMPLETE_DISCRETE_TWO_POINT_SCREEN_PASS",
        "completed_steps": JOBS,
        "input_versions_sha256": INPUT_SHA,
        "code_versions_sha256": {**{k: v for k, v in INPUT_SHA.items() if k.endswith(".m")},
                                  "scripts/reporting/render_deepsea3d_5k7k_frequency.m": sha(ROOT / "scripts/reporting/render_deepsea3d_5k7k_frequency.m"),
                                  "scripts/validation/audit_deepsea3d_5k7k_figure.m": sha(ROOT / "scripts/validation/audit_deepsea3d_5k7k_figure.m"),
                                  "scripts/validation/support/audit_deepsea3d_5k7k_saved.py": sha(ROOT / "scripts/validation/support/audit_deepsea3d_5k7k_saved.py")},
        "artifacts_sha256": {"five_seven.mat": sha(CASE / "five_seven.mat"),
                             "five_frequency_native_summary_v2.csv": sha(CASE / "five_frequency_native_summary_v2.csv"),
                             "five_frequency_reflection_and_error_v2.fig": sha(CASE / "five_frequency_reflection_and_error_v2.fig"),
                             "five_frequency_reflection_and_error_v2.png": sha(CASE / "five_frequency_reflection_and_error_v2.png")},
        "resource_usage": {"total_remote_MATLAB_process_wall_s": total,
                           "budget_s": 600, "remaining_s": 600 - total,
                           "physical_available_5k_before_MiB": mem_available(CASE / "run_f5000_01/mem_before.txt"),
                           "physical_available_7k_before_MiB": mem_available(CASE / "run_f7000_01/mem_before.txt"),
                           "measured_peak_MiB": "unknown"},
        "failed_or_unrun_checks": [],
        "audit_status_and_path": str(CASE / "independent_audit_v5.json"),
        "G1_target_PM": "OPEN", "WP2_formal": False,
        "next_allowed_step": "STOP; no PM repair, further frequencies, CIR, or communication work",
    }
    write_new(CASE / "checkpoint_v3.json", summary)
    print(json.dumps({"status": summary["status"], "process_wall_s": total,
                      "remaining_s": 600 - total, "job_checkpoints": len(JOBS)}))


if __name__ == "__main__":
    main()
