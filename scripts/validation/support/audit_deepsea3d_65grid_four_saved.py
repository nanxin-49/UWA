"""Independent read-only audit of four saved compact-dual FE/FK holdouts."""
import cmath
import csv
import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/validation/deepsea3d_65grid_four_20261027"
ORDER = [4750.0, 7250.0, 4125.0, 6562.5]
GROUP = ["routine", "routine", "cancellation_stress", "cancellation_stress"]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def table(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))

def ascomplex(row, stem):
    return complex(float(row[stem + "_real"]), float(row[stem + "_imag"]))

def phase_error(pred, ref):
    return cmath.phase(pred * ref.conjugate())

def main():
    lock = json.loads((OUT / "prediction_lock.json").read_text(encoding="utf-8"))
    lock_audit = json.loads((OUT / "prediction_lock_audit.json").read_text(encoding="utf-8"))
    assert lock_audit["status"] == "PASS"
    assert sha(OUT / "prediction_lock.json") == lock_audit["prediction_lock_sha256"]
    assert sha(OUT / "predictions_pre_FE.csv") == lock["prediction_csv_sha256"]
    assert all(sha(ROOT / r["path"]) == r["sha256"] for r in lock["inputs"])
    assertions = json.loads((OUT / "artifact_assertions.json").read_text(encoding="utf-8"))
    assert assertions["complete"] and assertions["frequency_hz"] == ORDER
    assert assertions["status"] == ["CONDITIONAL_G1_PASS"] * 4
    assert assertions["six_kHz_new_solves"] == 0 and assertions["frequency_only_change"] and assertions["no_phase_fit"]
    assert assertions["global_G1"] == "OPEN" and not assertions["WP2_allowed"]
    state = loadmat(OUT / "65grid_four.mat", simplify_cells=True)["state"]
    cases = np.atleast_1d(state["cases"])
    assert len(cases) == 4
    predictions = table(OUT / "predictions_pre_FE.csv")
    assert [float(r["frequency_hz"]) for r in predictions] == ORDER
    run_names = ["run_f4750_01", "run_f4750_02", "run_f7250_01", "run_f4125_01", "run_f6562.5_01", "run_audit_01"]
    starts = []
    process_wall = 95  # Conservative charge for the actual R2023b version probe.
    process_rows = []
    for name in run_names:
        folder = OUT / name
        start = (folder / "start_utc.txt").read_text().strip()
        end = (folder / "end_utc.txt").read_text().strip()
        starts.append(start)
        exit_status = int((folder / "exit_status.txt").read_text())
        wall = int((folder / "process_wall_s.txt").read_text())
        process_wall += wall
        rss_match = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", (folder / "stderr.log").read_text(errors="replace"))
        rss = int(rss_match.group(1)) if rss_match else None
        mem = (folder / "mem_before.txt").read_text()
        total = int(re.search(r"^MemTotal:\s*(\d+)", mem, re.M).group(1))
        available = int(re.search(r"^MemAvailable:\s*(\d+)", mem, re.M).group(1))
        assert available - 3145728 >= .30 * total
        disk_lines = (folder / "disk_before.txt").read_text(encoding="utf-8", errors="replace").splitlines()
        assert int(disk_lines[-1].split()[3]) >= 2097152
        process_rows.append({"job": name, "start_utc": start, "end_utc": end, "exit_status": exit_status,
                             "wall_s": wall, "max_rss_kB": rss, "mem_available_before_kB": available})
    assert starts == sorted(starts)
    assert datetime.fromisoformat(lock["generated_utc"]) < datetime.fromisoformat(starts[0].replace("Z", "+00:00"))
    assert datetime.fromisoformat(lock_audit["audited_utc"]) < datetime.fromisoformat(starts[0].replace("Z", "+00:00"))
    assert process_rows[0]["exit_status"] == 1 and "Required project directory is missing" in (OUT / run_names[0] / "stderr.log").read_text(errors="replace")
    assert all(x["exit_status"] == 0 for x in process_rows[1:])
    assert all(x["wall_s"] <= 150 for x in process_rows if x["job"].startswith("run_f"))
    assert process_wall <= 720
    all_rows = []
    for i, (f, group, c, pred) in enumerate(zip(ORDER, GROUP, cases, predictions)):
        assert float(c["frequency_hz"]) == f and bool(c["complete"]) and bool(c["reference_ready"])
        assert float(c["cfg"]["frequency_hz"]) == f
        folder = OUT / ("f%g" % f)
        guards = table(folder / "reference_guards.csv")[0]
        fk_guard = table(folder / "fk_guards.csv")[0]
        physical = table(folder / "physical_metrics.csv")[0]
        flat = table(folder / "flat_controls.csv")[0]
        eval_checks = table(folder / "evaluation_checks.csv")
        assert all(math.isfinite(float(v)) and float(v) < 1e-4 for v in guards.values())
        assert all(float(fk_guard[k]) < 1e-3 for k in ["grid", "gauss", "aperture_12_16", "aperture_16_24"])
        assert float(physical["reflected_complex"]) < .01 and float(physical["increment_complex"]) < .01
        assert all(float(v) < 1e-10 for k, v in flat.items() if "error" in k)
        # N12 preflight can fail forcing; the selected final reference must pass.
        assert float(eval_checks[int(c["reference_index"]) - 1]["forcing_error"]) < 1e-8
        observer = c["observer"]
        assert float(observer["term_match_error"]) < 1e-12 and float(observer["norm_match_error"]) < 1e-12
        assert len(np.atleast_1d(c["fk"])) == 5
        ref = c["runs"][int(c["reference_index"]) - 1]["result"]
        r = complex(ref["reflect_native"])
        delta = complex(ref["rough_correction"])
        h = complex(ref["total_native"])
        d = ascomplex(pred, "D")
        r0 = ascomplex(pred, "R0")
        assert abs(r - (r0 + delta)) < 1e-15 and abs(h - (d + r)) < 1e-15
        assert abs(r - ascomplex(table(folder / "reference_runs.csv")[-1], "reflect")) < 1e-15
        assert abs(delta - ascomplex(table(folder / "reference_runs.csv")[-1], "delta")) < 1e-15
        assert abs(h - ascomplex(table(folder / "reference_runs.csv")[-1], "total")) < 1e-15
        rpred = ascomplex(pred, "Rpred")
        dpred = ascomplex(pred, "Deltapred")
        hpred = ascomplex(pred, "Hpred")
        assert abs(hpred - d - rpred) < 1e-15 and abs(dpred - rpred + r0) < 1e-15
        gap = abs(rpred - r)
        variants = []
        for key in ("coarse_index", "fine_index", "reference_index"):
            part = np.atleast_1d(c["runs"][int(c[key]) - 1]["result"]["partials"])
            variants.extend(abs(complex(v) - delta) for v in part[-3:])
        sensitivity = max(variants)
        assert sensitivity > 0
        ld = float(c["flat_control"]["Ld_m"])
        l0 = float(c["flat_control"]["Lr_m"])
        tau = ld / 1500
        direct_dsp = d.conjugate() * cmath.exp(2j * math.pi * f * tau)
        flat_rel = r0.conjugate() * cmath.exp(2j * math.pi * f * tau)
        expected_flat_rel = -cmath.exp(-2j * math.pi * f * (l0 - ld) / 1500) / (4 * math.pi * l0)
        assert abs(direct_dsp - 1 / (4 * math.pi * ld)) < 1e-14
        assert abs(flat_rel - expected_flat_rel) < 1e-14
        er = gap / abs(r)
        ed = abs(dpred - delta) / abs(delta)
        eh = abs(hpred - h) / abs(h)
        ad = abs(hpred - h) / abs(d)
        assert er < 1e-3 and ed < 1e-3 and ad < 1e-3
        row = {"frequency_hz": f, "group": group, "FE_order": int(ref["order"]),
               "FE_L_m": float(c["runs"][int(c["reference_index"]) - 1]["settings"][0]),
               "FE_h_m": float(ref["spacing_m"]),
               "D_real": d.real, "D_imag": d.imag, "R0_real": r0.real, "R0_imag": r0.imag,
               "R_FE_real": r.real, "R_FE_imag": r.imag, "Delta_FE_real": delta.real, "Delta_FE_imag": delta.imag,
               "H_FE_real": h.real, "H_FE_imag": h.imag,
               "Rpred_real": rpred.real, "Rpred_imag": rpred.imag,
               "Deltapred_real": dpred.real, "Deltapred_imag": dpred.imag,
               "Hpred_real": hpred.real, "Hpred_imag": hpred.imag,
               "absolute_complex_gap_m_inv": gap, "e_R": er, "e_Delta": ed, "e_H": eh, "a_D": ad,
               "R_magnitude_relative": abs(abs(rpred) - abs(r)) / abs(r),
               "Delta_magnitude_relative": abs(abs(dpred) - abs(delta)) / abs(delta),
               "H_magnitude_relative": abs(abs(hpred) - abs(h)) / abs(h),
               "R_phase_error_rad": phase_error(rpred, r),
               "Delta_phase_error_rad": phase_error(dpred, delta),
               "H_phase_error_rad": phase_error(hpred, h),
               "FE_empirical_sensitivity_m_inv": sensitivity,
               "prediction_gap_over_FE_empirical": gap / sensitivity,
               "cancellation_C": (abs(d) + abs(r)) / abs(h),
               "FE_boundary_relative": float(guards["boundary"]),
               "FE_halfgrid_relative": float(guards["halfgrid"]),
               "FK_gauss_increment_relative": float(fk_guard["gauss"]),
               "FK_FE_reflection_complex": float(physical["reflected_complex"]),
               "FK_FE_increment_complex": float(physical["increment_complex"]),
               "prediction_gates_pass": True}
        all_rows.append(row)
    with (OUT / "four_point_blind_metrics.csv").open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)
    audit = {"status": "PASS", "audited_utc": datetime.now(timezone.utc).isoformat(),
             "prediction_lock_sha256": sha(OUT / "prediction_lock.json"),
             "first_new_MATLAB_start_utc": starts[0],
             "first_new_FE_start_utc": starts[1], "prediction_lock_preceded_first_FE": True,
             "frequencies_hz": ORDER, "reference_qualified_count": 4,
             "prediction_gate_pass_count": 4,
             "max_e_R": max(x["e_R"] for x in all_rows),
             "max_e_Delta": max(x["e_Delta"] for x in all_rows),
             "max_a_D": max(x["a_D"] for x in all_rows),
             "max_e_H": max(x["e_H"] for x in all_rows),
             "conservative_process_wall_s": process_wall,
             "processes": process_rows, "returned_file_sha_count": len((OUT / "output.sha256").read_text().splitlines()),
             "frequency_phase": "native→DSP analytic direct/flat closure; rough group-delay branch unresolved on sparse points",
             "new_physical_solver_frequency_count": 4,
             "matlab_release": "R2023b", "global_G1_target_PM": "OPEN", "WP2_formal": False}
    with (OUT / "independent_audit.json").open("x", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    print(json.dumps({k: audit[k] for k in ["status", "reference_qualified_count", "prediction_gate_pass_count",
                                              "max_e_R", "max_e_Delta", "max_a_D", "max_e_H",
                                              "conservative_process_wall_s", "returned_file_sha_count"]}))

if __name__ == "__main__":
    main()
