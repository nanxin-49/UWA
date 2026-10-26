"""Independent complex-field audit and plots for the four prelocked FE holdouts."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/validation/deepsea3d_four_holdout_20261026"
OLD = ROOT / "results/validation/deepsea3d_5k7k_20261026/five_frequency_native_summary_v2.csv"
FREQ = [4500, 5500, 6500, 7500]


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def z(row, name):
    return complex(float(row[name + "_real"]), float(row[name + "_imag"]))


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def native_background(f):
    tx = np.array([0., 150.5, 1500.])
    rx = np.array([0., 0., 5.])
    direct_length = np.linalg.norm(rx - tx)
    image = tx.copy()
    image[2] = -image[2]
    reflect_length = np.linalg.norm(rx - image)
    k = 2 * np.pi * f / 1500
    return (np.exp(1j * k * direct_length) / (4 * np.pi * direct_length),
            -np.exp(1j * k * reflect_length) / (4 * np.pi * reflect_length))


def main():
    lock = json.loads((OUT / "prediction_lock.json").read_text(encoding="utf-8"))
    check(lock["new_FE_reference_count_at_freeze"] == 0, "pre-FE freeze")
    check(digest(OUT / "predictions_pre_FE.csv") == lock["predictions_csv_sha256"], "prediction SHA")
    for row in lock["inputs"]:
        check(digest(ROOT / row["path"]) == row["sha256"], "anchor input SHA: " + row["path"])
    returned = (OUT / "remote_output.sha256").read_text().splitlines()
    check(len(returned) == 106, "returned file count")
    for line in returned:
        sha, rel = line.split("  ", 1)
        check(digest(OUT / rel) == sha, "return SHA: " + rel)
    assertions = json.loads((OUT / "artifact_assertions.json").read_text())
    check(assertions["complete"] and assertions["frequency_hz"] == FREQ, "MATLAB audit cases")
    check(assertions["status"] == ["CONDITIONAL_G1_PASS"] * 4, "MATLAB physical qualification")
    check(assertions["global_G1"] == "OPEN" and not assertions["WP2_allowed"], "scope")
    state = loadmat(OUT / "four_holdout.mat", simplify_cells=True)["state"]
    cases = list(np.atleast_1d(state["cases"]))
    check(len(cases) == 4 and "R2023b" in state["matlab_version"], "MATLAB state or release")
    predictions = {(int(r["frequency_hz"]), r["predictor"]): z(r, "Rpred")
                   for r in csv.DictReader((OUT / "predictions_pre_FE.csv").open(newline="", encoding="utf-8"))}
    check(len(predictions) == 8, "eight frozen predictions")
    anchors = list(csv.DictReader(OLD.open(newline="", encoding="utf-8")))
    check([int(r["frequency_hz"]) for r in anchors] == [4000, 5000, 6000, 7000, 8000], "anchors")

    total_wall = 0
    max_rss_kb = 0
    qualifier = []
    blind = []
    truth = {}
    for f, case in zip(FREQ, cases):
        check(int(case["cfg"]["frequency_hz"]) == f, "case frequency")
        cfg = case["cfg"]
        p = cfg["profile"]
        check(np.array_equal(cfg["tx_m"], [0, 150.5, 1500]) and np.array_equal(cfg["rx_m"], [0, 0, 5]), "geometry")
        check(float(cfg["c0_mps"]) == 1500 and np.array_equal(p["amplitudes_m"], [.05, .025]), "source or amplitudes")
        check(np.array_equal(p["wavevectors_radpm"], [[.4, 0], [.2, .6]]) and np.array_equal(p["phases_rad"], [.2, .7]), "surface spectrum")
        check(p["plateau_m"] == 8 and p["support_m"] == 12 and case["flat_control"]["passed"], "support or plane")
        check(case["reference_ready"] and case["complete"], "reference completion")
        runs = list(np.atleast_1d(case["runs"]))
        ref = runs[int(case["reference_index"]) - 1]["result"]
        settings = runs[int(case["reference_index"]) - 1]["settings"]
        check(settings[0] in (32, 48) and settings[1] == .0625 and int(ref["order"]) in (16, 20), "reference grid")
        guards = case["reference_guards"]
        check(all(np.isfinite(float(v)) and float(v) < 1e-4 for v in guards.values()), "FE reference guards")
        check(ref["boundary_residual_relative"] < 1e-4 and ref["offgrid_boundary_residual_relative"] < 1e-4, "true surface boundary")
        check(ref["forcing_series_relative"] < 1e-8, "forcing")
        L, spacing = map(float, settings[:2])
        grid_n = int(2 * np.ceil(L / spacing))
        kv = (2 * np.pi / (2 * L)) * np.r_[np.arange(grid_n // 2), np.arange(-grid_n // 2, 0)]
        k2 = (2 * np.pi * f / 1500) ** 2
        horizontal2 = kv[:, None] ** 2 + kv[None, :] ** 2
        kz = np.sqrt((k2 - horizontal2).astype(complex))
        modal_symbol_residual = float(np.max(np.abs(kz ** 2 + horizontal2 - k2)) / k2)
        check(modal_symbol_residual < 1e-12 and np.min(kz.imag) >= 0, "outgoing discrete PDE/DtN symbol")
        check(np.all(np.diff(ref["boundary_coefficient_norms"][-5:]) < 0), "coefficient tail")
        check(np.all(np.diff(np.abs(ref["terms"][-5:])) < 0), "pressure tail")
        check(abs(sum(ref["terms"]) - ref["rough_correction"]) < 1e-16, "FE complex series")
        ev = case["evaluation_controls"]
        check(max(float(ev[k]) for k in ("term_match", "norm_match")) < 1e-12, "M representation")
        check(max(float(ev[k]) for k in ("tail32_relative", "tail40_relative")) < 1e-4, "M tails")
        obs = case["observer"]
        check(max(float(obs[k]) for k in ("term_match_error", "norm_match_error")) < 1e-12, "observer")
        fk = list(np.atleast_1d(case["fk"]))
        check(len(fk) == 5, "full FK numerical set")
        def rel(a, b):
            return abs(a - b) / abs(b)
        fk_guard = [rel(fk[0]["rough_correction"], fk[1]["rough_correction"]),
                    rel(fk[2]["rough_correction"], fk[1]["rough_correction"]),
                    rel(fk[3]["rough_correction"], fk[1]["rough_correction"]),
                    rel(fk[4]["rough_correction"], fk[1]["rough_correction"])]
        check(max(fk_guard) < 1e-3, "FullFK numerical gate")
        check(max(abs(float(fk_guard[i] - case["fk_guards"][k])) for i, k in enumerate(
            ("grid", "gauss", "aperture_12_16", "aperture_16_24"))) < 1e-12, "FullFK guard receipt")
        R = complex(ref["reflect_native"])
        D, R0 = native_background(f)
        delta = R - R0
        H = D + R
        check(abs(complex(ref["rough_correction"]) - delta) < 2e-16, "rough increment identity")
        check(abs(complex(ref["total_native"]) - H) < 2e-16, "total field identity")
        fk_R = complex(fk[1]["reflect_native"])
        check(rel(fk_R, R) < .01 and rel(complex(fk[1]["rough_correction"]), delta) < .01, "FullFK physics")
        empirical = max(abs(complex(partial) - complex(ref["rough_correction"]))
                        for index in (case["coarse_index"], case["fine_index"], case["reference_index"])
                        for partial in np.asarray(runs[int(index) - 1]["result"]["partials"])[[-5, -3, -1]])
        check(empirical > 0, "empirical sensitivity nonzero")
        saved_rows = list(csv.DictReader((OUT / f"f{f}/reference_runs.csv").open(newline="", encoding="utf-8")))
        selected = saved_rows[int(case["reference_index"]) - 1]
        check(abs(z(selected, "reflect") - R) < 1e-16 and abs(z(selected, "delta") - delta) < 1e-16, "reference CSV/MAT")
        process = OUT / f"run_f{f}_01"
        check((process / "exit_status.txt").read_text().strip() == "0", "process exit")
        wall = int((process / "process_wall_s.txt").read_text())
        check(wall <= 85 and wall < 120, "frequency wall")
        total_wall += wall
        stderr = (process / "stderr.log").read_text(errors="replace")
        rss_line = next(line for line in stderr.splitlines() if "Maximum resident set size" in line)
        rss = int(rss_line.rsplit(":", 1)[1].strip())
        max_rss_kb = max(max_rss_kb, rss)
        mem = (process / "mem_before.txt").read_text()
        mem_values = {line.split(":", 1)[0]: int(line.split(":", 1)[1].split()[0])
                      for line in mem.splitlines() if line.startswith(("MemTotal:", "MemAvailable:"))}
        check(mem_values["MemAvailable"] - 1556480 >= .30 * mem_values["MemTotal"], "30 percent physical RAM reserve")
        qualifier.append({"frequency_hz": f, "origin": "new_independent_FE",
                          "reference_qualification": "CONDITIONAL_G1_PASS", "FE_order": int(ref["order"]),
                          "FE_L_m": float(settings[0]), "FE_h_m": float(settings[1]),
                          "R_real": R.real, "R_imag": R.imag, "R0_real": R0.real, "R0_imag": R0.imag,
                          "D_real": D.real, "D_imag": D.imag, "Delta_real": delta.real, "Delta_imag": delta.imag,
                          "H_real": H.real, "H_imag": H.imag,
                          "boundary_relative": float(ref["boundary_residual_relative"]),
                          "halfgrid_relative": float(ref["offgrid_boundary_residual_relative"]),
                          "order_Nminus4_relative": float(guards["order_Nminus4"]),
                          "grid_relative": float(guards["grid"]), "domain_relative": float(guards["domain"]),
                          "M32_tail_relative": float(ev["tail32_relative"]),
                          "outgoing_modal_symbol_relative": modal_symbol_residual,
                          "FK_max_numerical_relative": float(max(fk_guard)),
                          "FK_reflection_relative": float(rel(fk_R, R)),
                          "FE_empirical_abs_m_inv": float(empirical)})
        truth[f] = (R, R0, D, delta, H, empirical)
        for name in ("A_4_6_8", "B_4_5_6_7_8"):
            pred = predictions[(f, name)]
            pred_delta = pred - R0
            pred_H = D + pred
            gap = abs(pred - R)
            phase = float(np.angle(pred * np.conj(R)))
            ratio = gap / empirical
            blind.append({"frequency_hz": f, "predictor": name,
                          "R_pred_real": pred.real, "R_pred_imag": pred.imag,
                          "R_FE_real": R.real, "R_FE_imag": R.imag,
                          "absolute_complex_gap_m_inv": gap,
                          "R_relative_complex": gap / abs(R),
                          "Delta_relative_complex": gap / abs(delta),
                          "H_relative_complex": gap / abs(H),
                          "H_error_over_direct_free_field": gap / abs(D),
                          "R_relative_magnitude": abs(abs(pred) - abs(R)) / abs(R),
                          "R_phase_error_rad": phase,
                          "Delta_relative_magnitude": abs(abs(pred_delta) - abs(delta)) / abs(delta),
                          "Delta_phase_error_rad": float(np.angle(pred_delta * np.conj(delta))),
                          "H_relative_magnitude": abs(abs(pred_H) - abs(H)) / abs(H),
                          "H_phase_error_rad": float(np.angle(pred_H * np.conj(H))),
                          "FE_empirical_abs_m_inv": empirical,
                          "gap_over_FE_empirical": ratio,
                          "R_screen_pass": gap / abs(R) < 1e-3,
                          "difference_resolvable_5x": ratio > 5,
                          "near_cancellation_C": (abs(D) + abs(R)) / abs(H)})
    old_rows = []
    for row in anchors:
        f = int(row["frequency_hz"])
        old_rows.append({"frequency_hz": f, "origin": "existing_FE_anchor",
                         "reference_qualification": "CONDITIONAL_G1_PASS",
                         "FE_order": int(row["FE_order"]), "FE_L_m": "", "FE_h_m": "",
                         **{k: row[k] for k in ("R_real", "R_imag", "R0_real", "R0_imag", "D_real", "D_imag",
                                                    "Delta_real", "Delta_imag", "H_real", "H_imag")},
                         **{k: "" for k in ("boundary_relative", "halfgrid_relative", "order_Nminus4_relative",
                                              "grid_relative", "domain_relative", "M32_tail_relative",
                                              "outgoing_modal_symbol_relative", "FK_max_numerical_relative",
                                              "FK_reflection_relative", "FE_empirical_abs_m_inv")}})
    audit_wall = int((OUT / "run_audit_01/process_wall_s.txt").read_text())
    check((OUT / "run_audit_01/exit_status.txt").read_text().strip() == "0", "MATLAB audit exit")
    check(total_wall + audit_wall <= 400, "goal process budget")
    nine = sorted(old_rows + qualifier, key=lambda row: row["frequency_hz"])
    write_csv(OUT / "nine_frequency_qualification.csv", nine)
    write_csv(OUT / "holdout_blind_metrics.csv", blind)

    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": .25})
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(12, 4.8), layout="constrained")
    all_f = np.array([row["frequency_hz"] for row in nine])
    all_R = np.array([z(row, "R") for row in nine])
    ax0.scatter(all_R[[0, 2, 4, 6, 8]].real * 1e5, all_R[[0, 2, 4, 6, 8]].imag * 1e5,
                marker="o", color="#355C7D", label="Anchors (4–8 kHz)")
    ax0.scatter(all_R[[1, 3, 5, 7]].real * 1e5, all_R[[1, 3, 5, 7]].imag * 1e5,
                marker="s", color="#C06C35", label="Independent FE holdouts")
    for f in FREQ:
        R = truth[f][0]
        ax0.annotate(f"{f/1000:g}", (R.real * 1e5, R.imag * 1e5), xytext=(4, 4), textcoords="offset points")
    ax0.set(xlabel="Re(R) × 10⁻⁵ m⁻¹", ylabel="Im(R) × 10⁻⁵ m⁻¹", title="Native reflected complex response")
    ax0.legend(loc="best", fontsize=8)
    for name, marker, color in (("A_4_6_8", "o", "#355C7D"), ("B_4_5_6_7_8", "s", "#C06C35")):
        rows = [r for r in blind if r["predictor"] == name]
        ax1.scatter(np.array(FREQ) / 1000, [r["R_relative_complex"] for r in rows], marker=marker,
                    color=color, label=name)
    ax1.axhline(1e-3, color="black", linestyle="--", linewidth=1, label="Frozen 10⁻³ screen")
    ax1.set(xlabel="Frequency (kHz)", ylabel="|Rpred−RFE| / |RFE|", title="Four prelocked reconstruction tests", yscale="log")
    ax1.legend(fontsize=8)
    fig.savefig(OUT / "reflection_complex_and_reconstruction.png", dpi=180)
    fig.savefig(OUT / "reflection_complex_and_reconstruction.svg")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 4.8), layout="constrained")
    x = np.arange(4)
    width = .25
    for offset, name, color in ((-width, "A_4_6_8", "#355C7D"), (0, "B_4_5_6_7_8", "#C06C35")):
        vals = [r["absolute_complex_gap_m_inv"] for r in blind if r["predictor"] == name]
        ax.bar(x + offset, vals, width, label=name, color=color)
    ax.bar(x + width, [truth[f][-1] for f in FREQ], width, label="FE empirical sensitivity", color="#8F9678")
    ax.set(xticks=x, xticklabels=[f"{f/1000:g}" for f in FREQ], xlabel="Frequency (kHz)",
           ylabel="Absolute complex difference (m⁻¹)", title="Prediction gaps and FE numerical sensitivity")
    ax.legend(fontsize=8)
    fig.savefig(OUT / "prediction_vs_reference_sensitivity.png", dpi=180)
    fig.savefig(OUT / "prediction_vs_reference_sensitivity.svg")
    plt.close(fig)

    summary = {"status": "PASS", "reference_qualified_count": 4,
               "prediction_screen_pass_count": sum(r["R_screen_pass"] for r in blind),
               "total_new_MATLAB_process_wall_s": total_wall + audit_wall,
               "frequency_process_wall_s": {str(f): int((OUT / f"run_f{f}_01/process_wall_s.txt").read_text()) for f in FREQ},
               "audit_process_wall_s": audit_wall, "max_process_RSS_kB": max_rss_kb,
               "returned_SHA_count": len(returned), "prediction_SHA": lock["predictions_csv_sha256"],
               "maximum_R_relative_complex": {name: max(r["R_relative_complex"] for r in blind if r["predictor"] == name)
                                              for name in ("A_4_6_8", "B_4_5_6_7_8")},
               "maximum_gap_over_FE_empirical": {name: max(r["gap_over_FE_empirical"] for r in blind if r["predictor"] == name)
                                                 for name in ("A_4_6_8", "B_4_5_6_7_8")}}
    (OUT / "independent_audit.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
