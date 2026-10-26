"""Independent saved-data audit for the locked compact-dual 5/7 kHz test."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[3]
CASE = ROOT / "results/validation/deepsea3d_5k7k_20261026"
OLD = ROOT / "results/validation/deepsea3d_g1_frequency/frequency.mat"
LOCK = ROOT / "results/validation/deepsea_dual_broadband_decision_20261015/locked_candidate_and_gate.mat"
EVIDENCE = ROOT / "results/validation/deepsea_dual_broadband_decision_20261015/saved_evidence.mat"


def check(condition, label):
    if not condition:
        raise AssertionError(label)


def z(row, stem):
    return complex(float(row[f"{stem}_real"]), float(row[f"{stem}_imag"]))


def all_below(values, limit):
    vals = np.asarray(list(values), dtype=float)
    return bool(np.all(np.isfinite(vals)) and np.all(vals < limit))


def main():
    hashes = (CASE / "output.sha256").read_text().splitlines()
    for line in hashes:
        digest, relative = line.split("  ", 1)
        check(hashlib.sha256((CASE / relative).read_bytes()).hexdigest() == digest, f"SHA: {relative}")
    input_rows = json.loads((CASE / "input_manifest.json").read_text(encoding="utf-8-sig"))
    for row in input_rows:
        check(hashlib.sha256((ROOT / row["path"]).read_bytes()).hexdigest() == row["sha256"], f"input SHA: {row['path']}")
    check(hashlib.sha256(LOCK.read_bytes()).hexdigest() == "f13fed7f6a723bb3239a4726f442f4195491af5f4267703e2a51bcbe975c189d", "prediction lock changed")

    state = loadmat(CASE / "five_seven.mat", simplify_cells=True)["state"]
    old = loadmat(OLD, simplify_cells=True)["state"]
    locked = loadmat(LOCK, simplify_cells=True)
    evidence = loadmat(EVIDENCE, simplify_cells=True)
    cases = [old["cases"][0], state["cases"][0], old["baseline6"], state["cases"][1], old["cases"][1]]
    freqs = [4000, 5000, 6000, 7000, 8000]
    for case, freq in zip(cases, freqs):
        check(int(case["cfg"]["frequency_hz"]) == freq, "frequency or case order")
        check(np.array_equal(case["cfg"]["tx_m"], [0, 150.5, 1500]), "Tx")
        check(np.array_equal(case["cfg"]["rx_m"], [0, 0, 5]), "Rx")
        check(float(case["cfg"]["c0_mps"]) == 1500, "sound speed")
        p = case["cfg"]["profile"]
        check(np.array_equal(p["amplitudes_m"], [.05, .025]), "surface amplitudes")
        check(np.array_equal(p["wavevectors_radpm"], [[.4, 0], [.2, .6]]), "wavevectors")
        check(np.array_equal(p["phases_rad"], [.2, .7]), "phases")
        check(p["plateau_m"] == 8 and p["support_m"] == 12, "support")

    table = list(csv.DictReader((CASE / "five_frequency_native_summary_v2.csv").open(newline="", encoding="utf-8")))
    check(len(table) == 5, "five point table")
    results = []
    eta = evidence["output"]["geometry"]["surface_height_at_flat_specular_m"]
    for j, (case, freq, row) in enumerate(zip(cases, freqs, table)):
        check(int(row["frequency_hz"]) == freq, "table frequency")
        ref = case["runs"][int(case["reference_index"]) - 1]["result"]
        fk = case["fk"][1]
        R = complex(ref["reflect_native"])
        D = np.exp(1j * 2 * np.pi * freq / 1500 * np.linalg.norm(np.array([0, 0, 5.]) - np.array([0, 150.5, 1500.]))) / (4 * np.pi * np.linalg.norm(np.array([0, 0, 5.]) - np.array([0, 150.5, 1500.])))
        L0 = np.linalg.norm(np.array([0, 0, 5.]) - np.array([0, 150.5, -1500.]))
        R0 = -np.exp(1j * 2 * np.pi * freq / 1500 * L0) / (4 * np.pi * L0)
        delta = R - R0
        H = D + R
        for name, value in [("R", R), ("R0", R0), ("D", D), ("Delta", delta), ("H", H), ("FK_R", complex(fk["reflect_native"]))]:
            check(abs(z(row, name) - value) < 2e-16, f"saved complex {freq} {name}")
        check(abs(complex(ref["rough_correction"]) - delta) < 2e-16, "reflection identity")
        check(abs(complex(ref["total_native"]) - H) < 2e-16, "total identity")
        entry = {"frequency_hz": freq, "reference_R_real": R.real, "reference_R_imag": R.imag,
                 "FE_order": int(ref["order"]), "FK_reflection_relative": abs(complex(fk["reflect_native"]) - R) / abs(R)}
        if freq in (5000, 7000):
            n = 0 if freq == 5000 else 1
            check(case["reference_ready"] and case["complete"], "reference readiness")
            guards = case["reference_guards"]
            check(all_below(guards.values(), 1e-4), f"reference guards {freq}")
            check(float(ref["forcing_series_relative"]) < 1e-8, "forcing")
            check(np.all(np.diff(ref["boundary_coefficient_norms"][-5:]) < 0), "coefficient decay")
            check(np.all(np.diff(np.abs(ref["terms"][-5:])) < 0), "pressure term decay")
            check(abs(sum(ref["terms"]) - complex(ref["rough_correction"])) < 1e-16, "FE term representation")
            ev = case["evaluation_controls"]
            check(all_below([ev["term_match"], ev["norm_match"]], 1e-12), "M representation")
            check(all_below([ev["tail32_relative"], ev["tail40_relative"]], 1e-4), "M tails")
            obs = case["observer"]
            check(all_below([obs["term_match_error"], obs["norm_match_error"]], 1e-12), "observer")
            fg = case["fk_guards"]
            check(all_below([fg["grid"], fg["gauss"], fg["aperture_12_16"], fg["aperture_16_24"]], 1e-3), "FullFK numerical")
            metrics = case["metrics"]
            check(all_below([metrics["reflected_complex"], metrics["increment_complex"]], .01), "FullFK physics")
            locked_R = complex(locked["predictions"][n]["reflect_native_real"], locked["predictions"][n]["reflect_native_imag"])
            Lh = np.linalg.norm(np.array([0, 0, 5.]) - np.array([0, 150.5, 2 * eta - 1500.]))
            B = -np.exp(1j * 2 * np.pi * freq / 1500 * Lh) / (4 * np.pi * Lh)
            Q = locked["Q"]
            independently_rebuilt = .5 * (Q[n] + Q[n + 1]) * B
            check(abs(locked_R - independently_rebuilt) < 1e-15, "prelocked formula")
            errR = abs(locked_R - R) / abs(R)
            errDelta = abs(locked_R - R) / abs(delta)
            errH = abs(locked_R - R) / abs(H)
            variants = []
            for index in (case["coarse_index"], case["fine_index"], case["reference_index"]):
                partials = case["runs"][int(index) - 1]["result"]["partials"]
                variants.extend(complex(partials[k]) for k in (-5, -3, -1))
            empirical_FE_gap = max(abs(value - complex(ref["rough_correction"])) for value in variants)
            check(all_below([errR, errDelta, errH], 1e-3), "reconstruction screen")
            for name, val in [("reconstruction_R_relative", errR), ("reconstruction_Delta_relative", errDelta), ("reconstruction_H_relative", errH)]:
                check(abs(float(row[name]) - val) < 2e-15, f"table error {name}")
            entry.update(locked_prediction_real=locked_R.real, locked_prediction_imag=locked_R.imag,
                         locked_prediction_abs_gap=abs(locked_R - R),
                         empirical_FE_nine_variant_abs_gap=empirical_FE_gap,
                         prediction_gap_over_FE_empirical=abs(locked_R - R) / empirical_FE_gap,
                         reconstruction_R_relative=errR, reconstruction_Delta_relative=errDelta,
                         reconstruction_H_relative=errH,
                         boundary_relative=float(guards["boundary"]), halfgrid_relative=float(guards["halfgrid"]),
                         FK_gauss_relative=float(fg["gauss"]))
            # Every retained FE coefficient multiplies an outgoing half-space
            # Fourier mode. Independently check its exact PDE/DtN symbols on
            # the final grid; this structural test is distinct from the saved
            # true-surface boundary and half-grid residuals checked above.
            setting = case["runs"][int(case["reference_index"]) - 1]["settings"]
            L, spacing, _ = map(float, setting)
            grid_n = int(2 * np.ceil(L / spacing))
            kv = (2 * np.pi / (2 * L)) * np.r_[np.arange(grid_n // 2), np.arange(-grid_n // 2, 0)]
            k2 = (2 * np.pi * freq / 1500) ** 2
            horiz2 = kv[:, None] ** 2 + kv[None, :] ** 2
            kz = np.sqrt((k2 - horiz2).astype(complex))
            pde_symbol = float(np.max(np.abs(kz ** 2 + horiz2 - k2)) / k2)
            outgoing_min_imag = float(np.min(kz.imag))
            check(pde_symbol < 1e-12 and outgoing_min_imag >= 0, "outgoing PDE/DtN modal identity")
            entry.update(PDE_modal_symbol_relative=pde_symbol,
                         DtN_outgoing_branch_min_imag=outgoing_min_imag,
                         DtN_symbol="partial_z U_k = i*kz*U_k at z=0; exact for saved outgoing representation")
        results.append(entry)

    jobs = ["diag6000", "f5000", "f7000", "render", "render2", "figaudit", "figaudit2"]
    costs = {}
    for job in jobs:
        directory = CASE / f"run_{job}_01"
        check(int((directory / "exit_status.txt").read_text()) == 0, f"job exit {job}")
        costs[job] = int((directory / "process_wall_s.txt").read_text())
    check(sum(costs.values()) <= 600 and costs["f5000"] < 300 and costs["f7000"] < 300, "budget")
    image = Image.open(CASE / "five_frequency_reflection_and_error_v2.png")
    check(image.width > 1000 and image.height > 450, "figure dimensions")
    check((CASE / "five_frequency_reflection_and_error_v2.fig").stat().st_size > 10000, "editable figure")
    figure_audit = json.loads((CASE / "figure_object_audit_v3.json").read_text())
    check(figure_audit["status"] == "PASS" and figure_audit["axes"] == 2 and figure_audit["line_objects"] == 6,
          "reopened figure objects")
    check(figure_audit["maximum_source_coordinate_error"] < 1e-12, "figure source coordinates")
    report = {"status": "PASS", "output_sha_files": len(hashes), "input_sha_files": len(input_rows),
              "locked_mat_sha256": hashlib.sha256(LOCK.read_bytes()).hexdigest(),
              "frequencies": results, "remote_process_wall_s": costs,
              "remote_process_wall_total_s": sum(costs.values()),
              "figure_pixels": [image.width, image.height],
              "scope": "two new discrete points and two local reconstruction checks; no continuous frequency/CIR claim",
              "G1_target_PM": "OPEN", "WP2_formal": False}
    target = CASE / "independent_audit_v5.json"
    check(not target.exists(), "audit output exists")
    target.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"status": report["status"], "frequencies": results,
                      "remote_process_wall_total_s": report["remote_process_wall_total_s"]}, indent=2))


if __name__ == "__main__":
    main()
