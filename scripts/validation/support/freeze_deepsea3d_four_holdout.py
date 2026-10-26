"""Freeze two compact-dual complex predictors before new FE references exist."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/validation/deepsea3d_four_holdout_20261026"
OLD = ROOT / "results/validation/deepsea3d_5k7k_20261026"
LOCK = ROOT / "results/validation/deepsea_dual_broadband_decision_20261015/locked_candidate_and_gate.mat"
EVIDENCE = ROOT / "results/validation/deepsea_dual_broadband_decision_20261015/saved_evidence.mat"
INPUTS = [
    LOCK, EVIDENCE, OLD / "five_frequency_native_summary_v2.csv",
    OLD / "five_seven.mat",
    ROOT / "results/validation/deepsea3d_g1_frequency/frequency.mat",
    ROOT / "results/validation/deepsea3d_g1_azimuth/azimuth.mat",
    Path(__file__).resolve(),
]
FREQUENCIES = [4500, 5500, 6500, 7500]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    assert not OUT.exists(), "Freeze directory already exists; never overwrite a pre-FE lock"
    assert sha(LOCK) == "f13fed7f6a723bb3239a4726f442f4195491af5f4267703e2a51bcbe975c189d"
    assert not any(str(f) in str(p) for p in INPUTS for f in FREQUENCIES)
    locked = loadmat(LOCK, simplify_cells=True)
    evidence = loadmat(EVIDENCE, simplify_cells=True)
    eta = float(evidence["output"]["geometry"]["surface_height_at_flat_specular_m"])
    assert abs(eta - 0.06251088653876558) < 1e-15
    rows = list(csv.DictReader((OLD / "five_frequency_native_summary_v2.csv").open(newline="", encoding="utf-8")))
    anchors = [int(row["frequency_hz"]) for row in rows]
    assert anchors == [4000, 5000, 6000, 7000, 8000]

    def bstar(f):
        image = np.array([0., 150.5, 2 * eta - 1500.])
        length = np.linalg.norm(np.array([0., 0., 5.]) - image)
        return -np.exp(1j * 2 * np.pi * f / 1500 * length) / (4 * np.pi * length)

    def native(row):
        return complex(float(row["R_real"]), float(row["R_imag"]))

    q5 = [native(row) / bstar(f) for row, f in zip(rows, anchors)]
    q3 = [complex(q) for q in np.atleast_1d(locked["Q"])]
    assert len(q3) == 3
    assert max(abs(q3[i] - q5[2 * i]) for i in range(3)) < 1e-11
    pred = []
    for f in FREQUENCIES:
        a_left = (f - 4000) // 2000
        b_left = (f - 4000) // 1000
        qa = (1 - (f - [4000, 6000][a_left]) / 2000) * q3[a_left] + ((f - [4000, 6000][a_left]) / 2000) * q3[a_left + 1]
        qb = .5 * (q5[b_left] + q5[b_left + 1])
        for model, q in [("A_4_6_8", qa), ("B_4_5_6_7_8", qb)]:
            value = q * bstar(f)
            pred.append({"frequency_hz": f, "predictor": model,
                         "Q_real": float(q.real), "Q_imag": float(q.imag),
                         "Bstar_real": float(bstar(f).real), "Bstar_imag": float(bstar(f).imag),
                         "Rpred_real": float(value.real), "Rpred_imag": float(value.imag)})
    OUT.mkdir(parents=True)
    with (OUT / "predictions_pre_FE.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(pred[0]))
        writer.writeheader()
        writer.writerows(pred)
    manifest = {
        "status": "P1_PREDICTIONS_FROZEN_BEFORE_NEW_FE",
        "model_version": "fixed_eta_Q_piecewise_complex_linear_v1",
        "frequencies_hz": FREQUENCIES,
        "predictor_A_anchor_hz": [4000, 6000, 8000],
        "predictor_B_anchor_hz": anchors,
        "eta_star_m": eta,
        "geometry": {"tx_m": [0, 150.5, 1500], "rx_m": [0, 0, 5], "c_mps": 1500,
                     "source": "exp(ikR)/(4piR)", "time_convention": "exp(-i omega t)",
                     "surface": "frozen compact dual C2 plateau8 support12"},
        "formula": "Q=R/Bstar; piecewise complex linear Q(f); Rpred=Qpred*Bstar(f); Bstar=-exp(i k L_eta)/(4 pi L_eta)",
        "screens": {"reflection_relative_complex_strict_lt": 1e-3,
                    "FE_reference_guards_strict_lt": 1e-4,
                    "FE_replay_strict_lt": 1e-12,
                    "FE_forcing_strict_lt": 1e-8,
                    "FullFK_numerical_strict_lt": 1e-3,
                    "FullFK_reflection_and_increment_strict_lt": 0.01,
                    "empirical_sensitivity": "max absolute complex difference of N,N-2,N-4 partials at coarse/fine/domain vs final FE"},
        "inputs": [{"path": str(path.relative_to(ROOT)).replace("\\", "/"), "sha256": sha(path), "bytes": path.stat().st_size} for path in INPUTS],
        "predictions_csv_sha256": sha(OUT / "predictions_pre_FE.csv"),
        "new_FE_reference_count_at_freeze": 0,
    }
    (OUT / "prediction_lock.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": manifest["status"], "predictions": pred,
                      "csv_sha256": manifest["predictions_csv_sha256"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
