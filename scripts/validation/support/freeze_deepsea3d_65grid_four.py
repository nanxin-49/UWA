"""Freeze four fixed compact-dual predictions without reading new FE results."""
import cmath
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/validation/deepsea3d_65grid_four_20261027"
OLD = ROOT / "results/validation/deepsea3d_5k7k_20261026"
LOCKED = ROOT / "results/validation/deepsea_dual_broadband_decision_20261015/locked_candidate_and_gate.mat"
ANCHOR_CSV = OLD / "five_frequency_native_summary_v2.csv"
SCRIPT = Path(__file__).resolve()
FREQUENCIES = [(4750.0, "routine"), (7250.0, "routine"), (4125.0, "cancellation_stress"),
               (6562.5, "cancellation_stress")]
ETA = 0.06251088653876558
TX = (0.0, 150.5, 1500.0)
RX = (0.0, 0.0, 5.0)
C = 1500.0


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def distance(image_depth):
    return math.sqrt((TX[0] - RX[0]) ** 2 + (TX[1] - RX[1]) ** 2 + (image_depth - RX[2]) ** 2)


def field(f, length, sign):
    return sign * cmath.exp(2j * math.pi * f * length / C) / (4 * math.pi * length)


def main():
    assert sha(LOCKED) == "f13fed7f6a723bb3239a4726f442f4195491af5f4267703e2a51bcbe975c189d"
    assert sha(ANCHOR_CSV) == "36670224216b524dbb4f8abb13fa6ddda8bf331c65dd7ac43f0b69a805f49ac8"
    with ANCHOR_CSV.open(newline="", encoding="utf-8") as handle:
        anchors = list(csv.DictReader(handle))
    af = [float(r["frequency_hz"]) for r in anchors]
    assert af == [4000, 5000, 6000, 7000, 8000]
    ld, l0, le = distance(TX[2]), distance(-TX[2]), distance(2 * ETA - TX[2])
    q = [complex(float(r["R_real"]), float(r["R_imag"])) / field(f, le, -1)
         for f, r in zip(af, anchors)]
    rows = []
    for f, group in FREQUENCIES:
        i = next(i for i in range(4) if af[i] < f < af[i + 1])
        w = (f - af[i]) / (af[i + 1] - af[i])
        qp = (1 - w) * q[i] + w * q[i + 1]
        values = {"D": field(f, ld, 1), "R0": field(f, l0, -1),
                  "Bstar": field(f, le, -1), "Qpred": qp}
        values["Rpred"] = values["Bstar"] * qp
        values["Deltapred"] = values["Rpred"] - values["R0"]
        values["Hpred"] = values["D"] + values["Rpred"]
        row = {"frequency_hz": f, "group": group, "left_anchor_hz": af[i],
               "right_anchor_hz": af[i + 1], "weight": w}
        for name, value in values.items():
            row[name + "_real"] = value.real
            row[name + "_imag"] = value.imag
        rows.append(row)
    OUT.mkdir(parents=True, exist_ok=False)
    csv_path = OUT / "predictions_pre_FE.csv"
    with csv_path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    inputs = [LOCKED, ANCHOR_CSV, OLD / "five_seven.mat",
              ROOT / "results/validation/deepsea3d_g1_frequency/frequency.mat",
              ROOT / "results/validation/deepsea3d_g1_azimuth/azimuth.mat", SCRIPT]
    lock = {
        "status": "PREDICTIONS_FROZEN_BEFORE_NEW_FE",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "result_dir": str(OUT),
        "model": "five_anchor_fixed_eta_Bstar_piecewise_complex_linear_v1",
        "anchor_hz": af, "frequencies_hz_in_execution_order": [f for f, _ in FREQUENCIES],
        "eta_star_m": ETA, "tx_m": TX, "rx_m": RX, "c_mps": C,
        "time_convention": "native exp(-i omega t)",
        "source": "unit 3D point Green exp(ikR)/(4piR)",
        "surface": "fixed compact-dual C2 plateau8/support12",
        "predictor": "Q=R/Bstar; complex linear Q between nearest 1kHz anchors; Rpred=Bstar*Q",
        "gates": {"reference_relative_strict_lt": 1e-4,
                  "FullFK_increment_numerical_strict_lt": 1e-3,
                  "FullFK_reflection_and_increment_strict_lt": 0.01,
                  "e_R_strict_lt": 1e-3, "e_Delta_strict_lt": 1e-3,
                  "a_D_strict_lt": 1e-3,
                  "FE_replay_strict_lt": 1e-12, "FE_forcing_strict_lt": 1e-8},
        "limits": {"remote_matlab_wall_total_s": 720, "physical_frequency_wall_s": 150,
                   "ram_reserve_gib": 3, "min_physical_ram_fraction": 0.30,
                   "disk_reserve_gib": 2},
        "prediction_rows": rows,
        "prediction_csv_sha256": sha(csv_path),
        "inputs": [{"path": str(p.relative_to(ROOT)).replace("\\", "/"),
                    "sha256": sha(p), "bytes": p.stat().st_size} for p in inputs],
        "new_FE_reference_count_at_freeze": 0,
    }
    lock_path = OUT / "prediction_lock.json"
    with lock_path.open("x", encoding="utf-8") as handle:
        json.dump(lock, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    print(json.dumps({"status": lock["status"], "generated_utc": lock["generated_utc"],
                      "result_dir": str(OUT), "prediction_csv_sha256": lock["prediction_csv_sha256"],
                      "lock_sha256": sha(lock_path)}))


if __name__ == "__main__":
    main()
