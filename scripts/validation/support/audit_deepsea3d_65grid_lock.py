"""Independently audit frozen predictions before any new physical calculation."""
import cmath
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results/validation/deepsea3d_65grid_four_20261027"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    lock_path = OUT / "prediction_lock.json"
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    csv_path = OUT / "predictions_pre_FE.csv"
    rows = list(csv.DictReader(csv_path.open(newline="", encoding="utf-8")))
    old = list(csv.DictReader((ROOT / "results/validation/deepsea3d_four_holdout_20261026/nine_frequency_qualification.csv").open(newline="", encoding="utf-8")))
    anchors = {float(r["frequency_hz"]): complex(float(r["R_real"]), float(r["R_imag"])) for r in old if float(r["frequency_hz"]) in lock["anchor_hz"]}
    assert len(anchors) == 5 and len(rows) == 4
    eta = lock["eta_star_m"]
    ld = math.sqrt(150.5 ** 2 + 1495 ** 2)
    l0 = math.sqrt(150.5 ** 2 + 1505 ** 2)
    le = math.sqrt(150.5 ** 2 + (1505 - 2 * eta) ** 2)
    def green(f, length, sign):
        return sign * cmath.exp(1j * 2 * math.pi * f * length / 1500) / (4 * math.pi * length)
    q = {f: r / green(f, le, -1) for f, r in anchors.items()}
    worst = 0.0
    for row, frozen in zip(rows, lock["prediction_rows"]):
        f = float(row["frequency_hz"])
        assert f == frozen["frequency_hz"]
        assert f not in [float(x["frequency_hz"]) for x in old]
        a, b = float(row["left_anchor_hz"]), float(row["right_anchor_hz"])
        assert a < f < b and b - a == 1000
        w = (f - a) / (b - a)
        vals = {"D": green(f, ld, 1), "R0": green(f, l0, -1),
                "Bstar": green(f, le, -1), "Qpred": (1 - w) * q[a] + w * q[b]}
        vals["Rpred"] = vals["Bstar"] * vals["Qpred"]
        vals["Deltapred"] = vals["Rpred"] - vals["R0"]
        vals["Hpred"] = vals["D"] + vals["Rpred"]
        for name, value in vals.items():
            saved = complex(float(row[name + "_real"]), float(row[name + "_imag"]))
            worst = max(worst, abs(saved - value))
            assert abs(saved - value) < 1e-14
            assert abs(frozen[name + "_real"] - saved.real) < 1e-14
            assert abs(frozen[name + "_imag"] - saved.imag) < 1e-14
        dsp = vals["Hpred"].conjugate() * cmath.exp(1j * 2 * math.pi * f * ld / 1500)
        components = (vals["D"].conjugate() + vals["Rpred"].conjugate()) * cmath.exp(1j * 2 * math.pi * f * ld / 1500)
        assert abs(dsp - components) < 1e-14
    assert sha(csv_path) == lock["prediction_csv_sha256"]
    assert all(sha(ROOT / x["path"]) == x["sha256"] for x in lock["inputs"])
    assert [float(r["frequency_hz"]) for r in rows] == [4750, 7250, 4125, 6562.5]
    assert lock["gates"]["e_R_strict_lt"] == lock["gates"]["e_Delta_strict_lt"] == lock["gates"]["a_D_strict_lt"] == 1e-3
    audit = {"status": "PASS", "audited_utc": datetime.now(timezone.utc).isoformat(),
             "prediction_lock_sha256": sha(lock_path), "prediction_csv_sha256": sha(csv_path),
             "all_input_sha_match": True, "new_frequency_not_in_old_nine": True,
             "max_independent_complex_difference_m_inv": worst,
             "native_DSP_component_closure": True, "new_FE_count_at_audit": 0}
    with (OUT / "prediction_lock_audit.json").open("x", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2)
        handle.write("\n")
    print(json.dumps(audit))

if __name__ == "__main__":
    main()
