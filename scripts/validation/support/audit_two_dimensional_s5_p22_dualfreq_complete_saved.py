"""Independent saved-data closeout for the two-frequency P22 physical pilot.

No MATLAB process or acoustic solve is started. The input mirror is assembled
from SHA-verified remote snapshots; the original snapshots remain unchanged.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path

from audit_two_dimensional_s5_p22_saved import SOURCE_SHA, TAU_REF, audit, sha256


def verify_transfer(manifest: Path, root: Path) -> int:
    lines = manifest.read_text(encoding="utf-8").splitlines()
    assert lines
    for line in lines:
        digest, separator, relative = line.partition("  ./")
        assert separator and len(digest) == 64 and ".." not in Path(relative).parts
        assert "cash" not in Path(relative).parts
        target = root / relative
        assert target.is_file() and sha256(target) == digest, relative
    return len(lines)


def complex_fields(row: dict, prefix: str) -> complex:
    return complex(float(row[f"{prefix}_real"]), float(row[f"{prefix}_imag"]))


def read_kv(path: Path) -> dict[str, str]:
    return dict(line.split("=", 1) for line in path.read_text().splitlines())


def main() -> None:
    assert len(sys.argv) == 3, "usage: audit_complete_saved.py RUN_DIR NEW_OUTPUT"
    run = Path(sys.argv[1]).resolve(strict=True)
    output = Path(sys.argv[2]).resolve()
    assert not output.exists(), "Refuse to overwrite audit"
    mirror = run / "audit_mirror"
    addendum = run / "remote_addendum"
    transfer_count = verify_transfer(run / "remote_addendum.sha256", addendum)
    transfer_receipt = json.loads((run / "transfer_audit.json").read_text())
    assert transfer_receipt["status"] == "PASS"
    assert transfer_receipt["manifest_rows"] == transfer_count == 53

    base_metrics, _, base_response = audit(mirror, 4062.5)
    half_metrics, layers, half_response = audit(mirror, 4031.25)
    assert base_metrics["reference_gate_1e_4"]
    assert half_metrics["reference_gate_1e_4"]
    assert base_metrics["source_sha256"] == half_metrics["source_sha256"] == SOURCE_SHA
    assert base_response["tau_ref_s"] == half_response["tau_ref_s"] == TAU_REF

    prediction_dir = mirror / "prediction"
    prediction_path = prediction_dir / "prediction.json"
    prediction_sha = sha256(prediction_path)
    assert prediction_sha == (prediction_dir / "prediction_sha256.txt").read_text().strip()
    prediction = json.loads(prediction_path.read_text())
    assert prediction["status"] == "FROZEN_BEFORE_HALFSTEP_BIE"
    assert prediction["source_sha256"] == SOURCE_SHA
    assert prediction["frequency_hz"] == 4031.25
    assert prediction["tau_ref_s"] == TAU_REF
    reference_csv = mirror / "input" / "fields_core.csv"
    assert sha256(reference_csv) == prediction["reference_csv_sha256"]
    with reference_csv.open(newline="", encoding="utf-8") as stream:
        matches = [row for row in csv.DictReader(stream)
                   if row["case_name"] == "weak_low_K" and row["method"] == "BIE"
                   and row["receiver_id"] == "P22" and float(row["frequency_hz"]) == 4000]
    assert len(matches) == 1
    r4000 = complex_fields(matches[0], "u_native")
    r4062 = complex_fields(base_response, "H_bie_reflect_native")
    assert r4000 == complex_fields(prediction, "R4000")
    assert r4062 == complex_fields(prediction, "R4062p5")
    base_layers = ("main96", "self64", "space14", "window60", "window70")
    for layer, digest in zip(base_layers, prediction["base_layer_sha256"], strict=True):
        assert sha256(mirror / "jobs" / "f4062p5" / layer / "result.mat") == digest
    flat = [complex(re, im) for re, im in zip(
        prediction["F_real"], prediction["F_imag"], strict=True)]
    predicted = complex_fields(prediction, "native")
    reconstructed = flat[2] * (r4000 / flat[0] + r4062 / flat[1]) / 2
    formula_error = abs(predicted - reconstructed) / abs(predicted)
    assert formula_error <= 1e-14
    first_half_bie_start = (mirror / "jobs/f4031p25/main96/start_utc.txt").read_text().strip()
    assert prediction["created_utc"] < first_half_bie_start
    actual = complex_fields(half_response, "H_bie_reflect_native")
    prediction_error = abs(predicted - actual) / abs(actual)
    magnitude_error = abs(abs(predicted) - abs(actual)) / abs(actual)
    phase_error_rad = math.atan2((predicted / actual).imag, (predicted / actual).real)
    native_to_dsp = complex_fields(half_response, "H_bie_total_dsp")
    total_native = complex_fields(half_response, "H_bie_total_native")
    expected_dsp = total_native.conjugate() * complex(
        math.cos(2 * math.pi * 4031.25 * TAU_REF),
        math.sin(2 * math.pi * 4031.25 * TAU_REF),
    )
    assert abs(native_to_dsp - expected_dsp) <= 1e-12 * max(1, abs(expected_dsp))

    probe = addendum / "jobs/license_probe_20261031_01"
    assert (probe / "qualification.txt").read_text().strip() == "PASS"
    assert (probe / "exit_code.txt").read_text().strip() == "0"
    probe_seconds = int((probe / "process_elapsed_seconds.txt").read_text())
    assert 0 < probe_seconds <= 90
    probe_log = (probe / "stdout.log").read_text(encoding="utf-8")
    assert "LICENSE_PROBE_OK" in probe_log and "R2023b" in probe_log
    assert (probe / "stderr.log").stat().st_size == 0
    controller = addendum / "jobs/continuation_20261031_01"
    assert (controller / "exit_code.txt").read_text().strip() == "0"
    assert "ALL_LAYERS_COMPLETE" in (controller / "controller.log").read_text()
    spent = probe_seconds
    for layer, upper in (("space14", 350), ("window60", 300), ("window70", 350)):
        root = addendum / "jobs/f4031p25" / layer
        attempt = root / "attempt_20261031_01"
        preflight = read_kv(attempt / "preflight.txt")
        assert preflight["layer"] == layer
        assert preflight["prediction_sha256"] == prediction_sha
        assert int(preflight["spent_before_seconds"]) == spent
        assert int(preflight["process_upper_seconds"]) == upper
        assert int(preflight["active_matlab_before"]) == 0
        available = int(preflight["ram_available_bytes"])
        total = int(preflight["ram_total_bytes"])
        increment = int(preflight["estimated_increment_bytes"])
        assert (available - increment) * 10 >= total * 3
        assert int(preflight["disk_available_bytes"]) >= 10 * 2**30
        assert int(preflight["tmp_available_bytes"]) >= 10 * 2**30
        assert (attempt / "exit_code.txt").read_text().strip() == "0"
        seconds = int((attempt / "process_elapsed_seconds.txt").read_text())
        assert 0 < seconds <= upper and seconds < 900
        assert seconds == next(row["process_seconds"] for row in layers if row["layer"] == layer)
        assert prediction["created_utc"] < (attempt / "start_utc.txt").read_text().strip()
        assert (attempt / "stderr.log").stat().st_size == 0
        assert "登录失败" not in (attempt / "stdout.log").read_text(encoding="utf-8")
        spent += seconds
    assert spent == int((controller / "new_budget_spent_seconds.txt").read_text())
    assert spent <= 1200

    metrics = {
        "status": "DUAL_FREQUENCY_P22_REFERENCE_PASS_MODEL_RESOLUTION_LIMITED",
        "frequency_hz": [4062.5, 4031.25],
        "remote_addendum_sha_verified_files": transfer_count,
        "combined_mirror_sha_verified_files": base_metrics["mirror_file_sha_verified_count"],
        "source_sha256": SOURCE_SHA,
        "prediction_json_sha256": prediction_sha,
        "prediction_frozen_before_halfstep_bie": True,
        "prediction_formula_relative_error": formula_error,
        "prediction_complex_relative_error": prediction_error,
        "prediction_magnitude_relative_error": magnitude_error,
        "prediction_phase_error_rad": phase_error_rad,
        "single_holdout_under_1e_3_diagnostic": prediction_error <= 1e-3,
        "full_band_65_129_frequency_gate": "NOT_RUN",
        "bie_reference_4062p5_relative": base_metrics["reference_sensitivity_relative"],
        "bie_reference_4031p25_relative": half_metrics["reference_sensitivity_relative"],
        "bie_reference_limit": 1e-4,
        "fk_model_4031p25_relative": half_metrics["fk_model_relative_to_bie"],
        "fk_model_resolved_over_5x": half_metrics["fk_model_resolved_over_5x"],
        "pe_model_4031p25_relative": half_metrics["pe_model_relative_to_bie"],
        "pe_model_resolution": "NOT_QUALIFIED_NO_INDEPENDENT_PE_STEP_SENSITIVITY",
        "prior_budget_conservative_seconds": 1780,
        "prior_budget_limit_seconds": 1800,
        "new_budget_actual_seconds": spent,
        "new_budget_limit_seconds": 1200,
        "scope": "P22 weak-sine two-frequency pilot only; no CIR/PDP/BER",
    }
    output.mkdir(parents=True)
    result = output / "audit.json"
    result.write_text(json.dumps(metrics, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    with (output / "output_sha256.txt").open("w", encoding="utf-8") as stream:
        stream.write(f"{hashlib.sha256(result.read_bytes()).hexdigest()}  audit.json\n")
    print(json.dumps({"status": metrics["status"], "prediction_error": prediction_error,
                      "bie_reference_4031p25": metrics["bie_reference_4031p25_relative"],
                      "new_budget_seconds": spent, "output": str(output)}))


if __name__ == "__main__":
    main()
