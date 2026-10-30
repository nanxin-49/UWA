"""Read-only audit of the incomplete 4031.25 Hz P22 pilot evidence."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.io import loadmat

from audit_two_dimensional_s5_p22_saved import (
    SOURCE_SHA,
    TAU_REF,
    close,
    complex_array,
    scalar,
    sha256,
    source_at,
)


def verify_manifest(root: Path, relative_manifest: str, prefix: str) -> int:
    lines = (root / relative_manifest).read_text(encoding="utf-8").splitlines()
    seen: set[str] = set()
    for line in lines:
        digest, separator, name = line.partition("  ")
        assert separator and len(digest) == 64 and name.startswith(prefix), line
        name = name[len(prefix) :]
        assert name not in seen and ".." not in Path(name).parts
        assert "cash" not in Path(name).parts
        seen.add(name)
        target = root / name
        assert target.is_file() and sha256(target) == digest, name
    return len(lines)


def job_status(jobdir: Path) -> tuple[int, dict]:
    assert int((jobdir / "exit_code.txt").read_text().strip()) == 0
    seconds = int((jobdir / "process_elapsed_seconds.txt").read_text().strip())
    assert 0 < seconds < 900
    assert (jobdir / "stderr.log").stat().st_size == 0
    assert (jobdir / "stdout.log").is_file()
    receipt = json.loads((jobdir / "receipt.json").read_text(encoding="utf-8"))
    assert receipt["status"] == "COMPLETE_UNAUDITED"
    for line in (jobdir / "output_sha256.txt").read_text().splitlines():
        digest, separator, path = line.partition("  ")
        name = Path(path).name
        assert separator and name in ("result.mat", "receipt.json")
        assert sha256(jobdir / name) == digest
    return seconds, receipt


def audit(original: Path, continuation: Path, output: Path) -> dict:
    assert not output.exists(), "Audit output must be new"
    original_count = verify_manifest(
        original, "audit/remote_files.sha256", "./"
    )
    continuation_count = verify_manifest(
        continuation, "audit/continuation_files.sha256", ""
    )
    source_path = original / "input" / "frozen_source.mat"
    assert sha256(source_path) == SOURCE_SHA
    frozen = loadmat(source_path, squeeze_me=True, struct_as_record=False)["src"]
    assert int(frozen.n_original) == 984 and float(frozen.z_tx_m) == 100

    prediction_dir = continuation / "prediction"
    prediction = json.loads(
        (prediction_dir / "prediction.json").read_text(encoding="utf-8")
    )
    prediction_sha = sha256(prediction_dir / "prediction.json")
    assert prediction_sha == (prediction_dir / "prediction_sha256.txt").read_text().strip()
    assert prediction["status"] == "FROZEN_BEFORE_HALFSTEP_BIE"
    assert prediction["source_sha256"] == SOURCE_SHA
    assert prediction["frequency_hz"] == 4031.25
    assert prediction["tau_ref_s"] == TAU_REF
    F = [complex(re, im) for re, im in zip(prediction["F_real"], prediction["F_imag"])]
    R4000 = complex(prediction["R4000_real"], prediction["R4000_imag"])
    R4062 = complex(prediction["R4062p5_real"], prediction["R4062p5_imag"])
    predicted = complex(prediction["native_real"], prediction["native_imag"])
    reconstructed = F[2] * (R4000 / F[0] + R4062 / F[1]) / 2
    predictor_formula_error = abs(reconstructed - predicted) / abs(predicted)
    assert predictor_formula_error <= 1e-14

    jobs = continuation / "jobs"
    retry_dir = jobs / "freeze_retry1"
    assert int((retry_dir / "exit_code.txt").read_text().strip()) == 0
    assert int((retry_dir / "process_elapsed_seconds.txt").read_text().strip()) == 12
    assert (retry_dir / "stderr.log").stat().st_size == 0
    assert "outer_shell_exit_code=1" in (
        retry_dir / "wrapper_check_failure.txt"
    ).read_text()
    frequency_root = jobs / "f4031p25"
    assert prediction["created_utc"] < (
        frequency_root / "main96" / "start_utc.txt"
    ).read_text().strip()

    model_dir = frequency_root / "model"
    model_seconds, model_receipt = job_status(model_dir)
    assert model_receipt["job_id"] == "model"
    assert model_receipt["frequency_hz"] == 4031.25
    assert prediction_sha in (
        model_dir / "prediction_sha256_at_start.txt"
    ).read_text()
    with h5py.File(model_dir / "result.mat", "r") as model:
        assert scalar(model, "cfg/frequency_hz") == 4031.25
        assert scalar(model, "tau_ref") == TAU_REF
        direct = scalar(model, "H_direct_native")
        pe = scalar(model, "H_pe_reflect_native")
        fk = scalar(model, "H_fk_reflect_native")
        assert close(direct, source_at(frozen, 0, 3, 4031.25), 1e-12)
        assert close(scalar(model, "H_flat_native"), -source_at(frozen, 0, -3, 4031.25), 1e-12)
        assert close(scalar(model, "H_pe_total_native"), direct + pe)
        assert close(scalar(model, "H_fk_total_native"), direct + fk)
        phase = np.exp(1j * 2 * math.pi * 4031.25 * TAU_REF)
        assert close(scalar(model, "H_pe_dsp"), np.conj(direct + pe) * phase)
        assert close(scalar(model, "H_fk_dsp"), np.conj(direct + fk) * phase)
        assert np.array_equal(
            complex_array(model, "src/coeff"), np.asarray(frozen.coeff).reshape(-1)
        )
        assert np.array_equal(
            np.asarray(model["src/kx"][()]).reshape(-1), np.asarray(frozen.kx).reshape(-1)
        )
        fk400 = scalar(model, "fk400/receiver_field")
        fk800 = scalar(model, "fk800/receiver_field")
        fk_quadrature = abs(fk400 - fk800) / abs(fk800)
        assert close(fk_quadrature, model_receipt["fk_400_800_relative"])

    rows: list[dict] = []
    fields: list[complex] = []
    seconds: dict[str, int] = {"model": model_seconds}
    for layer in ("main96", "self64"):
        jobdir = frequency_root / layer
        seconds[layer], receipt = job_status(jobdir)
        assert receipt["job_id"] == layer and receipt["frequency_hz"] == 4031.25
        assert prediction_sha in (jobdir / "prediction_sha256_at_start.txt").read_text()
        with h5py.File(jobdir / "result.mat", "r") as result:
            boundary = result["boundary"]
            assert scalar(result, "pre/frequency_hz") == 4031.25
            assert scalar(result, "pre/tau_ref_s") == TAU_REF
            assert scalar(result, "pre/receiver_x_m") == 0
            assert scalar(result, "pre/receiver_z_m") == 3
            assert scalar(result, "pre/surface_A_m") == 0.01
            assert scalar(result, "pre/surface_K_radpm") == 0.1
            assert scalar(result, "pre/surface_phi_rad") == 0
            assert scalar(boundary, "finite") == 1
            nodes = int(scalar(boundary, "node_count"))
            assert nodes == int(scalar(result, "pre/node_estimate"))
            field = scalar(result, "H_reflect_native")
            assert close(field, scalar(boundary, "receiver_field"), 1e-14)
            assert close(field, complex(receipt["receiver_native_real"],
                                        receipt["receiver_native_imag"]), 1e-14)
            assert close(scalar(result, "H_direct_native"), direct)
            assert close(scalar(result, "H_total_native"), direct + field)
            assert close(scalar(result, "H_reflect_dsp"), np.conj(field) * phase)
            assert close(scalar(result, "H_total_dsp"), np.conj(direct + field) * phase)
            assert np.array_equal(
                complex_array(result, "src/coeff"), np.asarray(frozen.coeff).reshape(-1)
            )
            linear = scalar(boundary, "linear_residual")
            offgrid = scalar(boundary, "offgrid_boundary_residual")
            assert math.isfinite(linear) and math.isfinite(offgrid)
        fields.append(field)
        rows.append({
            "layer": layer, "nodes": nodes, "native_real": field.real,
            "native_imag": field.imag, "linear_residual": linear,
            "offgrid_boundary_residual": offgrid, "process_seconds": seconds[layer],
        })

    for layer in ("space14", "window60", "window70"):
        assert not (frequency_root / layer / "result.mat").exists()
        assert not (frequency_root / layer / "receipt.json").exists()
    no_output = frequency_root / "group_remaining3"
    assert int((no_output / "exit_code.txt").read_text().strip()) == 0
    assert int((no_output / "process_elapsed_seconds.txt").read_text().strip()) == 54
    assert (no_output / "stdout.log").stat().st_size == 0
    assert "FAIL_MISSING_LAYER_OUTPUTS" in (
        no_output / "wrapper_postcheck_failure.txt"
    ).read_text()
    license_failure = frequency_root / "group_remaining3_retry1"
    assert not (license_failure / "exit_code.txt").exists()
    assert not (license_failure / "process_elapsed_seconds.txt").exists()
    prompt = (license_failure / "stdout.log").read_text(encoding="utf-8")
    assert "MathWorks" in prompt and "登录失败" in prompt
    assert "GROUP_START" not in prompt
    status = json.loads(
        (continuation / "audit" / "continuation_status.json").read_text()
    )
    assert status["known_prior_process_seconds"] == 1270
    assert status["unobserved_attempt_upper_seconds"] == 510
    assert status["conservative_total_seconds"] == 1780 <= 1800
    assert 829 + 12 + model_seconds + seconds["main96"] + seconds["self64"] + 54 == 1270

    main = fields[0]
    partial_sensitivity = abs(fields[1] - main) / abs(main)
    metrics = {
        "status": "PARTIAL_NUMERICS_VALIDATED_FULL_BIE_REFERENCE_BLOCKED_LICENSE",
        "frequency_hz": 4031.25,
        "original_file_sha_verified_count": original_count,
        "continuation_file_sha_verified_count": continuation_count,
        "prediction_json_sha256": prediction_sha,
        "prediction_formula_relative_error": predictor_formula_error,
        "prediction_relative_to_main96_error": abs(predicted - main) / abs(main),
        "two_layer_change_relative": partial_sensitivity,
        "full_five_layer_reference_gate": "NOT_EVALUABLE",
        "pe_relative_to_main96": abs(pe - main) / abs(main),
        "fk_relative_to_main96": abs(fk - main) / abs(main),
        "fk_quadrature_relative": fk_quadrature,
        "known_process_seconds": 1270,
        "unobserved_attempt_max_charge_seconds": 510,
        "conservative_total_seconds": 1780,
        "budget_seconds": 1800,
        "time_recovery_requires_new_authorization": True,
        "remaining_layers": ["space14", "window60", "window70"],
        "scope": "P22 two-frequency pilot incomplete; no qualified CIR/PDP/BER",
    }
    response = {
        "frequency_hz": 4031.25,
        "direct_native_real": direct.real, "direct_native_imag": direct.imag,
        "predicted_reflect_real": predicted.real,
        "predicted_reflect_imag": predicted.imag,
        "bie_main_reflect_real": main.real,
        "bie_main_reflect_imag": main.imag,
        "bie_self_reflect_real": fields[1].real,
        "bie_self_reflect_imag": fields[1].imag,
        "pe_reflect_real": pe.real, "pe_reflect_imag": pe.imag,
        "fk_reflect_real": fk.real, "fk_reflect_imag": fk.imag,
        "tau_ref_s": TAU_REF,
    }
    output.mkdir(parents=True)
    (output / "audit.json").write_text(
        json.dumps(metrics, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    for name, records in (("layers.csv", rows), ("response.csv", [response])):
        with (output / name).open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    with (output / "output_sha256.txt").open("w", encoding="utf-8") as stream:
        for name in ("audit.json", "layers.csv", "response.csv"):
            stream.write(f"{sha256(output / name)}  {name}\n")
    return metrics


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit("usage: partial_saved.py ORIGINAL_MIRROR CONTINUATION_MIRROR NEW_OUTPUT")
    print(json.dumps(audit(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))))
