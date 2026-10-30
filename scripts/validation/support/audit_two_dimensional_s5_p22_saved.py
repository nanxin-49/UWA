"""Independent read-only audit of saved P22 pilot MAT files (no acoustic solves)."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
import traceback
from pathlib import Path

import h5py
import numpy as np
from scipy.io import loadmat


LAYERS = ("main96", "self64", "space14", "window60", "window70")
TAU_REF = 96.875 / 1500
SOURCE_SHA = "8d60d3d9d6aa3aa16110ac33d81e0d193a61526aef2c0c2e9edd531cae531f55"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def scalar(group: h5py.Group, key: str) -> float | complex:
    value = np.asarray(group[key][()]).reshape(-1)
    assert value.size == 1, key
    if value.dtype.names == ("real", "imag"):
        return complex(float(value["real"][0]), float(value["imag"][0]))
    return float(value[0])


def complex_array(group: h5py.Group, key: str) -> np.ndarray:
    value = np.asarray(group[key][()]).reshape(-1)
    if value.dtype.names == ("real", "imag"):
        return value["real"].astype(float) + 1j * value["imag"].astype(float)
    return value.astype(complex)


def close(a: complex, b: complex, tolerance: float = 1e-12) -> bool:
    return abs(a - b) <= tolerance * max(1.0, abs(b))


def verify_manifest(mirror: Path) -> int:
    manifest = mirror / "audit" / "remote_files.sha256"
    assert manifest.is_file()
    rows = manifest.read_text(encoding="utf-8").splitlines()
    seen: set[str] = set()
    for line in rows:
        expected, marker, relative = line.partition("  ")
        assert marker and len(expected) == 64 and relative.startswith("./")
        name = relative[2:]
        assert name not in seen and ".." not in Path(name).parts and "cash" not in Path(name).parts
        seen.add(name)
        path = mirror / name
        assert path.is_file() and sha256(path) == expected, name
    return len(rows)


def source_at(src: object, x: float, z: float, frequency_hz: float) -> complex:
    kx = np.asarray(src.kx).reshape(-1)
    coeff = np.asarray(src.coeff).reshape(-1)
    k0 = 2 * math.pi * frequency_hz / 1500
    kz = np.sqrt((k0 * k0 - kx * kx).astype(complex))
    phase = np.exp(1j * (x - float(src.x0_m)) * kx) * np.exp(
        -1j * (z - float(src.z_tx_m)) * kz
    )
    return complex(np.sum(phase * coeff) / int(src.n_original))


def job_evidence(jobdir: Path) -> tuple[int, dict]:
    assert int((jobdir / "exit_code.txt").read_text().strip()) == 0
    seconds = int((jobdir / "process_elapsed_seconds.txt").read_text().strip())
    assert seconds > 0 and seconds < 900
    receipt = json.loads((jobdir / "receipt.json").read_text(encoding="utf-8"))
    assert (jobdir / "stdout.log").is_file() and (jobdir / "stderr.log").is_file()
    assert (jobdir / "result.mat").is_file()
    return seconds, receipt


def audit(mirror: Path, frequency_hz: float) -> tuple[dict, list[dict], dict]:
    assert frequency_hz in (4062.5, 4031.25)
    file_count = verify_manifest(mirror)
    source_path = mirror / "input" / "frozen_source.mat"
    assert sha256(source_path) == SOURCE_SHA
    frozen = loadmat(source_path, squeeze_me=True, struct_as_record=False)["src"]
    assert int(frozen.n_original) == 984 and float(frozen.z_tx_m) == 100
    label = "f4062p5" if frequency_hz == 4062.5 else "f4031p25"
    jobroot = mirror / "jobs" / label
    modeldir = jobroot / "model"
    model_seconds, model_receipt = job_evidence(modeldir)
    assert model_receipt["frequency_hz"] == frequency_hz
    with h5py.File(modeldir / "result.mat", "r") as model:
        assert scalar(model, "cfg/frequency_hz") == frequency_hz
        assert scalar(model, "cfg/step_m") == 0.05
        assert scalar(model, "cfg/nx") == 984
        assert scalar(model, "tau_ref") == TAU_REF
        direct = scalar(model, "H_direct_native")
        flat = scalar(model, "H_flat_native")
        pe_reflect = scalar(model, "H_pe_reflect_native")
        fk_reflect = scalar(model, "H_fk_reflect_native")
        assert close(direct, source_at(frozen, 0, 3, frequency_hz), 1e-12)
        assert close(flat, -source_at(frozen, 0, -3, frequency_hz), 1e-12)
        assert close(pe_reflect, scalar(model, "H_pe_total_native") - direct)
        assert close(fk_reflect, scalar(model, "H_fk_total_native") - direct)
        phase = np.exp(1j * 2 * math.pi * frequency_hz * TAU_REF)
        assert close(scalar(model, "H_pe_dsp"), np.conj(direct + pe_reflect) * phase)
        assert close(scalar(model, "H_fk_dsp"), np.conj(direct + fk_reflect) * phase)
        fk400 = scalar(model, "fk400/receiver_field")
        fk800 = scalar(model, "fk800/receiver_field")
        assert close(fk_reflect, fk800)
        fk_quadrature = abs(fk400 - fk800) / abs(fk800)
        assert close(fk_quadrature, model_receipt["fk_400_800_relative"])
        assert scalar(model, "pe/frequency_hz") == frequency_hz
        assert close(scalar(model, "pe/dx_m"), 192.1875 / 984)
        ix_rx = int(scalar(model, "pe/ix_rx")) - 1
        pe_saved = complex_array(model, "pe/reflected_native")
        assert close(pe_reflect, pe_saved[ix_rx])
        assert np.array_equal(complex_array(model, "src/coeff"), np.asarray(frozen.coeff).reshape(-1))
        assert np.array_equal(np.asarray(model["src/kx"][()]).reshape(-1), np.asarray(frozen.kx).reshape(-1))
    rows: list[dict] = []
    fields: list[complex] = []
    process_seconds = model_seconds
    for layer in LAYERS:
        jobdir = jobroot / layer
        seconds, receipt = job_evidence(jobdir)
        process_seconds += seconds
        assert receipt["status"] == "COMPLETE_UNAUDITED"
        assert receipt["job_id"] == layer and receipt["frequency_hz"] == frequency_hz
        with h5py.File(jobdir / "result.mat", "r") as result:
            boundary = result["boundary"]
            assert scalar(boundary, "frequency_hz") == frequency_hz
            assert scalar(result, "pre/frequency_hz") == frequency_hz
            assert scalar(result, "pre/tau_ref_s") == TAU_REF
            assert scalar(result, "pre/receiver_x_m") == 0
            assert scalar(result, "pre/receiver_z_m") == 3
            assert scalar(result, "pre/surface_A_m") == 0.01
            assert scalar(result, "pre/surface_K_radpm") == 0.1
            assert scalar(result, "pre/surface_phi_rad") == 0
            assert scalar(boundary, "finite") == 1
            assert scalar(boundary, "receiver_x_m") == 0
            assert scalar(boundary, "receiver_z_m") == 3
            nodes = int(scalar(boundary, "node_count"))
            assert nodes == int(scalar(result, "pre/node_estimate"))
            field = scalar(boundary, "receiver_field")
            assert close(field, scalar(result, "H_reflect_native"), 1e-14)
            assert close(field, complex(receipt["receiver_native_real"], receipt["receiver_native_imag"]), 1e-14)
            assert close(scalar(result, "H_total_native"), direct + field)
            assert close(scalar(result, "H_reflect_dsp"), np.conj(field) * phase)
            assert close(scalar(result, "H_total_dsp"), np.conj(direct + field) * phase)
            assert np.array_equal(complex_array(result, "src/coeff"), np.asarray(frozen.coeff).reshape(-1))
            linear = scalar(boundary, "linear_residual")
            offgrid = scalar(boundary, "offgrid_boundary_residual")
            assert math.isfinite(linear) and math.isfinite(offgrid)
            assert close(linear, receipt["linear_residual"], 1e-14)
            assert close(offgrid, receipt["offgrid_boundary_residual"], 1e-14)
        fields.append(field)
        rows.append({"layer": layer, "nodes": nodes, "native_real": field.real,
                     "native_imag": field.imag, "linear_residual": linear,
                     "offgrid_residual": offgrid, "process_seconds": seconds})
    reference = fields[0]
    sensitivity = max(abs(value - reference) for value in fields[1:]) / abs(reference)
    assert sensitivity <= 1e-4
    for row, value in zip(rows, fields):
        row["change_from_main_relative"] = abs(value - reference) / abs(reference)
    pe_model = abs(pe_reflect - reference) / abs(reference)
    fk_model = abs(fk_reflect - reference) / abs(reference)
    fk_resolution = sensitivity + fk_quadrature
    metrics = {
        "status": "BIE_EMPIRICAL_REFERENCE_PASS_MODEL_LIMITED",
        "frequency_hz": frequency_hz,
        "mirror_file_sha_verified_count": file_count,
        "reference_sensitivity_relative": sensitivity,
        "reference_gate_1e_4": True,
        "pe_model_relative_to_bie": pe_model,
        "fk_model_relative_to_bie": fk_model,
        "fk_quadrature_relative": fk_quadrature,
        "fk_empirical_resolution_relative": fk_resolution,
        "fk_model_resolved_over_5x": fk_model > 5 * fk_resolution,
        "pe_model_resolved_over_5x": None,
        "pe_resolution_limitation": "Independent PE step/grid sensitivity at this frequency not yet measured",
        "process_seconds_model_plus_five_bie": process_seconds,
        "source_sha256": SOURCE_SHA,
        "scope": "P22 single weak-sine frequency; no CIR/PDP/BER",
    }
    response = {
        "frequency_hz": frequency_hz,
        "H_direct_native_real": direct.real, "H_direct_native_imag": direct.imag,
        "H_bie_reflect_native_real": reference.real, "H_bie_reflect_native_imag": reference.imag,
        "H_pe_reflect_native_real": pe_reflect.real, "H_pe_reflect_native_imag": pe_reflect.imag,
        "H_fk_reflect_native_real": fk_reflect.real, "H_fk_reflect_native_imag": fk_reflect.imag,
        "H_bie_total_native_real": (direct + reference).real,
        "H_bie_total_native_imag": (direct + reference).imag,
        "H_bie_total_dsp_real": (np.conj(direct + reference) * phase).real,
        "H_bie_total_dsp_imag": (np.conj(direct + reference) * phase).imag,
        "tau_ref_s": TAU_REF,
    }
    return metrics, rows, response


def main() -> None:
    assert len(sys.argv) == 4, "usage: audit_saved.py MIRROR FREQ OUTPUT_DIR"
    mirror = Path(sys.argv[1]).resolve(strict=True)
    frequency_hz = float(sys.argv[2])
    output = Path(sys.argv[3]).resolve()
    assert not output.exists(), "Refuse to overwrite audit"
    output.mkdir(parents=True)
    try:
        metrics, layers, response = audit(mirror, frequency_hz)
        (output / "audit.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
        for name, records in (("layers.csv", layers), ("response.csv", [response])):
            with (output / name).open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(records[0]))
                writer.writeheader()
                writer.writerows(records)
        with (output / "output_sha256.txt").open("w", encoding="utf-8") as stream:
            for name in ("audit.json", "layers.csv", "response.csv"):
                stream.write(f"{sha256(output / name)}  {name}\n")
        print(json.dumps({"status": metrics["status"],
                          "sensitivity": metrics["reference_sensitivity_relative"],
                          "fk_model": metrics["fk_model_relative_to_bie"],
                          "pe_model": metrics["pe_model_relative_to_bie"],
                          "output": str(output)}))
    except Exception:
        (output / "failure.txt").write_text(traceback.format_exc(), encoding="utf-8")
        raise


if __name__ == "__main__":
    main()
