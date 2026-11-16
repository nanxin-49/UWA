"""Strict runtime configuration serialization/validation; no numerical solve."""
import json
import os
import re
import sys
from pathlib import Path


def uint(text):
    if not re.fullmatch(r"[0-9]+", text):
        raise ValueError("empty or non-integer memory value")
    value = int(text)
    if value > 2**53 - 1:
        raise ValueError("memory value exceeds exact integer range")
    return value


def validate(record):
    expected = {"RAM_baseline_bytes", "starting_available_bytes", "initial_available_bytes",
                "reserve_fraction", "cumulative_time_cap", "cost_prediction_gate"}
    if not isinstance(record, dict) or set(record) != expected:
        raise ValueError("invalid runtime schema")
    for key in ("RAM_baseline_bytes", "starting_available_bytes", "initial_available_bytes"):
        if type(record[key]) is not int or not 0 <= record[key] <= 2**53 - 1:
            raise ValueError("invalid integer field: " + key)
    if record["RAM_baseline_bytes"] != max(record["starting_available_bytes"], record["initial_available_bytes"]):
        raise ValueError("baseline must equal max(starting, initial)")
    if record["reserve_fraction"] != 0.30 or record["cumulative_time_cap"] is not False or record["cost_prediction_gate"] is not False:
        raise ValueError("execution policy changed")
    return record


def main():
    mode, filename = sys.argv[1:3]
    path = Path(filename)
    if mode == "write":
        baseline, starting, initial = map(uint, sys.argv[3:6])
        record = validate(dict(RAM_baseline_bytes=baseline, starting_available_bytes=starting,
                               initial_available_bytes=initial, reserve_fraction=0.30,
                               cumulative_time_cap=False, cost_prediction_gate=False))
        if path.exists():
            raise FileExistsError("refusing to replace existing runtime configuration")
        temp = path.with_suffix(path.suffix + ".tmp")
        with temp.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(record, stream, allow_nan=False)
            stream.write("\n")
        validate(json.loads(temp.read_text(encoding="utf-8")))
        os.replace(temp, path)
    elif mode != "verify":
        raise ValueError("unknown runtime operation")
    validate(json.loads(path.read_text(encoding="utf-8")))


if __name__ == "__main__":
    main()
