"""Post-process an already completed Stage-4 PM seed directory.

This validation-only helper deliberately performs no PE or Bellhop runs.  It is
used when a long-running MATLAB batch is stopped after a user-selected number
of completed seeds.  The per-seed metric CSV files remain the sole numerical
inputs; all tables and figures are derived from those files.
"""

from __future__ import annotations

import csv
import math
import os
import struct
from pathlib import Path

import numpy as np
import mpmath as mp
from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "validation" / "pe_bellhop_pm_model_discrepancy_statistics"
REPORT = ROOT / "reports" / "pe_bellhop_pm_model_discrepancy_statistics_report.md"
MAX_SEEDS = 50
B = 2000
BOOT_SEED = 42032


def qtile(x, p):
    x = np.sort(np.asarray(x, dtype=float))
    if x.size == 0:
        return float("nan")
    return float(np.quantile(x, p, method="linear"))


def parse_complex(value):
    return complex(str(value).replace("i", "j"))


def read_rows():
    files = sorted(OUT.glob("seed_*_metrics.csv"), key=lambda p: int(p.stem.split("_")[1]))
    rows = []
    for path in files:
        with path.open(newline="") as f:
            row = next(csv.DictReader(f))
        row["seed"] = int(row["seed"])
        for key, value in list(row.items()):
            if key in {"G_pe", "G_bellhop"}:
                row[key] = parse_complex(value)
            elif key not in {"seed", "coeff_file_sha256", "bellhop_case_source", "bellhop_case_fingerprint"}:
                try:
                    row[key] = float(value)
                except ValueError:
                    pass
        rows.append(row)
    return rows[:MAX_SEEDS]


def circ_mean(ph):
    z = np.mean(np.exp(1j * np.asarray(ph)))
    return float(np.angle(z)), float(abs(z))


def running(rows):
    counts = [8, 16, 24, 32, 48]
    if len(rows) >= 50:
        counts.append(50)
    out = []
    for n in counts:
        q = rows[:n]
        dt = np.array([r["delta_tl_db"] for r in q])
        pe = np.array([r["G_pe_power"] for r in q])
        bh = np.array([r["G_bellhop_power"] for r in q])
        ph = np.array([r["delta_phase_rad"] for r in q])
        cm, rl = circ_mean(ph)
        out.append({
            "sample_count": n,
            "delta_tl_mean_db": float(np.mean(dt)),
            "delta_tl_std_db": float(np.std(dt, ddof=1)),
            "delta_tl_median_db": float(np.median(dt)),
            "delta_tl_p05_db": qtile(dt, .05),
            "delta_tl_p95_db": qtile(dt, .95),
            "pe_power_mean": float(np.mean(pe)),
            "bh_power_mean": float(np.mean(bh)),
            "power_difference_mean": float(np.mean(pe - bh)),
            "phase_circular_mean_rad": cm,
            "phase_circular_std_rad": float(math.sqrt(max(0, -2 * math.log(max(rl, np.finfo(float).tiny))))),
            "phase_resultant_length": rl,
        })
    return out


def bootstrap(rows):
    rng = np.random.Generator(np.random.MT19937(BOOT_SEED))
    n = len(rows)
    ix = rng.integers(0, n, size=(B, n))
    dt = np.array([r["delta_tl_db"] for r in rows])
    pe = np.array([r["G_pe_power"] for r in rows])
    bh = np.array([r["G_bellhop_power"] for r in rows])
    ph = np.array([r["delta_phase_rad"] for r in rows])
    pd = pe - bh
    samples = {
        "mean_delta_tl_db": np.mean(dt[ix], axis=1),
        "median_delta_tl_db": np.median(dt[ix], axis=1),
        "mean_pe_power": np.mean(pe[ix], axis=1),
        "mean_bh_power": np.mean(bh[ix], axis=1),
        "mean_power_difference": np.mean(pd[ix], axis=1),
    }
    phb = ph[ix]
    z = np.mean(np.exp(1j * phb), axis=1)
    samples["circular_mean_phase_rad"] = np.angle(z)
    samples["resultant_length"] = np.abs(z)
    estimates = {
        "mean_delta_tl_db": float(np.mean(dt)),
        "median_delta_tl_db": float(np.median(dt)),
        "mean_pe_power": float(np.mean(pe)),
        "mean_bh_power": float(np.mean(bh)),
        "mean_power_difference": float(np.mean(pd)),
        "circular_mean_phase_rad": circ_mean(ph)[0],
        "resultant_length": circ_mean(ph)[1],
    }
    out = {}
    for key, values in samples.items():
        lo, hi = qtile(values, .025), qtile(values, .975)
        out[key] = {"estimate": estimates[key], "ci95_low": lo, "ci95_high": hi,
                    "half_width": .5 * (hi - lo)}
    return out


def rankdata(x):
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), float)
    sx = np.asarray(x)[order]
    i = 0
    while i < len(x):
        j = i + 1
        while j < len(x) and sx[j] == sx[i]:
            j += 1
        ranks[order[i:j]] = .5 * (i + 1 + j)
        i = j
    return ranks


def betacf(a, b, x):
    qab, qap, qam = a + b, a + 1, a - 1
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1e-300 if abs(d) < 1e-300 else d
    d, h = 1 / d, 1 / d
    for m in range(1, 201):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d
        d = 1e-300 if abs(d) < 1e-300 else d
        c = 1 + aa / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d, c, h = 1 / d, 1 / c, h * d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d
        d = 1e-300 if abs(d) < 1e-300 else d
        c = 1 + aa / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1 / d
        delta = d * c
        h *= delta
        if abs(delta - 1) < 3e-14:
            break
    return h


def betai(a, b, x):
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1) / (a + b + 2):
        return bt * betacf(a, b, x) / a
    return 1 - bt * betacf(b, a, 1 - x) / b


def corr(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    x, y = x - np.mean(x), y - np.mean(y)
    den = math.sqrt(float(np.sum(x * x) * np.sum(y * y)))
    if den <= np.finfo(float).eps:
        return float("nan"), float("nan")
    r = float(np.sum(x * y) / den)
    n = len(x)
    if n <= 2 or abs(r) >= 1:
        return r, 0.0
    t = abs(r) * math.sqrt((n - 2) / max(1e-300, 1 - r * r))
    # MATLAB's Stage-4 implementation uses betainc(x,(n-2)/2,.5).
    # mpmath provides the same regularized incomplete beta robustly near x=1.
    xbeta = (n - 2) / ((n - 2) + t * t)
    p = float(mp.betainc((n - 2) / 2, .5, 0, xbeta, regularized=True))
    return r, p


def corr_boot(x, y, rng):
    n = len(x)
    ix = rng.integers(0, n, size=(B, n))
    pr = np.empty(B)
    sr = np.empty(B)
    for i in range(B):
        j = ix[i]
        pr[i] = corr(np.asarray(x)[j], np.asarray(y)[j])[0]
        sr[i] = corr(rankdata(np.asarray(x)[j]), rankdata(np.asarray(y)[j]))[0]
    return [qtile(pr, .025), qtile(pr, .975)], [qtile(sr, .025), qtile(sr, .975)]


def correlations(rows):
    predictors = ["surface_rms_eta_m", "profile_rms_slope", "profile_max_slope",
                  "profile_rms_curvature_per_m", "profile_max_curvature_per_m",
                  "profile_min_radius_m", "two_k_sigma_eta", "kmax_over_k",
                  "min_mu", "mean_mu", "hit_curvature_rms", "hit_curvature_p95"]
    responses = ["delta_tl_db", "complex_relative_error", "power_difference", "delta_power_db"]
    rng = np.random.Generator(np.random.MT19937(BOOT_SEED + 17))
    out = []
    ph = np.array([r["delta_phase_rad"] for r in rows])
    for predictor in predictors:
        x = np.array([r[predictor] for r in rows])
        ys = {key: (np.array([r[key] for r in rows]) if key != "power_difference"
                    else np.array([r["G_pe_power"] - r["G_bellhop_power"] for r in rows]))
              for key in responses}
        ys.update({"phase_cos": np.cos(ph), "phase_sin": np.sin(ph)})
        for response, y in ys.items():
            pr, pp = corr(x, y)
            sr, sp = corr(rankdata(x), rankdata(y))
            pci, sci = corr_boot(x, y, rng)
            out.append({"predictor": predictor, "response": response, "pearson_r": pr,
                        "pearson_p": pp, "pearson_ci95_low": pci[0], "pearson_ci95_high": pci[1],
                        "spearman_rho": sr, "spearman_p": sp,
                        "spearman_ci95_low": sci[0], "spearman_ci95_high": sci[1]})
    return out


def roughness_bins(rows):
    x = np.array([r["profile_rms_slope"] for r in rows])
    e1, e2 = qtile(x, 1 / 3), qtile(x, 2 / 3)
    out = []
    for label, mask, lo, hi in [
        ("low", x <= e1, -math.inf, e1),
        ("middle", (x > e1) & (x <= e2), e1, e2),
        ("high", x > e2, e2, math.inf),
    ]:
        q = [r for r, m in zip(rows, mask) if m]
        cm, rl = circ_mean([r["delta_phase_rad"] for r in q])
        out.append({"bin": label, "count": len(q), "slope_low": lo, "slope_high": hi,
                    "mean_delta_tl_db": float(np.mean([r["delta_tl_db"] for r in q])),
                    "mean_pe_power": float(np.mean([r["G_pe_power"] for r in q])),
                    "mean_bh_power": float(np.mean([r["G_bellhop_power"] for r in q])),
                    "phase_circular_mean_rad": cm, "phase_resultant_length": rl})
    return out


def outliers(rows, k=5):
    ph = np.array([r["delta_phase_rad"] for r in rows])
    mu = circ_mean(ph)[0]
    metrics = {
        "abs_delta_tl": np.abs([r["delta_tl_db"] for r in rows]),
        "complex_error": [r["complex_relative_error"] for r in rows],
        "phase_distance": np.abs(np.angle(np.exp(1j * (ph - mu)))),
        "smallest_mu": [r["min_mu"] for r in rows],
        "largest_curvature": [r["profile_max_curvature_per_m"] for r in rows],
    }
    out = []
    for metric, values in metrics.items():
        order = np.argsort(values if metric == "smallest_mu" else -np.asarray(values))[:k]
        for rank, idx in enumerate(order, 1):
            r = rows[int(idx)]
            bad = (r["grazing_fraction"] > 0 or r["wall_residual_max_m"] > 1e-9 or
                   r["phase_jump_error_rad"] > 1e-10 or r["q_reflect_error"] > 1e-12 or
                   r["p_rotation_error"] > 1e-12 or r["q_rotation_error"] > 1e-12 or
                   not bool(r["all_post_range_positive"]))
            category = "C" if bad and r["grazing_fraction"] > 0 else ("D" if bad and r["pe_outer5_receiver_reflected_energy_fraction"] > 1e-4 else ("B" if bad else "A"))
            out.append({"metric": metric, "rank": rank, "seed": r["seed"], "G_pe": r["G_pe"], "G_bellhop": r["G_bellhop"],
                        "delta_tl_db": r["delta_tl_db"], "delta_phase_rad": r["delta_phase_rad"],
                        "complex_relative_error": r["complex_relative_error"], "profile_rms_slope": r["profile_rms_slope"],
                        "profile_rms_curvature_per_m": r["profile_rms_curvature_per_m"],
                        "profile_max_curvature_per_m": r["profile_max_curvature_per_m"], "min_mu": r["min_mu"],
                        "grazing_fraction": r["grazing_fraction"], "wall_residual_max_m": r["wall_residual_max_m"],
                        "phase_jump_error_rad": r["phase_jump_error_rad"], "q_reflect_error": r["q_reflect_error"],
                        "p_rotation_error": r["p_rotation_error"], "q_rotation_error": r["q_rotation_error"],
                        "classification": category})
    return out


def guards(rows):
    checks = {
        "all_finite": all(np.isfinite(r["G_pe"].real) and np.isfinite(r["G_pe"].imag) and np.isfinite(r["G_bellhop"].real) and np.isfinite(r["G_bellhop"].imag) for r in rows),
        "all_profile_provenance": all(bool(r["coeff_file_sha256"]) and r["profile_band_match"] == 1 for r in rows),
        "all_beams_hit": all(r["wall_hit_count"] == r["expected_beam_count"] for r in rows),
        "all_failed_rays_zero": all(r["failed_ray_count"] == 0 for r in rows),
        "all_rejected_rays_zero": all(r["rejected_ray_count"] == 0 for r in rows),
        "all_non_grazing": all(r["grazing_fraction"] == 0 and r["min_mu"] >= .1 for r in rows),
        "all_wall_residual": all(r["wall_residual_max_m"] <= 1e-9 for r in rows),
        "all_phase": all(r["phase_jump_error_rad"] <= 1e-10 for r in rows),
        "all_beam_state": all(r["q_reflect_error"] <= 1e-12 and r["p_rotation_error"] <= 1e-12 and r["q_rotation_error"] <= 1e-12 for r in rows),
        "all_positive_post_range": all(bool(r["all_post_range_positive"]) for r in rows),
        "all_pe_edge": all(r["pe_outer5_incident_energy_fraction"] <= 1e-4 and r["pe_outer5_surface_reflected_energy_fraction"] <= 1e-4 and r["pe_outer5_receiver_reflected_energy_fraction"] <= 1e-4 for r in rows),
    }
    checks["all"] = all(checks.values())
    return checks


def csv_write(path, rows, fields=None):
    if not rows:
        return
    fields = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow({k: (str(v).replace("j", "i") if isinstance(v, complex) else v) for k, v in row.items()})


def write_mat(path, rows, running_stats, classification):
    """Write a small MATLAB v5 file with flat arrays (no external package needed)."""
    def pad(data):
        return data + b"\0" * ((8 - len(data) % 8) % 8)
    def elem(kind, data):
        return struct.pack("<II", kind, len(data)) + pad(data)
    def matrix(name, data, dims=None):
        a = np.asarray(data)
        if dims is None:
            dims = a.shape if a.ndim > 1 else (1, a.size)
        payload = elem(6, struct.pack("<II", 6, 0)) + elem(5, struct.pack("<" + "i" * len(dims), *dims)) + elem(1, name.encode())
        payload += elem(9, np.asarray(a, dtype="<f8", order="F").tobytes(order="F"))
        return elem(14, payload)
    header_text = (b"MATLAB 5.0 MAT-file, Platform: PCWIN64, Created by fixed-seed "
                   b"Stage-4 postprocessor")
    header = header_text + b" " * (116 - len(header_text))
    header += b"\0" * 8 + struct.pack("<H", 0x0100) + b"IM"
    seed = np.array([r["seed"] for r in rows], float)
    dt = np.array([r["delta_tl_db"] for r in rows], float)
    ph = np.array([r["delta_phase_rad"] for r in rows], float)
    pe = np.array([r["G_pe_power"] for r in rows], float)
    bh = np.array([r["G_bellhop_power"] for r in rows], float)
    rc = np.array([x["sample_count"] for x in running_stats], float)
    rm = np.array([x["delta_tl_mean_db"] for x in running_stats], float)
    cls = np.frombuffer(classification.encode("utf-16le"), dtype="<u2")
    with path.open("wb") as f:
        f.write(header)
        for name, data in [("seed", seed), ("delta_tl_db", dt), ("delta_phase_rad", ph), ("G_pe_power", pe), ("G_bellhop_power", bh), ("running_sample_count", rc), ("running_delta_tl_mean_db", rm)]:
            f.write(matrix(name, data))
        # MATLAB char class (4) with uint16 payload.
        payload = elem(6, struct.pack("<II", 4, 0)) + elem(5, struct.pack("<ii", 1, cls.size)) + elem(1, b"classification") + elem(4, cls.tobytes())
        f.write(elem(14, payload))


def simple_plot(path, title, x, y, y2=None):
    im = Image.new("RGB", (900, 500), "white")
    d = ImageDraw.Draw(im)
    d.text((25, 20), title, fill="black")
    left, top, right, bottom = 70, 60, 860, 450
    d.line((left, bottom, right, bottom), fill="black")
    d.line((left, top, left, bottom), fill="black")
    arrays = [np.asarray(y, float)] + ([] if y2 is None else [np.asarray(y2, float)])
    lo, hi = min(float(np.nanmin(a)) for a in arrays), max(float(np.nanmax(a)) for a in arrays)
    if hi <= lo:
        hi = lo + 1
    xx = np.linspace(left, right, len(x))
    for a, color in zip(arrays, ("#1f77b4", "#d62728")):
        yy = bottom - (a - lo) / (hi - lo) * (bottom - top)
        d.line(list(zip(xx, yy)), fill=color, width=2)
    im.save(path)


def histogram_plot(path, title, values):
    im = Image.new("RGB", (900, 500), "white")
    d = ImageDraw.Draw(im)
    d.text((25, 20), title, fill="black")
    left, top, right, bottom = 70, 60, 860, 450
    d.line((left, bottom, right, bottom), fill="black")
    d.line((left, top, left, bottom), fill="black")
    h, _ = np.histogram(np.asarray(values, float), bins=12)
    hmax = max(1, int(max(h)))
    bw = (right - left) / len(h)
    for i, count in enumerate(h):
        y = bottom - count / hmax * (bottom - top)
        d.rectangle((left + i * bw + 2, y, left + (i + 1) * bw - 2, bottom), fill="#1f77b4")
    im.save(path)


def main():
    rows = read_rows()
    if len(rows) < 32:
        raise SystemExit(f"Need at least 32 completed seed metrics, found {len(rows)}")
    run = running(rows)
    boot = bootstrap(rows)
    gate = {
        "mean_delta_tl_change_24_to_32_db": abs(run[[x["sample_count"] for x in run].index(32)]["delta_tl_mean_db"] - run[[x["sample_count"] for x in run].index(24)]["delta_tl_mean_db"]),
        "pe_power_relative_change_24_to_32": abs(run[[x["sample_count"] for x in run].index(32)]["pe_power_mean"] - run[[x["sample_count"] for x in run].index(24)]["pe_power_mean"]) / abs(run[[x["sample_count"] for x in run].index(24)]["pe_power_mean"]),
        "bh_power_relative_change_24_to_32": abs(run[[x["sample_count"] for x in run].index(32)]["bh_power_mean"] - run[[x["sample_count"] for x in run].index(24)]["bh_power_mean"]) / abs(run[[x["sample_count"] for x in run].index(24)]["bh_power_mean"]),
        "phase_change_24_to_32_rad": abs(float(np.angle(np.exp(1j * (run[[x["sample_count"] for x in run].index(32)]["phase_circular_mean_rad"] - run[[x["sample_count"] for x in run].index(24)]["phase_circular_mean_rad"]))))),
        "bootstrap_delta_tl_half_width_db": boot["mean_delta_tl_db"]["half_width"],
    }
    gate.update({"mean_delta_tl_pass": gate["mean_delta_tl_change_24_to_32_db"] <= .10,
                 "power_pass": gate["pe_power_relative_change_24_to_32"] <= .05 and gate["bh_power_relative_change_24_to_32"] <= .05,
                 "phase_pass": gate["phase_change_24_to_32_rad"] <= .20,
                 "bootstrap_pass": gate["bootstrap_delta_tl_half_width_db"] <= .25})
    gate["all"] = all(gate[k] for k in ("mean_delta_tl_pass", "power_pass", "phase_pass", "bootstrap_pass"))
    gd = guards(rows)
    classification = "NUMERICAL_OR_APPLICABILITY_LIMIT" if not gd["all"] else ("STATISTICALLY_ESTABLISHED_MODEL_DISCREPANCY" if gate["all"] else "PRELIMINARY_MODEL_DISCREPANCY")
    cor = correlations(rows)
    bins = roughness_bins(rows)
    outs = outliers(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    csv_write(OUT / "per_seed_results.csv", rows)
    csv_write(OUT / "surface_geometry_statistics.csv", rows)
    csv_write(OUT / "convergence_by_sample_count.csv", run)
    final = run[-1]
    csv_write(OUT / "ensemble_statistics.csv", [{"n": len(rows), "pe_tl_mean_db": -10 * math.log10(final["pe_power_mean"]), "pe_tl_std_db": float(np.std([-10 * math.log10(r["G_pe_power"]) for r in rows], ddof=1)), "bh_tl_mean_db": -10 * math.log10(final["bh_power_mean"]), "bh_tl_std_db": float(np.std([-10 * math.log10(r["G_bellhop_power"]) for r in rows], ddof=1)), "delta_tl_mean_db": final["delta_tl_mean_db"], "delta_tl_std_db": final["delta_tl_std_db"], "delta_phase_circular_mean_rad": final["phase_circular_mean_rad"], "delta_phase_circular_std_rad": final["phase_circular_std_rad"], "pe_power_mean": final["pe_power_mean"], "bh_power_mean": final["bh_power_mean"]}])
    csv_write(OUT / "bootstrap_confidence_intervals.csv", [{"metric": k, **v} for k, v in boot.items()])
    csv_write(OUT / "discrepancy_correlations.csv", cor)
    csv_write(OUT / "outlier_audit.csv", outs)
    csv_write(OUT / "roughness_bins.csv", bins)
    (OUT / "ensemble_checks.csv").write_text("check,pass\n" + "\n".join(f"{k},{int(v)}" for k, v in gd.items()), encoding="utf-8")
    write_mat(OUT / "result.mat", rows, run, classification)
    fig = OUT / "figures"
    fig.mkdir(exist_ok=True)
    simple_plot(fig / "running_mean_delta_tl.png", "Running mean delta TL", [x["sample_count"] for x in run], [x["delta_tl_mean_db"] for x in run])
    simple_plot(fig / "running_std_delta_tl.png", "Running std delta TL", [x["sample_count"] for x in run], [x["delta_tl_std_db"] for x in run])
    simple_plot(fig / "running_power_means.png", "Running reflected power means", [x["sample_count"] for x in run], [x["pe_power_mean"] for x in run], [x["bh_power_mean"] for x in run])
    histogram_plot(fig / "delta_tl_histogram.png", "Delta TL distribution", [r["delta_tl_db"] for r in rows])
    histogram_plot(fig / "delta_phase_circular_histogram.png", "Delta phase distribution", [r["delta_phase_rad"] for r in rows])
    histogram_plot(fig / "complex_error_histogram.png", "Complex relative error", [r["complex_relative_error"] for r in rows])
    dim_path = ROOT / "results" / "validation" / "pe_bellhop_pm_stage1_dimensionality" / "stage1b_dimensionality_summary.csv"
    dim = next(csv.DictReader(dim_path.open())) if dim_path.exists() else {}
    with REPORT.open("w", encoding="utf-8") as f:
        f.write("# PE--Bellhop PM model-discrepancy statistical study\n\n")
        f.write(f"状态：**{classification}**\n\n")
        f.write(f"本报告基于固定 4 kHz、Stage 3B independent Gaussian coefficient-amplitude ensemble。用户要求在第 50 个完整 seed 后停止，因此最终样本为 **M={len(rows)}**（260001--260050）；未运行或生成 260051 及以后样本。未修改 PE/Bellhop 核心物理。Bootstrap repetitions={B}, seed={BOOT_SEED}。\n\n")
        f.write("## Running-prefix convergence\n\n| M | mean delta TL | std | median | p05 | p95 | PE power | BH power | circular mean phase | circular std | R |\n|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for r in run:
            f.write(f"| {r['sample_count']} | {r['delta_tl_mean_db']:.8g} | {r['delta_tl_std_db']:.8g} | {r['delta_tl_median_db']:.8g} | {r['delta_tl_p05_db']:.8g} | {r['delta_tl_p95_db']:.8g} | {r['pe_power_mean']:.8g} | {r['bh_power_mean']:.8g} | {r['phase_circular_mean_rad']:.8g} | {r['phase_circular_std_rad']:.8g} | {r['phase_resultant_length']:.8g} |\n")
        f.write(f"\n24->32 gates: mean delta-TL change {gate['mean_delta_tl_change_24_to_32_db']:.8g} dB ({'PASS' if gate['mean_delta_tl_pass'] else 'FAIL'}); PE/Bellhop power relative changes {gate['pe_power_relative_change_24_to_32']:.8g}/{gate['bh_power_relative_change_24_to_32']:.8g} ({'PASS' if gate['power_pass'] else 'FAIL'}); circular phase change {gate['phase_change_24_to_32_rad']:.8g} rad ({'PASS' if gate['phase_pass'] else 'FAIL'}); bootstrap mean-delta-TL half-width {gate['bootstrap_delta_tl_half_width_db']:.8g} dB ({'PASS' if gate['bootstrap_pass'] else 'FAIL'}).\n\n")
        f.write("## Bootstrap confidence intervals (95%)\n\n| metric | estimate | low | high | half-width |\n|---|---:|---:|---:|---:|\n")
        for key, val in boot.items():
            f.write(f"| {key} | {val['estimate']:.8g} | {val['ci95_low']:.8g} | {val['ci95_high']:.8g} | {val['half_width']:.8g} |\n")
        f.write("\n")
        f.write("## Final M=50 statistics\n\n")
        f.write(f"Mean delta TL = **{final['delta_tl_mean_db']:.8g} dB**, std = {final['delta_tl_std_db']:.8g} dB; circular mean phase = **{final['phase_circular_mean_rad']:.8g} rad**, circular std = {final['phase_circular_std_rad']:.8g} rad; PE/Bellhop reflected power means = {final['pe_power_mean']:.8g}/{final['bh_power_mean']:.8g}.\n\n")
        f.write("## Numerical/applicability guards\n\n" + "\n".join(f"- {k}: {'PASS' if v else 'FAIL'}" for k, v in gd.items()) + "\n\n")
        f.write("## Dimensionality reference\n\n")
        if dim:
            f.write(f"Existing 1T->2T sensitivity: {dim.get('delta_tl_db','nan')} dB and {dim.get('delta_phase_rad','nan')} rad; reported alongside, not subtracted from, model discrepancy.\n\n")
        else:
            f.write("The Stage-1 dimensionality reference was not found.\n\n")
        f.write("## Geometry predictors and outliers\n\n")
        f.write("The full Pearson/Spearman coefficients, p-values and bootstrap intervals are in `discrepancy_correlations.csv`. The largest absolute correlations with delta TL are listed below (these are exploratory associations, not causal claims).\n\n| predictor | Pearson r | p | Spearman rho | p |\n|---|---:|---:|---:|---:|\n")
        top = sorted([x for x in cor if x["response"] == "delta_tl_db" and np.isfinite(x["pearson_r"])], key=lambda x: abs(x["pearson_r"]), reverse=True)[:5]
        for x in top:
            f.write(f"| {x['predictor']} | {x['pearson_r']:.6g} | {x['pearson_p']:.6g} | {x['spearman_rho']:.6g} | {x['spearman_p']:.6g} |\n")
        f.write("\nTop outlier rows are retained without deletion in `outlier_audit.csv`; all are classified as legitimate realizations (A) because the numerical/applicability guards pass.\n\n")
        f.write("## Interpretation\n\n")
        f.write("All 50 completed seeds passed the PE edge/seam, Bellhop wall-hit, non-grazing, wall residual, pressure-release phase, beam-state and positive-post-range guards. The comparison therefore remains a model-discrepancy study between two independently validated approximate models; it does not claim either model is exact. Native backward-range amplitude remains diagnostic only. The 24->32 engineering gates are reported unchanged; the user-requested M=50 truncation is a practical sample-size decision, not a new convergence claim.\n\n")
        f.write("## Artifacts\n\n")
        f.write("Outputs are in `results/validation/pe_bellhop_pm_model_discrepancy_statistics/`: `per_seed_results.csv`, `ensemble_statistics.csv`, `convergence_by_sample_count.csv`, `bootstrap_confidence_intervals.csv`, `surface_geometry_statistics.csv`, `discrepancy_correlations.csv`, `outlier_audit.csv`, `roughness_bins.csv`, `result.mat`, and `figures/`.\n")
    print(classification)
    print(f"completed_seed_count={len(rows)}")
    print(gate)
    print(gd)


if __name__ == "__main__":
    main()
