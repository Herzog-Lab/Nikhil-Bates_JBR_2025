# -*- coding: utf-8 -*-
"""
Created on Thu Sep 25 18:33:09 2025

@author: Nikhil KL
"""

#!/usr/bin/env python3
"""
daily_peaks_findpeaks_py39.py

Detect *daily* peaks from embryonic-day time series using scipy.signal.find_peaks
with adaptive prominence and a minimum distance ~ one day.

- First column: embryonic day (float; e.g., 13.5, 13.667, ...)
- Remaining columns: samples (bioluminescence)
- Output CSV columns: Sample, Day, Peak_Embryonic_Day, Peak_Value, Series_Mean, Prominence_Threshold
- Output PNG: overview with each sample as a subplot; peak markers overlaid

CONFIG: edit the block below if running in Spyder/IDLE. No CLI required.
"""

import math
import os
import re
from typing import List, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks   # pip/conda: scipy

# ===================== CONFIG (edit these in Spyder) ===================== #
INPUT_CSV      = "Detrended_data.csv"       # path to your input file
OUTPUT_CSV     = "test.csv"       # path to save results
OUTPUT_PNG     = "daily_peaks_overview.png"     # path to save plot

TIME_COL_NAME: Optional[str] = None             # None = use first column; or e.g., "Embryonic day"
DAY_OFFSET     = 0.0                             # shift day boundary (0.5 = 12h)
MAX_SUBPLOT_COLS = 3
PNG_DPI        = 300

# Peak detector tuning (adaptive defaults work well; tweak if needed)
PROMINENCE_FACTOR = 0.00001   # fraction of robust amplitude (P95 - P05)
MIN_DISTANCE_FACTOR = 0.0 # fraction of samples-per-day (0.85*SPD ensures ≤1/day)
INTERPOLATE_MISSING = True # interpolate NaNs before peak finding
# ======================================================================== #


def sanitize(name: str) -> str:
    return re.sub(r"[^\w\-.]+", "_", str(name)).strip("_")[:150]


def load_table(path: str, time_col_name: Optional[str]) -> (pd.DataFrame, str, List[str]):
    df = pd.read_csv(path)
    if df.empty:
        raise ValueError("Input CSV is empty.")
    time_col = df.columns[0] if time_col_name is None else time_col_name
    if time_col not in df.columns:
        raise ValueError(f"Time column '{time_col}' not found. Available: {list(df.columns)}")
    df[time_col] = pd.to_numeric(df[time_col], errors="coerce")
    data_cols = [c for c in df.columns if c != time_col]
    if not data_cols:
        raise ValueError("No sample columns found.")
    return df, time_col, data_cols


def infer_samples_per_day(time_vals: np.ndarray) -> int:
    """Infer samples per (embryonic) day from median time step."""
    diffs = np.diff(time_vals.astype(float))
    diffs = diffs[np.isfinite(diffs)]
    if diffs.size == 0:
        return 6
    dt = float(np.median(diffs)) if np.isfinite(np.median(diffs)) else 1.0/6.0
    spd = int(round(1.0 / max(dt, 1e-6)))
    return max(spd, 2)


def adaptive_prominence(y: np.ndarray, factor: float) -> float:
    """Use robust amplitude (P95 - P05) to set a prominence threshold."""
    yy = y[np.isfinite(y)]
    if yy.size == 0:
        return 0.0
    p5, p95 = np.percentile(yy, [5, 95])
    rob_amp = max(p95 - p5, 0.0)
    if rob_amp <= 0:
        # fallback to std
        rob_amp = float(np.nanstd(yy))
    return max(factor * rob_amp, 0.0)


def find_daily_peaks_for_series(
    time_arr: np.ndarray,
    y_arr: np.ndarray,
    day_offset: float,
    spd: int,
    prominence_factor: float,
    min_distance_factor: float,
    interpolate_missing: bool
) -> pd.DataFrame:
    """Find local peaks with find_peaks, then keep at most one per day (best peak)."""
    # Interpolate missing if requested (does not alter original for reporting)
    if interpolate_missing:
        y_interp = pd.Series(y_arr).interpolate(limit_direction="both").to_numpy()
    else:
        y_interp = y_arr.copy()

    prom = adaptive_prominence(y_interp, prominence_factor)
    min_distance = max(1, int(round(min_distance_factor * spd)))

    # Peak detection (on interpolated series to avoid NaNs)
    peaks_idx, props = find_peaks(y_interp, prominence=prom if prom > 0 else None, distance=min_distance)

    # Assign each peak to an integer day bucket
    day_bucket = np.floor(time_arr + day_offset).astype("float")
    unique_days = np.array(sorted({int(d) for d in day_bucket[np.isfinite(day_bucket)]}))

    rows = []
    for d in unique_days:
        # peak indices belonging to this day
        mask = (day_bucket[peaks_idx] == d)
        idxs = peaks_idx[mask]
        if idxs.size == 0:
            rows.append({"Day": d, "Peak_Embryonic_Day": np.nan, "Peak_Value": np.nan,
                         "Prominence_Threshold": prom})
            continue
        # choose the peak with highest height within the day
        best = int(idxs[np.nanargmax(y_interp[idxs])])
        peak_time = float(time_arr[best])
        # prefer the original value if finite; else use interpolated
        peak_val = float(y_arr[best]) if np.isfinite(y_arr[best]) else float(y_interp[best])
        rows.append({"Day": d, "Peak_Embryonic_Day": peak_time, "Peak_Value": peak_val,
                     "Prominence_Threshold": prom})

    return pd.DataFrame(rows)


def process_file(
    input_csv: str,
    output_csv: str,
    output_png: str,
    time_col_name: Optional[str],
    day_offset: float,
    max_subplot_cols: int,
    png_dpi: int,
    prominence_factor: float,
    min_distance_factor: float,
    interpolate_missing: bool
) -> pd.DataFrame:
    df, time_col, data_cols = load_table(input_csv, time_col_name)
    t = df[time_col].to_numpy()
    spd = infer_samples_per_day(t)

    all_rows = []
    for col in data_cols:
        y = pd.to_numeric(df[col], errors="coerce").to_numpy()
        peaks_df = find_daily_peaks_for_series(
            t, y, day_offset, spd, prominence_factor, min_distance_factor, interpolate_missing
        )
        # add sample and series mean
        series_mean = float(np.nanmean(y)) if np.isfinite(y).any() else np.nan
        peaks_df.insert(0, "Sample", col)
        peaks_df["Series_Mean"] = series_mean
        all_rows.append(peaks_df)

    out = pd.concat(all_rows, ignore_index=True)
    out = out[["Sample", "Day", "Peak_Embryonic_Day", "Peak_Value", "Series_Mean", "Prominence_Threshold"]]
    out.sort_values(["Sample", "Day"], inplace=True, kind="stable")
    out.to_csv(output_csv, index=False)

    # Plot overview
    plot_overview(df, time_col, data_cols, out, output_png, max_subplot_cols, png_dpi)

    return out


def plot_overview(
    df: pd.DataFrame,
    time_col: str,
    data_cols: List[str],
    peaks_table: pd.DataFrame,
    out_png: str,
    max_cols: int,
    dpi: int
) -> None:
    n = len(data_cols)
    ncols = min(max_cols, n) if n > 0 else 1
    nrows = int(math.ceil(float(n) / float(ncols))) if n > 0 else 1

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(5 * ncols, 3.2 * nrows))
    if n == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = np.array([axes])
    elif ncols == 1:
        axes = np.array([[ax] for ax in axes])

    t = df[time_col].to_numpy()

    for i, col in enumerate(data_cols):
        r, c = divmod(i, ncols)
        ax = axes[r, c]

        y = pd.to_numeric(df[col], errors="coerce").to_numpy()
        ax.plot(t, y, linewidth=1.5)
        sub = peaks_table[peaks_table["Sample"] == col].dropna(subset=["Peak_Embryonic_Day"])
        if not sub.empty:
            ax.scatter(sub["Peak_Embryonic_Day"], sub["Peak_Value"], s=40, zorder=3)

        ax.set_title(str(col))
        if r == nrows - 1:
            ax.set_xlabel("Embryonic day")
        ax.set_ylabel("Bioluminescence (a.u.)")

    # Hide unused axes
    for j in range(n, nrows * ncols):
        r, c = divmod(j, ncols)
        axes[r, c].axis("off")

    plt.tight_layout()
    plt.show()
    os.makedirs(os.path.dirname(out_png) or ".", exist_ok=True)
    plt.savefig(out_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    results = process_file(
        input_csv=INPUT_CSV,
        output_csv=OUTPUT_CSV,
        output_png=OUTPUT_PNG,
        time_col_name=TIME_COL_NAME,
        day_offset=DAY_OFFSET,
        max_subplot_cols=MAX_SUBPLOT_COLS,
        png_dpi=PNG_DPI,
        prominence_factor=PROMINENCE_FACTOR,
        min_distance_factor=MIN_DISTANCE_FACTOR,
        interpolate_missing=INTERPOLATE_MISSING,
    )
    print("Saved CSV:", os.path.abspath(OUTPUT_CSV))
    print("Saved PNG:", os.path.abspath(OUTPUT_PNG))
