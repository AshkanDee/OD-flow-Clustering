"""Tidy HDBSCAN pipeline for per-city OD-flow clustering.

This script refactors the notebook-style workflow into reusable functions:
1) Load and clean OD data from RDS files.
2) Build feature matrices and z-score scaling.
3) Optional diagnostics (kNN-distance summaries, SNN overlap, HDBSCAN bounds).
4) HDBSCAN grid evaluation over `min_samples`.
5) Final per-city HDBSCAN fit and cluster assignment.
6) Post-fit validation/summary (cluster sizes, dominant cluster, direction vectors).
7) Optional NSGA-II multi-objective search on subsamples.
"""

from __future__ import annotations

from dataclasses import dataclass
from os import PathLike
import os
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence
import importlib

import numpy as np
import pandas as pd
from pyproj import Transformer
from sklearn.neighbors import NearestNeighbors


NEEDED_COLS = ("start_loc_lon", "start_loc_lat", "dest_loc_lon", "dest_loc_lat")
FEATURE_COLS = ("x_o", "y_o", "dx", "dy")




def _require_package(package_name: str, purpose: str):
    if importlib.util.find_spec(package_name) is None:
        raise ModuleNotFoundError(
            f"Missing optional dependency '{package_name}' required for {purpose}. "
            f"Install it first (e.g., call install_required_packages_colab())."
        )
    return importlib.import_module(package_name)


@dataclass
class CityPrepResult:
    city_id: str
    data: pd.DataFrame
    info: Dict[str, int | str]


@dataclass
class CityFeatureResult:
    city_id: str
    city_clean: pd.DataFrame
    X: np.ndarray
    X_scaled: np.ndarray
    ok_mask: np.ndarray


def read_rds_as_df(rds_path: Path) -> pd.DataFrame:
    if not rds_path.exists():
        raise FileNotFoundError(f"File not found: {rds_path}")
    pyreadr = _require_package("pyreadr", "reading .rds files")
    res = pyreadr.read_r(str(rds_path))
    if not res:
        raise ValueError(f"No objects found in RDS: {rds_path}")
    obj = next(iter(res.values()))
    if not isinstance(obj, pd.DataFrame):
        raise TypeError(f"Expected dataframe in {rds_path}, got {type(obj)}")
    return obj


def prep_city_od(
    city_id: str,
    rds_path: Path,
    needed_cols: Sequence[str] = NEEDED_COLS,
    crs_in: int = 4326,
    crs_out: int = 25832,
    drop_zero_len: bool = True,
) -> CityPrepResult:
    dt = read_rds_as_df(rds_path).copy()

    missing = [c for c in needed_cols if c not in dt.columns]
    if missing:
        raise ValueError(f"{city_id}: missing columns {missing}")

    n_before = len(dt)
    dt = dt.dropna(subset=list(needed_cols))
    n_after_na = len(dt)

    dt = dt[
        dt["start_loc_lat"].between(-90, 90)
        & dt["dest_loc_lat"].between(-90, 90)
        & dt["start_loc_lon"].between(-180, 180)
        & dt["dest_loc_lon"].between(-180, 180)
    ].copy()
    n_after_range = len(dt)

    transformer = Transformer.from_crs(f"EPSG:{crs_in}", f"EPSG:{crs_out}", always_xy=True)
    dt["x_o"], dt["y_o"] = transformer.transform(dt["start_loc_lon"].to_numpy(), dt["start_loc_lat"].to_numpy())
    dt["x_d"], dt["y_d"] = transformer.transform(dt["dest_loc_lon"].to_numpy(), dt["dest_loc_lat"].to_numpy())

    dt["dx"] = dt["x_d"] - dt["x_o"]
    dt["dy"] = dt["y_d"] - dt["y_o"]
    dt["angle"] = np.arctan2(dt["dy"].to_numpy(), dt["dx"].to_numpy())
    dt["len"] = np.hypot(dt["dx"].to_numpy(), dt["dy"].to_numpy())

    n_zero_len = int((dt["len"] == 0).sum())
    if drop_zero_len:
        dt = dt[dt["len"] > 0].copy()

    dt["city"] = city_id
    info = {
        "file": str(rds_path),
        "n_before": n_before,
        "n_after_na": n_after_na,
        "n_after_range": n_after_range,
        "n_zero_len": n_zero_len,
        "n_after_zero": len(dt),
    }
    return CityPrepResult(city_id=city_id, data=dt, info=info)


def build_features(city_id: str, city_df: pd.DataFrame, feature_cols: Sequence[str] = FEATURE_COLS) -> CityFeatureResult:
    missing = [c for c in feature_cols if c not in city_df.columns]
    if missing:
        raise ValueError(f"{city_id}: missing feature columns {missing}")

    X = city_df.loc[:, feature_cols].to_numpy()
    ok = np.isfinite(X).all(axis=1)
    X_ok = X[ok, :]
    mu = X_ok.mean(axis=0)
    sd = X_ok.std(axis=0, ddof=1)
    sd_safe = np.where(sd == 0, 1.0, sd)
    X_scaled = (X_ok - mu) / sd_safe

    return CityFeatureResult(
        city_id=city_id,
        city_clean=city_df.loc[ok, :].copy(),
        X=X_ok,
        X_scaled=X_scaled,
        ok_mask=ok,
    )


def cleaning_summary(prepared: Mapping[str, CityPrepResult]) -> pd.DataFrame:
    rows = []
    for city_id, res in prepared.items():
        rows.append(
            {
                "city": city_id,
                "file": res.info["file"],
                "n_before": res.info["n_before"],
                "n_after_na": res.info["n_after_na"],
                "n_after_range": res.info["n_after_range"],
                "zero_len_trips": res.info["n_zero_len"],
                "n_after_zero_len_removed": res.info["n_after_zero"],
            }
        )
    return pd.DataFrame(rows)


def knee_elbow(d_sorted: np.ndarray) -> float:
    d = np.asarray(d_sorted, dtype=float)
    n = d.size
    if n < 10:
        return float("nan")
    x = np.arange(1, n + 1, dtype=float)
    x1, y1 = 1.0, d[0]
    x2, y2 = float(n), d[-1]
    num = np.abs((y2 - y1) * x - (x2 - x1) * d + x2 * y1 - y2 * x1)
    den = np.sqrt((y2 - y1) ** 2 + (x2 - x1) ** 2)
    return float(d[int(np.argmax(num / den))])


def knn_distances(X_scaled: np.ndarray, k: int) -> np.ndarray:
    nn = NearestNeighbors(n_neighbors=k + 1, metric="euclidean").fit(X_scaled)
    dists, _ = nn.kneighbors(X_scaled, return_distance=True)
    return dists[:, k]


def knn_summary_one(
    X_scaled: np.ndarray,
    k: int,
    quantiles: Sequence[float] = (0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99),
) -> Dict[str, float | int]:
    d = knn_distances(X_scaled, k)
    out: Dict[str, float | int] = {"k": k, "n": int(d.size), "elbow_eps": knee_elbow(np.sort(d))}
    qs = np.quantile(d, q=quantiles)
    for q, val in zip(quantiles, qs):
        out[f"q{int(round(q * 100)):02d}"] = float(val)
    return out


def snn_sharednn_diagnostics(X_scaled: np.ndarray, k: int = 30, n_probe: int = 4000, seed: int = 1) -> Dict[str, float | int]:
    rng = np.random.default_rng(seed)
    n = X_scaled.shape[0]
    idx = rng.choice(n, size=min(n, n_probe), replace=False)

    nn = NearestNeighbors(n_neighbors=k + 1, metric="euclidean").fit(X_scaled)
    _, ids = nn.kneighbors(X_scaled, return_distance=True)
    ids = ids[:, 1 : k + 1]

    probe_sets = [set(ids[i].tolist()) for i in idx]
    checks = min(8000, len(probe_sets) * 2)

    shared = np.empty(checks, dtype=int)
    for t in range(checks):
        a = rng.integers(0, len(probe_sets))
        b = rng.integers(0, len(probe_sets))
        shared[t] = len(probe_sets[a].intersection(probe_sets[b]))

    return {
        "k": int(k),
        "probe_n": int(len(probe_sets)),
        "shared_mean": float(shared.mean()),
        "shared_median": float(np.median(shared)),
        "shared_q90": float(np.quantile(shared, 0.90)),
        "shared_q95": float(np.quantile(shared, 0.95)),
        "shared_q99": float(np.quantile(shared, 0.99)),
    }


def hdbscan_mcs_candidates(n_eff: int, n: int = 14, frac_min: float = 0.005, frac_max: float = 0.10) -> np.ndarray:
    lo = max(10, int(round(frac_min * n_eff)))
    hi = max(lo + 1, int(round(frac_max * n_eff)))
    vals = np.unique(np.round(np.exp(np.linspace(np.log(lo), np.log(hi), num=n))).astype(int))
    return np.sort(vals)


def evaluate_hdbscan_grid(
    features: Mapping[str, CityFeatureResult],
    min_samples_grid: Sequence[int],
    min_cluster_size: int,
    cluster_selection_method: str = "eom",
) -> pd.DataFrame:
    rows: List[Dict[str, float | int | str]] = []

    for city_id, res in features.items():
        X_scaled = res.X_scaled
        n_total = int(X_scaled.shape[0])
        for k in min_samples_grid:
            clusterer = _require_package("hdbscan", "HDBSCAN clustering").HDBSCAN(
                min_samples=int(k),
                min_cluster_size=int(min_cluster_size),
                metric="euclidean",
                cluster_selection_method=cluster_selection_method,
                core_dist_n_jobs=1,
            )
            labels = clusterer.fit_predict(X_scaled)
            is_noise = labels == -1
            labels_nonnoise = labels[~is_noise]
            n_clusters = int(np.unique(labels_nonnoise).size) if labels_nonnoise.size else 0

            if n_clusters:
                _, counts = np.unique(labels_nonnoise, return_counts=True)
                max_cluster_pct = 100.0 * float(counts.max()) / len(labels)
            else:
                max_cluster_pct = 0.0

            probs = getattr(clusterer, "probabilities_", None)
            mean_prob_nonnoise = float(probs[~is_noise].mean()) if probs is not None and (~is_noise).any() else np.nan

            rows.append(
                {
                    "city": city_id,
                    "algorithm": "HDBSCAN",
                    "minPts": int(k),
                    "eps": np.nan,
                    "n_clusters": n_clusters,
                    "noise_pct": 100.0 * float(is_noise.mean()),
                    "max_cluster_pct": max_cluster_pct,
                    "n_total": n_total,
                    "n_used": n_total,
                    "min_cluster_size": int(min_cluster_size),
                    "cluster_selection_method": cluster_selection_method,
                    "mean_prob_nonnoise": mean_prob_nonnoise,
                }
            )

    return pd.DataFrame(rows)


def fit_final_hdbscan(
    features: Mapping[str, CityFeatureResult],
    per_city_params: Mapping[str, Mapping[str, int]],
    min_cluster_size: int,
    cluster_selection_method: str = "eom",
) -> Dict[str, pd.DataFrame]:
    final_tables: Dict[str, pd.DataFrame] = {}

    for city_id, res in features.items():
        min_pts = int(per_city_params[city_id]["minPts"])
        clusterer = _require_package("hdbscan", "HDBSCAN clustering").HDBSCAN(
            min_samples=min_pts,
            min_cluster_size=int(min_cluster_size),
            metric="euclidean",
            cluster_selection_method=cluster_selection_method,
            core_dist_n_jobs=1,
        )
        labels = clusterer.fit_predict(res.X_scaled)
        labels_dbscan_style = np.where(labels == -1, 0, labels + 1)

        city_out = res.city_clean.copy()
        city_out["cluster"] = labels_dbscan_style
        final_tables[city_id] = city_out

    return final_tables


def verify_final_clustering(final_city_tables: Mapping[str, pd.DataFrame]) -> Dict[str, Dict[str, float | int | bool]]:
    checks: Dict[str, Dict[str, float | int | bool]] = {}
    for city_id, df in final_city_tables.items():
        if "cluster" not in df.columns:
            raise ValueError(f"{city_id}: missing `cluster` column")

        cl = df["cluster"]
        counts = cl.value_counts(dropna=False).sort_index()
        n_total = int(len(cl))
        n_noise = int((cl == 0).sum())
        noise_pct = 100.0 * n_noise / n_total if n_total > 0 else 0.0
        n_clusters = int(cl[cl != 0].nunique())
        max_cluster_size = int(counts.drop(index=0, errors="ignore").max()) if n_clusters > 0 else 0
        max_cluster_pct = 100.0 * max_cluster_size / n_total if n_clusters > 0 else 0.0

        checks[city_id] = {
            "n_total": n_total,
            "n_clusters": n_clusters,
            "noise_count": n_noise,
            "noise_pct": noise_pct,
            "max_cluster_pct": max_cluster_pct,
            "warn_high_noise": noise_pct > 90,
            "warn_no_clusters": n_clusters == 0,
            "warn_dominant_cluster": max_cluster_pct > 60,
        }
    return checks


def summarize_clusters(final_city_tables: Mapping[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    cluster_summaries: Dict[str, pd.DataFrame] = {}

    for city_id, city_dt in final_city_tables.items():
        if "cluster" not in city_dt.columns:
            raise ValueError(f"{city_id}: missing cluster column")

        has_len = "len" in city_dt.columns
        has_duration = "duration" in city_dt.columns
        has_dx = "dx" in city_dt.columns
        has_dy = "dy" in city_dt.columns

        missing_req = [c for c in ["x_o", "y_o"] if c not in city_dt.columns]
        if missing_req:
            raise ValueError(f"{city_id}: missing required columns: {missing_req}")

        df = city_dt.loc[city_dt["cluster"] != 0].copy()
        if df.empty:
            cluster_summaries[city_id] = pd.DataFrame()
            continue

        agg: Dict[str, tuple[str, str]] = {
            "size": ("cluster", "size"),
            "mean_x_o": ("x_o", "mean"),
            "mean_y_o": ("y_o", "mean"),
        }
        if has_dx and has_dy:
            agg["mean_dx"] = ("dx", "mean")
            agg["mean_dy"] = ("dy", "mean")
        if has_len:
            agg["mean_len"] = ("len", "mean")
        if has_duration:
            agg["mean_duration"] = ("duration", "mean")

        cluster_summary = (
            df.groupby("cluster", as_index=False)
            .agg(**agg)
            .sort_values("size", ascending=False)
            .reset_index(drop=True)
        )

        if not (has_dx and has_dy):
            cluster_summary["mean_dx"] = np.nan
            cluster_summary["mean_dy"] = np.nan
        if not has_len:
            cluster_summary["mean_len"] = np.nan
        if not has_duration:
            cluster_summary["mean_duration"] = np.nan

        cluster_summaries[city_id] = cluster_summary

    return cluster_summaries


def analyze_dominant_cluster(
    cluster_summaries: Mapping[str, pd.DataFrame],
    features: Mapping[str, CityFeatureResult],
) -> Dict[str, Dict[str, object]]:
    out: Dict[str, Dict[str, object]] = {}
    for city_id, cs in cluster_summaries.items():
        if isinstance(cs, pd.DataFrame) and (len(cs) > 0) and {"size", "cluster"}.issubset(cs.columns):
            dominant_cluster = int(cs.loc[cs["size"].idxmax(), "cluster"])
            cs_moo = cs.copy()
            k = int(len(cs_moo))
        else:
            dominant_cluster = None
            cs_moo = pd.DataFrame()
            k = 0

        x = features[city_id].X
        n_total = int(x.shape[0]) if hasattr(x, "shape") else int(len(x))

        out[city_id] = {
            "cs": cs_moo,
            "N_total": n_total,
            "dominant_cluster": dominant_cluster,
            "K": k,
        }
    return out


def build_trip_direction_vectors(final_city_tables: Mapping[str, pd.DataFrame]) -> Dict[str, np.ndarray]:
    trip_dirs_full: Dict[str, np.ndarray] = {}
    for city_id, df in final_city_tables.items():
        if ("dx" not in df.columns) or ("dy" not in df.columns):
            raise ValueError(f"{city_id}: dx/dy not found")

        dx = df["dx"].to_numpy(dtype=float)
        dy = df["dy"].to_numpy(dtype=float)
        norm = np.hypot(dx, dy)
        norm[norm == 0] = np.nan
        ux = dx / norm
        uy = dy / norm
        trip_dirs_full[city_id] = np.column_stack([ux, uy])
    return trip_dirs_full


def dir_cohesion_within_weighted(
    labels: np.ndarray,
    u: np.ndarray,
    noise_label: int = 0,
    min_cluster_size: int = 2,
    map_to_01: bool = True,
) -> float:
    labels = np.asarray(labels, dtype=int)
    u = np.asarray(u, dtype=float)
    if u.ndim != 2 or u.shape[1] != 2:
        raise ValueError("u must be Nx2")

    ux = u[:, 0]
    uy = u[:, 1]
    ok = (~np.isnan(labels)) & (labels != noise_label) & np.isfinite(ux) & np.isfinite(uy)
    if not np.any(ok):
        return 0.0

    labels = labels[ok]
    ux = ux[ok]
    uy = uy[ok]

    num = 0.0
    den = 0.0

    for cid in np.unique(labels):
        idx = labels == cid
        n = int(idx.sum())
        if n < int(min_cluster_size):
            continue

        mx = float(np.mean(ux[idx]))
        my = float(np.mean(uy[idx]))
        nm = np.hypot(mx, my)
        if (not np.isfinite(nm)) or (nm < 1e-12):
            continue

        mx /= nm
        my /= nm

        cos_i = np.clip(ux[idx] * mx + uy[idx] * my, -1.0, 1.0)
        coh_c = float(np.mean(cos_i))

        num += n * coh_c
        den += n

    if den == 0:
        return 0.0

    coh = num / den
    return (coh + 1.0) / 2.0 if map_to_01 else coh


def eval_solution(
    labels: np.ndarray,
    u_trips: np.ndarray,
    noise_label: int = 0,
    min_cluster_size: int = 2,
) -> Dict[str, float | int]:
    labels = np.asarray(labels, dtype=int)
    u_trips = np.asarray(u_trips, dtype=float)

    if u_trips.ndim != 2 or u_trips.shape[1] != 2:
        raise ValueError("u_trips must be Nx2")
    if labels.size <= 0:
        raise ValueError("empty labels")
    if u_trips.shape[0] != labels.size:
        raise ValueError("len(labels) must match nrow(u_trips)")

    clustered = (labels != noise_label) & np.isfinite(u_trips[:, 0]) & np.isfinite(u_trips[:, 1])
    n_clustered = int(clustered.sum())
    coverage = float(n_clustered / labels.size)
    k = int(np.unique(labels[clustered]).size) if n_clustered > 0 else 0

    coh = dir_cohesion_within_weighted(
        labels=labels,
        u=u_trips,
        noise_label=noise_label,
        min_cluster_size=min_cluster_size,
        map_to_01=True,
    )

    f1 = 0.0 if coverage == 0 else -coverage
    f2 = 0.0 if coverage == 0 else -float(coh)

    return {
        "f1_neg_coverage": f1,
        "f2_neg_dir_cohesion": f2,
        "coverage": coverage,
        "dir_cohesion": float(coh),
        "K": k,
        "clustered_n": n_clustered,
        "N_total": int(labels.size),
    }


def snap_scalar(v: float, allowed: Sequence[int]) -> int:
    vals = np.asarray(allowed, dtype=int)
    return int(vals[np.argmin(np.abs(vals - int(np.round(v))))])


def eval_hdbscan_params(
    min_cluster_size: int,
    min_samples: int,
    X_scaled: np.ndarray,
    u_trips: np.ndarray,
    min_cluster_size_for_cohesion: int = 2,
    cohesion_floor: Optional[float] = None,
    cluster_selection_method: str = "eom",
    core_dist_n_jobs: int = 1,
) -> np.ndarray:
    clusterer = _require_package("hdbscan", "HDBSCAN clustering").HDBSCAN(
        min_cluster_size=int(min_cluster_size),
        min_samples=int(min_samples),
        metric="euclidean",
        cluster_selection_method=cluster_selection_method,
        core_dist_n_jobs=int(core_dist_n_jobs),
    )
    labels = clusterer.fit_predict(X_scaled).astype(int)
    labels_dbscan_style = np.where(labels == -1, 0, labels + 1)

    res = eval_solution(
        labels=labels_dbscan_style,
        u_trips=u_trips,
        noise_label=0,
        min_cluster_size=min_cluster_size_for_cohesion,
    )

    if cohesion_floor is not None and (
        (not np.isfinite(res["dir_cohesion"])) or (res["dir_cohesion"] < float(cohesion_floor))
    ):
        return np.array([1.0, 1.0], dtype=float)

    return np.array([res["f1_neg_coverage"], res["f2_neg_dir_cohesion"]], dtype=float)


def make_objfun_city_cached(
    city_id: str,
    data_for_opt: Mapping[str, Mapping[str, np.ndarray]],
    trip_dirs_for_opt: Mapping[str, np.ndarray],
    mcs_allowed: Sequence[int],
    ms_allowed: Sequence[int],
    min_cluster_size_for_cohesion: int = 2,
    cohesion_floor: Optional[float] = None,
    cluster_selection_method: str = "eom",
    core_dist_n_jobs: int = 1,
):
    x_scaled = data_for_opt[city_id]["X_scaled"]
    u_trips = np.asarray(trip_dirs_for_opt[city_id])
    if x_scaled.shape[0] != u_trips.shape[0]:
        raise ValueError(f"{city_id}: X_scaled rows != trip_dirs rows")

    cache: Dict[tuple[int, int], np.ndarray] = {}
    mcs_allowed = np.asarray(mcs_allowed, dtype=int)
    ms_allowed = np.asarray(ms_allowed, dtype=int)

    def objfun(x: np.ndarray) -> np.ndarray:
        mcs = snap_scalar(float(x[0]), mcs_allowed)
        ms = snap_scalar(float(x[1]), ms_allowed)
        key = (int(mcs), int(ms))
        if key in cache:
            return cache[key]

        val = eval_hdbscan_params(
            min_cluster_size=mcs,
            min_samples=ms,
            X_scaled=x_scaled,
            u_trips=u_trips,
            min_cluster_size_for_cohesion=min_cluster_size_for_cohesion,
            cohesion_floor=cohesion_floor,
            cluster_selection_method=cluster_selection_method,
            core_dist_n_jobs=core_dist_n_jobs,
        )
        cache[key] = val
        return val

    info = {
        "mcs_allowed": mcs_allowed.tolist(),
        "ms_allowed": ms_allowed.tolist(),
        "mcs_min": int(mcs_allowed.min()),
        "mcs_max": int(mcs_allowed.max()),
        "ms_min": int(ms_allowed.min()),
        "ms_max": int(ms_allowed.max()),
        "N_eff": int(x_scaled.shape[0]),
    }
    return objfun, info


def subsample_for_optimization(
    features: Mapping[str, CityFeatureResult],
    trip_dirs_full: Mapping[str, np.ndarray],
    n_sub: int = 10_000,
    seed: int = 1,
) -> tuple[Dict[str, Dict[str, np.ndarray]], Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    rng = np.random.default_rng(seed)
    data_for_opt: Dict[str, Dict[str, np.ndarray]] = {}
    trip_dirs_for_opt: Dict[str, np.ndarray] = {}
    sub_idx_by_city: Dict[str, np.ndarray] = {}

    for city_id, feat in features.items():
        x = feat.X_scaled
        n = x.shape[0]
        idx = rng.choice(n, size=n_sub, replace=False) if n > n_sub else np.arange(n)
        data_for_opt[city_id] = {"X_scaled": x[idx, :]}
        trip_dirs_for_opt[city_id] = np.asarray(trip_dirs_full[city_id])[idx, :]
        sub_idx_by_city[city_id] = idx

    return data_for_opt, trip_dirs_for_opt, sub_idx_by_city


def run_nsga2_hdbscan(
    features: Mapping[str, CityFeatureResult],
    trip_dirs_full: Mapping[str, np.ndarray],
    n_sub: int = 10_000,
    ms_allowed: Sequence[int] = (10, 15, 20, 25, 30, 35, 40, 50, 60),
    popsize: int = 100,
    generations: int = 50,
    cohesion_floor: Optional[float] = None,
    cluster_selection_method: str = "eom",
    core_dist_n_jobs: int = 1,
) -> Dict[str, object]:
    try:
        from pymoo.algorithms.moo.nsga2 import NSGA2
        from pymoo.core.problem import ElementwiseProblem
        from pymoo.operators.crossover.sbx import SBX
        from pymoo.operators.mutation.pm import PM
        from pymoo.operators.sampling.rnd import FloatRandomSampling
        from pymoo.optimize import minimize
        from pymoo.termination import get_termination
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError("run_nsga2_hdbscan requires pymoo. Install with: pip install pymoo") from exc

    data_for_opt, trip_dirs_for_opt, sub_idx_by_city = subsample_for_optimization(
        features=features,
        trip_dirs_full=trip_dirs_full,
        n_sub=n_sub,
        seed=1,
    )

    bounds_used: Dict[str, Dict[str, object]] = {}
    nsga_results: Dict[str, object] = {}
    pareto_tables_hdb: Dict[str, pd.DataFrame] = {}
    knee_solutions_hdb: Dict[str, pd.DataFrame] = {}

    for city_id in features.keys():
        n_eff = data_for_opt[city_id]["X_scaled"].shape[0]
        mcs_allowed = hdbscan_mcs_candidates(n_eff=n_eff, n=14, frac_min=0.005, frac_max=0.10)

        objfun, info = make_objfun_city_cached(
            city_id=city_id,
            data_for_opt=data_for_opt,
            trip_dirs_for_opt=trip_dirs_for_opt,
            mcs_allowed=mcs_allowed,
            ms_allowed=ms_allowed,
            min_cluster_size_for_cohesion=2,
            cohesion_floor=cohesion_floor,
            cluster_selection_method=cluster_selection_method,
            core_dist_n_jobs=core_dist_n_jobs,
        )
        bounds_used[city_id] = info

        mcs_min, mcs_max = info["mcs_min"], info["mcs_max"]
        ms_min, ms_max = info["ms_min"], info["ms_max"]

        class CityProblem(ElementwiseProblem):
            def __init__(self):
                super().__init__(
                    n_var=2,
                    n_obj=2,
                    n_constr=0,
                    xl=np.array([mcs_min, ms_min], dtype=float),
                    xu=np.array([mcs_max, ms_max], dtype=float),
                )

            def _evaluate(self, x, out, *args, **kwargs):
                out["F"] = objfun(x)

        problem = CityProblem()
        algorithm = NSGA2(
            pop_size=popsize,
            sampling=FloatRandomSampling(),
            crossover=SBX(prob=0.9, eta=15),
            mutation=PM(eta=20),
            eliminate_duplicates=True,
        )

        res = minimize(
            problem,
            algorithm,
            get_termination("n_gen", generations),
            seed=1,
            save_history=False,
            verbose=False,
        )
        nsga_results[city_id] = res

        xvals = np.asarray(res.X)
        fvals = np.asarray(res.F)
        mcs_snap = np.array([snap_scalar(v, info["mcs_allowed"]) for v in xvals[:, 0]], dtype=int)
        ms_snap = np.array([snap_scalar(v, info["ms_allowed"]) for v in xvals[:, 1]], dtype=int)

        df = pd.DataFrame(
            {
                "city": city_id,
                "min_cluster_size": mcs_snap,
                "min_samples": ms_snap,
                "coverage": -fvals[:, 0],
                "dir_cohesion": -fvals[:, 1],
            }
        )
        df = df[(df["coverage"] >= 0) & (df["dir_cohesion"] >= 0)].copy()
        df = df.drop_duplicates(subset=["min_cluster_size", "min_samples"]).reset_index(drop=True)
        pareto_tables_hdb[city_id] = df

        if len(df) > 0:
            cov = df["coverage"].to_numpy(dtype=float)
            coh = df["dir_cohesion"].to_numpy(dtype=float)
            cov_n = (cov - cov.min()) / (cov.max() - cov.min() + 1e-12)
            coh_n = (coh - coh.min()) / (coh.max() - coh.min() + 1e-12)
            knee_idx = int(np.argmax(np.minimum(cov_n, coh_n)))
            knee_solutions_hdb[city_id] = df.iloc[[knee_idx]].copy()
        else:
            knee_solutions_hdb[city_id] = pd.DataFrame()

    return {
        "data_for_opt": data_for_opt,
        "trip_dirs_for_opt": trip_dirs_for_opt,
        "sub_idx_by_city": sub_idx_by_city,
        "bounds_used": bounds_used,
        "nsga_results": nsga_results,
        "pareto_tables_hdb": pareto_tables_hdb,
        "knee_solutions_hdb": knee_solutions_hdb,
    }


def clean_pareto_table_hdb(pareto_tables_hdb: Mapping[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    out: Dict[str, pd.DataFrame] = {}
    for city_id, df in pareto_tables_hdb.items():
        if df is None or len(df) == 0:
            out[city_id] = pd.DataFrame()
            continue

        df2 = df.loc[:, ["city", "min_cluster_size", "min_samples", "coverage", "dir_cohesion"]].copy()
        df2["min_cluster_size"] = df2["min_cluster_size"].astype(int)
        df2["min_samples"] = df2["min_samples"].astype(int)
        df2["coverage"] = pd.to_numeric(df2["coverage"], errors="coerce")
        df2["dir_cohesion"] = pd.to_numeric(df2["dir_cohesion"], errors="coerce")

        df2 = df2[np.isfinite(df2["coverage"]) & np.isfinite(df2["dir_cohesion"])].copy()
        df2 = df2[(df2["coverage"] >= 0) & (df2["dir_cohesion"] >= 0)].copy()
        out[city_id] = df2.drop_duplicates(subset=["min_cluster_size", "min_samples"]).reset_index(drop=True)

    return out


def select_configs_hdb(df: pd.DataFrame, knee_row: Optional[pd.DataFrame], n_each: int = 2) -> pd.DataFrame:
    if df is None or len(df) == 0:
        return pd.DataFrame()

    if knee_row is None or (not isinstance(knee_row, pd.DataFrame)) or len(knee_row) == 0:
        return df.sort_values(["coverage", "dir_cohesion"], ascending=[False, False]).head(min(5, len(df))).copy()

    kmcs = int(knee_row["min_cluster_size"].iloc[0])
    kms = int(knee_row["min_samples"].iloc[0])
    knee_match = df[(df["min_cluster_size"] == kmcs) & (df["min_samples"] == kms)].copy()

    if len(knee_match) == 0:
        kcov = float(knee_row["coverage"].iloc[0])
        kcoh = float(knee_row["dir_cohesion"].iloc[0])
        knee_match = df[(np.abs(df["coverage"] - kcov) < 1e-12) & (np.abs(df["dir_cohesion"] - kcoh) < 1e-12)].copy()

    if len(knee_match) == 0:
        knee_match = df.sort_values(["coverage", "dir_cohesion"], ascending=[False, False]).head(1).copy()
    else:
        knee_match = knee_match.head(1).copy()

    knee_cov = float(knee_match["coverage"].iloc[0])
    knee_coh = float(knee_match["dir_cohesion"].iloc[0])

    higher_cov = df[df["coverage"] > knee_cov].sort_values(["coverage", "dir_cohesion"], ascending=[True, False]).head(n_each)
    higher_coh = df[df["dir_cohesion"] > knee_coh].sort_values(["dir_cohesion", "coverage"], ascending=[True, False]).head(n_each)

    chosen = pd.concat([knee_match, higher_cov, higher_coh], ignore_index=True)
    return chosen.drop_duplicates(subset=["min_cluster_size", "min_samples"]).head(1 + 2 * n_each).reset_index(drop=True)


def full_eval_one_hdb(
    city_id: str,
    min_cluster_size: int,
    min_samples: int,
    features: Mapping[str, CityFeatureResult],
    trip_dirs_full: Mapping[str, np.ndarray],
    cluster_selection_method: str = "eom",
    core_dist_n_jobs: int = 1,
    min_cluster_size_for_cohesion: int = 2,
) -> pd.DataFrame:
    x_full = features[city_id].X_scaled
    u_full = np.asarray(trip_dirs_full[city_id], dtype=float)
    if u_full.shape[0] != x_full.shape[0]:
        raise ValueError(f"{city_id}: trip_dirs_full rows != X_scaled rows")

    clusterer = _require_package("hdbscan", "HDBSCAN clustering").HDBSCAN(
        min_cluster_size=int(min_cluster_size),
        min_samples=int(min_samples),
        metric="euclidean",
        cluster_selection_method=cluster_selection_method,
        core_dist_n_jobs=int(core_dist_n_jobs),
    )
    labels = clusterer.fit_predict(x_full).astype(int)
    labels_dbscan_style = np.where(labels == -1, 0, labels + 1)

    res = eval_solution(
        labels=labels_dbscan_style,
        u_trips=u_full,
        noise_label=0,
        min_cluster_size=int(min_cluster_size_for_cohesion),
    )

    return pd.DataFrame(
        [
            {
                "city": city_id,
                "min_cluster_size": int(min_cluster_size),
                "min_samples": int(min_samples),
                "coverage": float(res["coverage"]),
                "clustered_pct": 100.0 * float(res["coverage"]),
                "dir_cohesion": float(res["dir_cohesion"]),
                "n_clusters": int(res["K"]),
                "noise_pct": 100.0 * float((labels_dbscan_style == 0).mean()),
                "n_total": int(res["N_total"]),
            }
        ]
    )


def validate_selected_hdbscan_configs(
    features: Mapping[str, CityFeatureResult],
    trip_dirs_full: Mapping[str, np.ndarray],
    pareto_tables_hdb: Mapping[str, pd.DataFrame],
    knee_solutions_hdb: Mapping[str, pd.DataFrame],
    hdbscan_params: Optional[Mapping[str, Mapping[str, int]]] = None,
    min_cluster_size_fixed: int = 20,
) -> Dict[str, object]:
    pareto_unique_hdb = clean_pareto_table_hdb(pareto_tables_hdb)
    full_validated_hdb: Dict[str, pd.DataFrame] = {}

    for city_id, df in pareto_unique_hdb.items():
        if df is None or len(df) == 0:
            continue
        knee = knee_solutions_hdb.get(city_id, pd.DataFrame())
        chosen = select_configs_hdb(df, knee_row=knee, n_each=2)
        pareto_out = pd.concat(
            [
                full_eval_one_hdb(
                    city_id=city_id,
                    min_cluster_size=int(chosen.loc[i, "min_cluster_size"]),
                    min_samples=int(chosen.loc[i, "min_samples"]),
                    features=features,
                    trip_dirs_full=trip_dirs_full,
                )
                for i in range(len(chosen))
            ],
            ignore_index=True,
        )
        pareto_out["config"] = "pareto_local"

        expert_out = None
        if hdbscan_params is not None and city_id in hdbscan_params:
            expert_out = full_eval_one_hdb(
                city_id=city_id,
                min_cluster_size=int(min_cluster_size_fixed),
                min_samples=int(hdbscan_params[city_id]["minPts"]),
                features=features,
                trip_dirs_full=trip_dirs_full,
            )
            expert_out["config"] = "expert"

        full_validated_hdb[city_id] = pd.concat([pareto_out, expert_out], ignore_index=True) if expert_out is not None else pareto_out

    full_validated_hdb_df = pd.concat(full_validated_hdb.values(), ignore_index=True) if len(full_validated_hdb) > 0 else pd.DataFrame()
    knee_full_hdb = pd.DataFrame()

    if len(full_validated_hdb_df) > 0:
        dfp = full_validated_hdb_df[full_validated_hdb_df["config"] == "pareto_local"].copy()
        rows: List[pd.DataFrame] = []

        for city_id in dfp["city"].unique():
            d = dfp[dfp["city"] == city_id].copy()
            cov = d["coverage"].to_numpy(dtype=float)
            coh = d["dir_cohesion"].to_numpy(dtype=float)
            cov_n = (cov - cov.min()) / (cov.max() - cov.min() + 1e-12)
            coh_n = (coh - coh.min()) / (coh.max() - coh.min() + 1e-12)
            dist_ideal = np.sqrt((1 - cov_n) ** 2 + (1 - coh_n) ** 2)
            j = int(np.argmin(dist_ideal))
            r = d.iloc[[j]].copy()
            r["config"] = "knee_full"
            rows.append(r)

        knee_full_hdb = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()

    compare_df = pd.DataFrame()
    if len(full_validated_hdb_df) > 0 and len(knee_full_hdb) > 0:
        compare_df = full_validated_hdb_df[full_validated_hdb_df["config"] == "expert"].copy()
        compare_df = pd.concat([compare_df, knee_full_hdb], ignore_index=True)
        cols = ["city", "config", "min_cluster_size", "min_samples", "coverage", "clustered_pct", "dir_cohesion", "n_clusters", "noise_pct"]
        compare_df = compare_df.loc[:, [c for c in cols if c in compare_df.columns]]

    return {
        "pareto_unique_hdb": pareto_unique_hdb,
        "full_validated_hdb": full_validated_hdb,
        "full_validated_hdb_df": full_validated_hdb_df,
        "knee_full_hdb": knee_full_hdb,
        "expert_vs_knee_full": compare_df,
    }


def nondominated_min(f: np.ndarray) -> np.ndarray:
    f = np.asarray(f, dtype=float)
    n = f.shape[0]
    if n == 0:
        return f
    keep = np.ones(n, dtype=bool)
    for i in range(n):
        if not keep[i]:
            continue
        for j in range(n):
            if i == j or not keep[i]:
                continue
            if np.all(f[j, :] <= f[i, :]) and np.any(f[j, :] < f[i, :]):
                keep[i] = False
    return f[keep, :]


def hv2d_min(f: np.ndarray, ref: np.ndarray) -> float:
    f = np.asarray(f, dtype=float)
    ref = np.asarray(ref, dtype=float)
    if f.shape[0] == 0:
        return 0.0

    f = nondominated_min(f)
    order = np.lexsort((f[:, 1], f[:, 0]))
    f = f[order, :]

    hv = 0.0
    h = float(ref[1])
    for i in range(f.shape[0]):
        f1 = float(f[i, 0])
        f2 = float(f[i, 1])
        if f2 < h:
            hv += (float(ref[0]) - f1) * (h - f2)
            h = f2
    return float(hv)


def hv2d_contrib_min(f: np.ndarray, ref: np.ndarray) -> np.ndarray:
    hv_all = hv2d_min(f, ref)
    n = f.shape[0]
    if n == 0:
        return np.array([], dtype=float)
    return np.array([hv_all - hv2d_min(np.delete(f, i, axis=0), ref) for i in range(n)], dtype=float)


def safe_norm(z: np.ndarray, zmin: float, zmax: float) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    if (not np.isfinite(zmin)) or (not np.isfinite(zmax)) or np.isclose(zmax, zmin):
        return np.zeros_like(z)
    return (z - zmin) / (zmax - zmin)


def gini(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    x = np.sort(np.maximum(0.0, x))
    s = float(x.sum())
    if s <= 0:
        return 0.0
    n = x.size
    return float((2.0 * (np.arange(1, n + 1) * x).sum() / (n * s)) - (n + 1) / n)


def hypervolume_analysis(
    pareto_by_algo: Mapping[str, Mapping[str, pd.DataFrame]],
    ref_min: np.ndarray = np.array([1.05, 1.05], dtype=float),
) -> Dict[str, pd.DataFrame]:
    all_rows = []
    for algo, pt_list in pareto_by_algo.items():
        if pt_list is None:
            continue
        for city_id, df in pt_list.items():
            if df is None or len(df) == 0:
                continue
            if not all(c in df.columns for c in ["coverage", "dir_cohesion"]):
                raise ValueError(f"{algo}/{city_id}: missing required columns")
            df2 = df.loc[:, ["coverage", "dir_cohesion"]].copy()
            df2["algorithm"] = algo
            df2["city"] = city_id
            df2 = df2[np.isfinite(df2["coverage"]) & np.isfinite(df2["dir_cohesion"])].copy()
            if len(df2) > 0:
                all_rows.append(df2)

    all_pts = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    if len(all_pts) == 0:
        return {"hv_results": pd.DataFrame(), "hv_summary": pd.DataFrame()}

    cov_min = float(all_pts["coverage"].min())
    cov_max = float(all_pts["coverage"].max())
    coh_min = float(all_pts["dir_cohesion"].min())
    coh_max = float(all_pts["dir_cohesion"].max())

    hv_max = float(np.prod(ref_min - np.array([0.0, 0.0])))
    out_rows = []

    for algo, pt_list in pareto_by_algo.items():
        for city_id, df in pt_list.items():
            if df is None or len(df) == 0:
                out_rows.append({"algorithm": algo, "city": city_id, "n_points": 0, "hv_raw": np.nan, "hv_unit": np.nan, "contrib_gini": np.nan, "contrib_top10_share": np.nan, "contrib_cv": np.nan})
                continue

            df2 = df.loc[:, ["coverage", "dir_cohesion"]].copy()
            df2 = df2[np.isfinite(df2["coverage"]) & np.isfinite(df2["dir_cohesion"])].copy()
            if len(df2) == 0:
                out_rows.append({"algorithm": algo, "city": city_id, "n_points": 0, "hv_raw": np.nan, "hv_unit": np.nan, "contrib_gini": np.nan, "contrib_top10_share": np.nan, "contrib_cv": np.nan})
                continue

            cov_n = safe_norm(df2["coverage"].to_numpy(), cov_min, cov_max)
            coh_n = safe_norm(df2["dir_cohesion"].to_numpy(), coh_min, coh_max)
            fmin = np.column_stack([1.0 - cov_n, 1.0 - coh_n])
            hv_raw = hv2d_min(fmin, ref_min)
            hv_unit = hv_raw / hv_max if hv_max > 0 else np.nan

            contrib = hv2d_contrib_min(fmin, ref_min)
            contrib = contrib[np.isfinite(contrib)]
            if contrib.size == 0 or np.isclose(contrib.sum(), 0.0):
                top10_share = np.nan
                contrib_cv = np.nan
            else:
                n_top = max(1, int(np.floor(0.1 * contrib.size)))
                top10_share = float(np.sort(contrib)[::-1][:n_top].sum() / contrib.sum())
                contrib_cv = float(np.std(contrib) / np.mean(contrib)) if (contrib.size > 1 and np.mean(contrib) > 0) else np.nan

            out_rows.append(
                {
                    "algorithm": algo,
                    "city": city_id,
                    "n_points": int(len(df2)),
                    "hv_raw": float(hv_raw),
                    "hv_unit": float(hv_unit),
                    "contrib_gini": float(gini(contrib)),
                    "contrib_top10_share": float(top10_share) if np.isfinite(top10_share) else np.nan,
                    "contrib_cv": float(contrib_cv) if np.isfinite(contrib_cv) else np.nan,
                }
            )

    hv_results = pd.DataFrame(out_rows).sort_values(["algorithm", "hv_unit", "city"], ascending=[True, False, True]).reset_index(drop=True)
    hv_summary = hv_results.loc[:, ["algorithm", "city", "hv_unit", "contrib_gini"]].rename(columns={"contrib_gini": "gini"})
    return {"hv_results": hv_results, "hv_summary": hv_summary}


def get_hdbscan_selected_params(
    city_ids: Sequence[str],
    knee_full_hdb: Optional[pd.DataFrame] = None,
    knee_solutions_hdb: Optional[Mapping[str, pd.DataFrame]] = None,
    hdbscan_params: Optional[Mapping[str, Mapping[str, int]]] = None,
    min_cluster_size_fixed: int = 20,
) -> Dict[str, Dict[str, int]]:
    out: Dict[str, Dict[str, int]] = {}
    for city_id in city_ids:
        if isinstance(knee_full_hdb, pd.DataFrame) and len(knee_full_hdb) > 0:
            row = knee_full_hdb[knee_full_hdb["city"] == city_id].head(1)
            if len(row) == 1:
                out[city_id] = {
                    "min_cluster_size": int(row["min_cluster_size"].iloc[0]),
                    "min_samples": int(row["min_samples"].iloc[0]),
                }
                continue

        if knee_solutions_hdb is not None:
            knee = knee_solutions_hdb.get(city_id, None)
            if isinstance(knee, pd.DataFrame) and len(knee) > 0:
                out[city_id] = {
                    "min_cluster_size": int(knee["min_cluster_size"].iloc[0]),
                    "min_samples": int(knee["min_samples"].iloc[0]),
                }
                continue

        if hdbscan_params is not None and city_id in hdbscan_params:
            out[city_id] = {
                "min_cluster_size": int(min_cluster_size_fixed),
                "min_samples": int(hdbscan_params[city_id]["minPts"]),
            }
            continue

        raise ValueError(f"No selected parameters found for city: {city_id}")

    return out


def plot_flow_topn_clean(
    city_df: pd.DataFrame,
    labels_dbscan_style: np.ndarray,
    city_id: str,
    algo_name: str = "HDBSCAN",
    noise_label: int = 0,
    out_dir: str | PathLike[str] = "plots_flows_top3",
    top_n_clusters: int = 3,
    sample_n_segments_top: int = 25_000,
    seed: int = 1,
    show_plot: bool = True,
 ) -> Optional[Dict[str, object]]:
    plt = _require_package("matplotlib.pyplot", "plotting flow maps")
    line_collection = _require_package("matplotlib.collections", "plotting flow maps")
    LineCollection = line_collection.LineCollection

    os.makedirs(out_dir, exist_ok=True)
    req = ["x_o", "y_o", "x_d", "y_d"]
    missing = [c for c in req if c not in city_df.columns]
    if missing:
        raise ValueError(f"{city_id}: Missing columns in city_clean: {missing}")

    labels = np.asarray(labels_dbscan_style, dtype=int)
    if labels.shape[0] != city_df.shape[0]:
        raise ValueError(f"{city_id}: labels length != nrow(city_clean)")

    cl_nonnoise = labels[(labels != noise_label) & np.isfinite(labels)]
    if cl_nonnoise.size == 0:
        return None

    uniq, cnt = np.unique(cl_nonnoise, return_counts=True)
    top_ids = uniq[np.argsort(-cnt)][:top_n_clusters].astype(int)

    mask_top = np.isin(labels, top_ids)
    df_top = city_df.loc[mask_top, ["x_o", "y_o", "x_d", "y_d"]].copy()
    cl_top = labels[mask_top]

    rng = np.random.default_rng(seed)
    if len(df_top) > sample_n_segments_top:
        idx = rng.choice(len(df_top), size=sample_n_segments_top, replace=False)
        df_top = df_top.iloc[idx, :]
        cl_top = cl_top[idx]

    segs = np.stack([df_top[["x_o", "y_o"]].to_numpy(dtype=float), df_top[["x_d", "y_d"]].to_numpy(dtype=float)], axis=1)
    top_ids_list = top_ids.tolist()
    cl_idx = np.array([top_ids_list.index(int(c)) for c in cl_top], dtype=int)

    fig, ax = plt.subplots(figsize=(9, 7.5))
    lc = LineCollection(segs, linewidths=0.30, alpha=0.28)
    lc.set_array(cl_idx)
    ax.add_collection(lc)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(f"{algo_name} — Top {top_n_clusters} clusters — {city_id.upper()}", fontweight="bold")
    ax.set_xlabel("UTM Easting (m)")
    ax.set_ylabel("UTM Northing (m)")
    ax.autoscale()

    cbar = fig.colorbar(lc, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Cluster (top-N)")
    cbar.set_ticks(np.arange(len(top_ids)))
    cbar.set_ticklabels([str(int(x)) for x in top_ids])
    ax.grid(True, alpha=0.25)

    f_out = os.path.join(str(out_dir), f"flows_{algo_name}_{city_id}_top{top_n_clusters}.png")
    plt.savefig(f_out, dpi=300, bbox_inches="tight")
    if show_plot:
        plt.show()
    plt.close(fig)

    return {"file": f_out, "top_cluster_ids": top_ids.tolist()}



def install_required_packages_colab() -> None:
    """Install required packages in Colab/Jupyter runtime."""
    import subprocess
    import sys

    commands = [
        [sys.executable, "-m", "pip", "install", "-q", "pyreadr", "pyproj", "pandas", "numpy"],
        [sys.executable, "-m", "pip", "install", "-q", "scikit-learn", "hdbscan"],
        [sys.executable, "-m", "pip", "install", "-q", "matplotlib"],
        [sys.executable, "-m", "pip", "install", "-q", "pymoo"],
    ]

    for cmd in commands:
        subprocess.check_call(cmd)

    print("All required packages installed successfully.")



def run_pipeline(
    city_files: Mapping[str, str],
    k_values_dbscan: Sequence[int] = (10, 20, 30, 50),
    k_allowed_snn: Sequence[int] = (15, 20, 25, 30, 35),
    min_samples_grid: Sequence[int] = (20, 25, 30, 35, 40, 45, 50, 60, 70),
    min_cluster_size_fixed: int = 20,
    n_sub: int = 10_000,
    final_hdbscan_params: Optional[Mapping[str, Mapping[str, int]]] = None,
    run_nsga2: bool = False,
) -> Dict[str, object]:
    if final_hdbscan_params is None:
        final_hdbscan_params = {
            "berlin": {"minPts": 50},
            "munich": {"minPts": 45},
            "cologne": {"minPts": 50},
        }

    prepared = {
        city_id: prep_city_od(city_id=city_id, rds_path=Path(path))
        for city_id, path in city_files.items()
    }
    summary_df = cleaning_summary(prepared)

    features = {
        city_id: build_features(city_id=city_id, city_df=prep.data)
        for city_id, prep in prepared.items()
    }

    knn_rows: List[Dict[str, float | int | str]] = []
    for city_id, feat in features.items():
        for k in k_values_dbscan:
            row = knn_summary_one(feat.X_scaled, int(k))
            row["city"] = city_id
            knn_rows.append(row)
    knn_df = pd.DataFrame(knn_rows)

    snn_rows: List[Dict[str, float | int | str]] = []
    for city_id, feat in features.items():
        for k in k_allowed_snn:
            row = snn_sharednn_diagnostics(feat.X_scaled, k=int(k), n_probe=4000, seed=1)
            row["city"] = city_id
            snn_rows.append(row)
    snn_df = pd.DataFrame(snn_rows)

    mcs_vals = hdbscan_mcs_candidates(n_eff=n_sub)
    hdb_bounds_df = pd.DataFrame(
        [
            {
                "N_eff": n_sub,
                "mcs_min": int(mcs_vals.min()),
                "mcs_max": int(mcs_vals.max()),
                "mcs_values": ", ".join(map(str, mcs_vals.tolist())),
            }
        ]
    )

    grid_df = evaluate_hdbscan_grid(
        features=features,
        min_samples_grid=min_samples_grid,
        min_cluster_size=min_cluster_size_fixed,
        cluster_selection_method="eom",
    )

    final_city_tables = fit_final_hdbscan(
        features=features,
        per_city_params=final_hdbscan_params,
        min_cluster_size=min_cluster_size_fixed,
        cluster_selection_method="eom",
    )

    clustering_checks = verify_final_clustering(final_city_tables)
    cluster_summaries = summarize_clusters(final_city_tables)
    moo_clusters = analyze_dominant_cluster(cluster_summaries=cluster_summaries, features=features)
    trip_dirs_full = build_trip_direction_vectors(final_city_tables)

    sanity_eval_rows: List[Dict[str, float | int | str]] = []
    rng = np.random.default_rng(1)
    for city_id, u_trips in trip_dirs_full.items():
        n = int(u_trips.shape[0])
        labels_all_noise = np.zeros(n, dtype=int)
        labels_one_cluster = np.ones(n, dtype=int)
        labels_toy = np.zeros(n, dtype=int)
        clustered_idx = rng.choice(n, size=int(np.floor(0.7 * n)), replace=False)
        labels_toy[clustered_idx] = rng.integers(1, 21, size=clustered_idx.size)

        for name, labels in [
            ("all_noise", labels_all_noise),
            ("one_cluster", labels_one_cluster),
            ("toy_70pct", labels_toy),
        ]:
            row = eval_solution(labels=labels, u_trips=u_trips, noise_label=0)
            row["city"] = city_id
            row["scenario"] = name
            sanity_eval_rows.append(row)

    sanity_eval_df = pd.DataFrame(sanity_eval_rows)

    outputs: Dict[str, object] = {
        "prepared": prepared,
        "summary": summary_df,
        "features": features,
        "knn_df": knn_df,
        "snn_df": snn_df,
        "hdb_bounds_df": hdb_bounds_df,
        "grid_df_hdbscan": grid_df,
        "final_city_tables": final_city_tables,
        "clustering_checks": clustering_checks,
        "cluster_summaries": cluster_summaries,
        "moo_clusters": moo_clusters,
        "trip_dirs_full": trip_dirs_full,
        "sanity_eval_df": sanity_eval_df,
    }

    if run_nsga2:
        nsga2_out = run_nsga2_hdbscan(features=features, trip_dirs_full=trip_dirs_full, n_sub=n_sub)
        outputs["nsga2"] = nsga2_out

        step15 = validate_selected_hdbscan_configs(
            features=features,
            trip_dirs_full=trip_dirs_full,
            pareto_tables_hdb=nsga2_out["pareto_tables_hdb"],
            knee_solutions_hdb=nsga2_out["knee_solutions_hdb"],
            hdbscan_params=final_hdbscan_params,
            min_cluster_size_fixed=min_cluster_size_fixed,
        )
        outputs["step15_validation"] = step15

        hv_out = hypervolume_analysis({"HDBSCAN": nsga2_out["pareto_tables_hdb"]})
        outputs["step16_hv"] = hv_out

    return outputs


if __name__ == "__main__":
    CITY_FILES = {
        "berlin": "dt_bolt_berlin_06_05.rds",
        "cologne": "dt_voi_cologne_06_05.rds",
        "munich": "dt_voi_munich_06_05.rds",
    }

    print("\n=== Colab install helper ===")
    print("Call install_required_packages_colab() once in Colab before running pipeline.")

    plt = _require_package("matplotlib.pyplot", "plotting diagnostics in __main__")

    outputs = run_pipeline(city_files=CITY_FILES, run_nsga2=True)

    print("\n=== Cleaning summary ===")
    print(outputs["summary"].to_string(index=False))

    print("\n=== kNN diagnostics ===")
    print(outputs["knn_df"].to_string(index=False))

    print("\n=== SNN diagnostics ===")
    print(outputs["snn_df"].to_string(index=False))

    print("\n=== HDBSCAN bounds ===")
    print(outputs["hdb_bounds_df"].to_string(index=False))

    print("\n=== HDBSCAN grid (full) ===")
    print(outputs["grid_df_hdbscan"].to_string(index=False))

    print("\n=== Final clustering checks ===")
    print(pd.DataFrame(outputs["clustering_checks"]).T.to_string())

    print("\n=== Objective sanity checks ===")
    print(outputs["sanity_eval_df"].to_string(index=False))

    for city_id, tbl in outputs["final_city_tables"].items():
        n_clusters = int(tbl.loc[tbl["cluster"] != 0, "cluster"].nunique())
        noise_pct = 100.0 * float((tbl["cluster"] == 0).mean())
        print(f"FINAL | {city_id} | clusters={n_clusters} | noise%={noise_pct:.2f}")

    if "nsga2" in outputs:
        nsga2_out = outputs["nsga2"]
        print("\n=== Pareto tables (HDBSCAN) ===")
        for city_id, df in nsga2_out["pareto_tables_hdb"].items():
            print(f"\nCITY: {city_id.upper()} | N Pareto points: {len(df)}")
            print(df.sort_values(["coverage", "dir_cohesion"]).to_string(index=False) if len(df) else "No Pareto rows")

            if len(df):
                plt.figure()
                plt.scatter(df["coverage"], df["dir_cohesion"], s=20, alpha=0.8)
                knee = nsga2_out["knee_solutions_hdb"].get(city_id, pd.DataFrame())
                if isinstance(knee, pd.DataFrame) and len(knee):
                    plt.scatter(knee["coverage"], knee["dir_cohesion"], s=120, marker="x")
                plt.title(f"HDBSCAN Pareto — {city_id.upper()}")
                plt.xlabel("Coverage")
                plt.ylabel("Directional cohesion")
                plt.grid(True)
                plt.show()

    if "step15_validation" in outputs:
        step15 = outputs["step15_validation"]
        print("\n=== FULL validation of selected Pareto configs ===")
        print(step15["full_validated_hdb_df"].to_string(index=False) if len(step15["full_validated_hdb_df"]) else "No full validation rows")
        print("\n=== EXPERT vs KNEE_FULL ===")
        print(step15["expert_vs_knee_full"].to_string(index=False) if len(step15["expert_vs_knee_full"]) else "No comparison rows")

    if "step16_hv" in outputs:
        hv = outputs["step16_hv"]
        print("\n=== Hypervolume results ===")
        print(hv["hv_results"].to_string(index=False) if len(hv["hv_results"]) else "No HV rows")

        if len(hv["hv_summary"]):
            hv_summary = hv["hv_summary"]
            cities = sorted(hv_summary["city"].unique().tolist())
            algos = sorted(hv_summary["algorithm"].unique().tolist())
            x = np.arange(len(cities))
            width = 0.8 / max(1, len(algos))

            plt.figure()
            for i, algo in enumerate(algos):
                y = []
                for c in cities:
                    r = hv_summary[(hv_summary["city"] == c) & (hv_summary["algorithm"] == algo)]
                    y.append(float(r["hv_unit"].iloc[0]) if len(r) else np.nan)
                plt.bar(x + i * width, y, width=width, label=algo)
            plt.xticks(x + width * (len(algos) - 1) / 2, [c.upper() for c in cities])
            plt.ylabel("Hypervolume (unit-scaled)")
            plt.title("Unit-scaled hypervolume by city and algorithm")
            plt.legend()
            plt.grid(axis="y", alpha=0.3)
            plt.show()
