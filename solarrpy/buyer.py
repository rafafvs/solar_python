"""Phase 5: Buyer/seller cash-flow analysis for the SoRadIDX hedge.

Computes per-year optimal hedge ratio N_c and MAPE metrics, matching
the paper's §4.3 perspective.

Model parameters (paper notation):
    E_f   = 0.15    fraction of capacity factor
    N_p   = 10 000  number of panels
    E_bar = 0.08    energy output per panel per unit GHI  [kWh / (panel·kWh/m²)]
    tick  = 0.10    contract tick size  [€ / kWh/m²]

Cash flows (daily, day n):
    unhedged  : CF_u(R_n)  = E_f · N_p · R_n · E_bar
    benchmark : CF_b(K_n)  = E_f · N_p · K_n · E_bar   (at seasonal strike)
    payoff    : P_n        = max(K_n − R_n, 0) · tick
    V0 cost   : C_n        = V0_n · tick                (premium paid at inception)
    hedged    : CF_h(R_n)  = CF_u(R_n) + N_c · (P_n − C_n)

Optimal N_c minimises ∑(CF_h − CF_b)² — closed-form quadratic solution,
clamped to [0, ∞) (no short selling of protection).

MAPE:
    unhedged MAPE = 100/N · ∑ |CF_u(R_n) − CF_b(K_n)| / CF_b(K_n)
    hedged   MAPE = 100/N · ∑ |CF_h(R_n) − CF_b(K_n)| / CF_b(K_n)

Seller annual net return (in € units):
    seller_net = N_c · tick · ∑_n (V0_n − P_n/tick)
               = N_c · tick · ∑_n (V0_n − max(K_n − R_n, 0))
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from typing import Sequence

from .sorad import sorad_price, daily_strike

# ── Contract / buyer parameters ────────────────────────────────────────────

E_f   = 0.15       # capacity-factor fraction
N_p   = 10_000     # number of panels
E_bar = 0.08       # output per panel per unit GHI
tick  = 0.10       # €/kWh/m²  — contract tick


# ── Cash-flow helpers ─────────────────────────────────────────────────────

def unhedged_cf(R: float | np.ndarray) -> float | np.ndarray:
    """CF_u(R) = E_f · N_p · R · E_bar."""
    return E_f * N_p * np.asarray(R, dtype=float) * E_bar


def benchmark_cf(K: float | np.ndarray) -> float | np.ndarray:
    """CF_b(K) = E_f · N_p · K · E_bar  (unhedged evaluated at seasonal strike)."""
    return E_f * N_p * np.asarray(K, dtype=float) * E_bar


def hedged_cf(
    R: np.ndarray,
    K: np.ndarray,
    V0: np.ndarray,
    N_c: float,
) -> np.ndarray:
    """CF_h = CF_u(R) + N_c · (P_n − V0_n·tick).

    P_n = max(K−R, 0)·tick   (payoff in €)
    V0_n·tick                 (premium cost in €)
    """
    payoff = np.maximum(K - R, 0.0) * tick
    cost   = V0 * tick
    return unhedged_cf(R) + N_c * (payoff - cost)


# ── Optimal N_c (closed form) ─────────────────────────────────────────────

def optimal_Nc(
    R: np.ndarray,
    K: np.ndarray,
    V0: np.ndarray,
) -> float:
    """Minimise ∑(CF_h − CF_b)² over N_c ≥ 0.

    Closed-form solution of the quadratic:
        N_c* = −∑(d_n · h_n) / ∑(h_n²),  clamped to 0

    where d_n = CF_u(R_n) − CF_b(K_n)  and  h_n = P_n − C_n.
    """
    d = unhedged_cf(R) - benchmark_cf(K)    # deviation from benchmark
    h = np.maximum(K - R, 0.0) * tick - V0 * tick   # hedge net P&L per contract
    denom = float(np.dot(h, h))
    if denom < 1e-30:
        return 0.0
    return float(max(0.0, -np.dot(d, h) / denom))


# ── MAPE ──────────────────────────────────────────────────────────────────

def mape(cf_actual: np.ndarray, cf_benchmark: np.ndarray) -> float:
    """Mean Absolute Percentage Error = 100/N · ∑ |cf − cf_b| / |cf_b|."""
    mask = np.abs(cf_benchmark) > 1e-12
    return float(100.0 / mask.sum() * np.sum(np.abs(cf_actual[mask] - cf_benchmark[mask])
                                              / np.abs(cf_benchmark[mask])))


# ── Per-year computation ──────────────────────────────────────────────────

def _buyer_year(
    eval_year: int,
    cal,
    df_full: pd.DataFrame,
    date_col: str = 'date',
    ghi_col: str = 'GHI',
) -> dict:
    """Compute buyer/seller metrics for one evaluation year.

    V0 pricing: each daily contract priced from t = last day of prior year
    (static, matching the V0 mode in soradidx.py).
    """
    df_year = df_full[df_full[date_col].dt.year == eval_year].sort_values(date_col)

    # Conditioning date: last observed day of prior year
    prev_rows = df_full[df_full[date_col].dt.year == (eval_year - 1)].sort_values(date_col)
    if len(prev_rows) == 0:
        t0_row  = df_year.iloc[0]
        R_t0    = float(t0_row[ghi_col])
        t0      = t0_row[date_col] - pd.Timedelta(days=1)
    else:
        t0_row  = prev_rows.iloc[-1]
        R_t0    = float(t0_row[ghi_col])
        t0      = t0_row[date_col]

    t0_str = str(t0.date())

    K_arr  = np.empty(len(df_year))
    V0_arr = np.empty(len(df_year))
    R_arr  = np.empty(len(df_year))

    for i, (_, row) in enumerate(df_year.iterrows()):
        t_hor = row[date_col]
        K_arr[i]  = daily_strike(t_hor, cal)
        V0_arr[i] = sorad_price(K_arr[i], R_t0, t0_str, str(t_hor.date()), cal)
        R_arr[i]  = float(row[ghi_col])

    N_c   = optimal_Nc(R_arr, K_arr, V0_arr)
    cf_u  = unhedged_cf(R_arr)
    cf_b  = benchmark_cf(K_arr)
    cf_h  = hedged_cf(R_arr, K_arr, V0_arr, N_c)

    mape_u = mape(cf_u, cf_b)
    mape_h = mape(cf_h, cf_b)

    # Seller net return (€): collects premium, pays payoffs
    payoff_arr       = np.maximum(K_arr - R_arr, 0.0)
    seller_net_eur   = N_c * tick * float(np.sum(V0_arr - payoff_arr))

    return {
        'eval_year':       eval_year,
        'train_end_year':  eval_year - 1,
        'N_c_opt':         N_c,
        'unhedged_MAPE':   mape_u,
        'hedged_MAPE':     mape_h,
        'seller_net_eur':  seller_net_eur,
        'V0_total':        float(V0_arr.sum()),
        'payoff_total':    float(payoff_arr.sum()),
        'n_days':          len(df_year),
    }


# ── Train-test loop ───────────────────────────────────────────────────────

def compute_buyer_table(
    df_full: pd.DataFrame,
    eval_years: Sequence[int] = (2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023),
    coords: dict | None = None,
    target_col: str = 'GHI',
    clearsky_col: str = 'clearsky',
    date_col: str = 'date',
    progress: bool = True,
) -> pd.DataFrame:
    """Run Phase 5 train-test loop, producing Table 2.

    For each eval year Y: train on 2005..(Y-1), compute buyer metrics for year Y.

    Returns
    -------
    pd.DataFrame with columns:
        eval_year, train_end_year, N_c_opt,
        unhedged_MAPE, hedged_MAPE, seller_net_eur,
        V0_total, payoff_total, n_days
    """
    from .calibration import calibrate_window

    if coords is None:
        coords = {'lat': 44.5}

    records = []
    for eval_year in eval_years:
        train_end = eval_year - 1
        if progress:
            print(f"  [{eval_year}] Training 2005–{train_end}...", end=' ', flush=True)

        df_train = df_full[df_full[date_col].dt.year <= train_end].copy()
        cal = calibrate_window(
            df_train,
            target_col=target_col,
            clearsky_col=clearsky_col,
            date_col=date_col,
            coords=coords,
            train_end_year=train_end,
        )

        if progress:
            print(f"buyer {eval_year}...", end=' ', flush=True)

        rec = _buyer_year(eval_year, cal, df_full, date_col=date_col, ghi_col=target_col)
        records.append(rec)

        if progress:
            print(
                f"N_c={rec['N_c_opt']:.1f}  "
                f"MAPE_u={rec['unhedged_MAPE']:.2f}%  "
                f"MAPE_h={rec['hedged_MAPE']:.2f}%  "
                f"seller={rec['seller_net_eur']:.2f}€"
            )

    return pd.DataFrame(records, columns=[
        'eval_year', 'train_end_year', 'N_c_opt',
        'unhedged_MAPE', 'hedged_MAPE', 'seller_net_eur',
        'V0_total', 'payoff_total', 'n_days',
    ])
