"""Phase 4: SoRadIDX — annual solar radiation index pricing.

Computes the SoRadIDX price as the sum of daily SoRad puts over an evaluation year:

    SoRadIDX_year = sum_{n in year} sorad_price(K_n, ...)

Two pricing modes:
  V0  — static: all contracts priced at t = Jan 1 of eval year (or last day of prior year)
        p_T uses the unconditional monthly prior from cal.monthly_mixture
  Vt  — dynamic: each contract priced at t = n-1 (previous day)
        p_T uses the unconditional monthly prior (λ^R = 0, so Q = P)

Under λ^R = 0: V^P = V^Q*, so V0_P = V0_Q and Vt_P = Vt_Q (same object).

Output: results/Table1_SoRadIDX.csv with columns:
    eval_year, train_end_year,
    V0, Vt,
    V_hist (sum of realized payoffs max(K_n - R_n, 0) over the year)
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from typing import Sequence

from .sorad import sorad_price, daily_strike
from .forecastDensity import clearsky_at


# ---------------------------------------------------------------------------
# Single-year pricing
# ---------------------------------------------------------------------------

def _soradidx_year(
    eval_year: int,
    cal,
    df_full: pd.DataFrame,
    mode: str = 'both',
    date_col: str = 'date',
    ghi_col: str = 'GHI',
) -> dict:
    """Price the SoRadIDX for eval_year using a given CalibrationWindowResult.

    Parameters
    ----------
    eval_year : int  — year whose daily contracts are being priced
    cal       : CalibrationWindowResult  — trained on data up to eval_year - 1
    df_full   : pd.DataFrame  — full dataset (needs prior-year end + eval-year rows)
    mode      : 'V0', 'Vt', 'both', or 'gamma' (realised payoff only — skips pricing)
    date_col  : column name for dates
    ghi_col   : column name for observed GHI

    Returns
    -------
    dict with keys: eval_year, V0, Vt, Gamma, n_days
        Gamma = sum_n tick · max(K_n − R_n, 0) — single-year realised payoff.
    """
    df_year = df_full[df_full[date_col].dt.year == eval_year].copy()
    df_year = df_year.sort_values(date_col).reset_index(drop=True)

    if len(df_year) == 0:
        raise ValueError(f"No data for eval_year={eval_year}")

    # Conditioning date for V0: last observed day before eval_year
    prev_rows = df_full[df_full[date_col].dt.year == (eval_year - 1)]
    if len(prev_rows) == 0:
        t0_row = df_year.iloc[0]
        R_t0   = float(t0_row[ghi_col])
        t0     = t0_row[date_col] - pd.Timedelta(days=1)
    else:
        t0_row = prev_rows.sort_values(date_col).iloc[-1]
        R_t0   = float(t0_row[ghi_col])
        t0     = t0_row[date_col]

    V0    = 0.0
    Vt    = 0.0
    Gamma = 0.0
    tick  = float(getattr(cal, 'tick', 1.0))

    for _, row in df_year.iterrows():
        t_hor  = row[date_col]
        R_real = float(row[ghi_col])

        K_n    = daily_strike(t_hor, cal)
        Gamma += max(K_n - R_real, 0.0)

        if mode in ('V0', 'both'):
            V0 += sorad_price(K_n, R_t0, str(t0.date()), str(t_hor.date()), cal)

        if mode in ('Vt', 'both'):
            t_prev = t_hor - pd.Timedelta(days=1)
            prev   = df_full[df_full[date_col] == t_prev]
            if len(prev) == 0:
                # Fall back to t0 if previous day not in dataset
                R_prev = R_t0
                t_prev_use = t0
            else:
                R_prev = float(prev[ghi_col].iloc[0])
                t_prev_use = prev[date_col].iloc[0]
            Vt += sorad_price(K_n, R_prev, str(t_prev_use.date()), str(t_hor.date()), cal)

    return {
        'eval_year': eval_year,
        'V0':        tick * V0,
        'Vt':        tick * Vt,
        'Gamma':     tick * Gamma,
        'n_days':    len(df_year),
    }


# ---------------------------------------------------------------------------
# Train-test loop over all evaluation years
# ---------------------------------------------------------------------------

def daily_vt_records(
    eval_year: int,
    cal,
    df_full: pd.DataFrame,
    date_col: str = 'date',
    ghi_col: str = 'GHI',
) -> pd.DataFrame:
    """Return a per-day DataFrame of dynamic Vt prices and payoffs for eval_year.

    Columns: date, day_of_year (1-based), K_n, Vt_n, payoff_n, R_n,
             cum_PL, cum_return_pct  (normalized by Vt_total).

    Used to construct Figure 5 left panel (cumulative net return for seller).
    """
    df_year = df_full[df_full[date_col].dt.year == eval_year].sort_values(date_col).reset_index(drop=True)

    rows = []
    for i, (_, row) in enumerate(df_year.iterrows()):
        t_hor  = row[date_col]
        t_prev = t_hor - pd.Timedelta(days=1)
        prev   = df_full[df_full[date_col] == t_prev]
        if len(prev) == 0:
            prev_fallback = df_full[df_full[date_col] < t_hor].sort_values(date_col).iloc[-1:]
            R_prev = float(prev_fallback[ghi_col].iloc[0])
            t_prev_use = prev_fallback[date_col].iloc[0]
        else:
            R_prev = float(prev[ghi_col].iloc[0])
            t_prev_use = prev[date_col].iloc[0]

        K_n    = daily_strike(t_hor, cal)
        tick   = float(getattr(cal, 'tick', 1.0))
        Vt_n   = tick * sorad_price(K_n, R_prev, str(t_prev_use.date()), str(t_hor.date()), cal)
        pay_n  = tick * max(K_n - float(row[ghi_col]), 0.0)

        rows.append({'date': t_hor, 'day': i + 1, 'K_n': K_n,
                     'Vt_n': Vt_n, 'payoff_n': pay_n, 'R_n': float(row[ghi_col])})

    df = pd.DataFrame(rows)
    Vt_total = df['Vt_n'].sum()
    df['cum_PL']         = (df['Vt_n'] - df['payoff_n']).cumsum()
    df['cum_return_pct'] = df['cum_PL'] / Vt_total * 100.0
    df['Vt_total']       = Vt_total
    return df


def compute_soradidx_table(
    df_full: pd.DataFrame,
    eval_years: Sequence[int] = (2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023),
    coords: dict | None = None,
    target_col: str = 'GHI',
    clearsky_col: str = 'clearsky',
    date_col: str = 'date',
    tick: float = 0.10,
    hist_start_year: int = 2006,
    progress: bool = True,
) -> pd.DataFrame:
    """Run the train-test loop producing Table 1 of the paper.

    For each year h in [hist_start_year, max(eval_years)]:
      - train on 2005..(h-1)
      - compute Gamma_h (single-year realised payoff) using the contemporaneous model
      - if h in eval_years: additionally compute V0, Vt
    Then V_hist[Y] = mean(Gamma_h for h in hist_start_year..Y).

    A final 'Average {hist_start_year}-{max(eval_years)}' row averages each column
    over the prediction years.

    Returns
    -------
    pd.DataFrame  with columns:
        eval_year, train_end_year, V_hist, V0, Vt, Gamma, n_days
    """
    from .calibration import calibrate_window

    if coords is None:
        coords = {'lat': 44.5}

    eval_years = list(eval_years)
    eval_set   = set(eval_years)
    max_eval   = max(eval_years)
    all_years  = list(range(hist_start_year, max_eval + 1))

    records: list[dict] = []
    gamma_history: dict[int, float] = {}

    for year in all_years:
        train_end = year - 1
        if progress:
            print(f"  [{year}] Training on 2005–{train_end}...", end=' ', flush=True)

        df_train = df_full[df_full[date_col].dt.year <= train_end].copy()
        cal = calibrate_window(
            df_train,
            target_col=target_col,
            clearsky_col=clearsky_col,
            date_col=date_col,
            coords=coords,
            train_end_year=train_end,
            tick=tick,
        )

        is_eval = year in eval_set
        if progress:
            print(("pricing " if is_eval else "gamma   ") + f"{year}...", end=' ', flush=True)

        result = _soradidx_year(
            year, cal, df_full,
            mode='both' if is_eval else 'gamma',
            date_col=date_col, ghi_col=target_col,
        )
        gamma_history[year] = result['Gamma']

        if not is_eval:
            if progress:
                print(f"Gamma={result['Gamma']:.3f}")
            continue

        # Running V_hist: mean of Gamma_h for h in [hist_start_year, year]
        running = [gamma_history[h] for h in range(hist_start_year, year + 1)]
        V_hist  = float(np.mean(running))

        records.append({
            'eval_year':      year,
            'train_end_year': train_end,
            'V_hist':         V_hist,
            'V0':             result['V0'],
            'Vt':             result['Vt'],
            'Gamma':          result['Gamma'],
            'n_days':         result['n_days'],
        })

        if progress:
            print(f"V_hist={V_hist:.3f}  V0={result['V0']:.3f}  Vt={result['Vt']:.3f}  Gamma={result['Gamma']:.3f}")

    df_out = pd.DataFrame(records, columns=[
        'eval_year', 'train_end_year', 'V_hist', 'V0', 'Vt', 'Gamma', 'n_days'
    ])

    # Average row (over the prediction years)
    avg_row = {
        'eval_year':      f'Average {hist_start_year}-{max_eval}',
        'train_end_year': '',
        'V_hist':         float(df_out['V_hist'].mean()),
        'V0':             float(df_out['V0'].mean()),
        'Vt':             float(df_out['Vt'].mean()),
        'Gamma':          float(df_out['Gamma'].mean()),
        'n_days':         '',
    }
    df_out = pd.concat([df_out, pd.DataFrame([avg_row])], ignore_index=True)
    return df_out
