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
    mode      : 'V0', 'Vt', or 'both'
    date_col  : column name for dates
    ghi_col   : column name for observed GHI

    Returns
    -------
    dict with keys: eval_year, V0, Vt, V_hist, n_days
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
    V_hist = 0.0

    for _, row in df_year.iterrows():
        t_hor  = row[date_col]
        R_real = float(row[ghi_col])

        K_n    = daily_strike(t_hor, cal)
        payoff = max(K_n - R_real, 0.0)
        V_hist += payoff

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
        'V0':        V0,
        'Vt':        Vt,
        'V_hist':    V_hist,
        'n_days':    len(df_year),
    }


# ---------------------------------------------------------------------------
# Train-test loop over all evaluation years
# ---------------------------------------------------------------------------

def compute_soradidx_table(
    df_full: pd.DataFrame,
    eval_years: Sequence[int] = (2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023),
    coords: dict | None = None,
    target_col: str = 'GHI',
    clearsky_col: str = 'clearsky',
    date_col: str = 'date',
    progress: bool = True,
) -> pd.DataFrame:
    """Run the train-test loop producing Table 1 of the paper.

    For each eval_year Y:
      - train on 2005..(Y-1)
      - compute V0, Vt, V_hist for year Y

    Parameters
    ----------
    df_full    : full dataset (all years)
    eval_years : evaluation years (paper: 2014–2023)
    coords     : {'lat': ...} passed to calibrate_window
    progress   : print progress messages

    Returns
    -------
    pd.DataFrame  with columns: eval_year, train_end_year, V0, Vt, V_hist, n_days
    """
    from .calibration import calibrate_window

    if coords is None:
        coords = {'lat': 44.5}

    records = []
    for eval_year in eval_years:
        train_end = eval_year - 1
        if progress:
            print(f"  [{eval_year}] Training on 2005–{train_end}...", end=' ', flush=True)

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
            print(f"pricing {eval_year}...", end=' ', flush=True)

        result = _soradidx_year(eval_year, cal, df_full, mode='both',
                                date_col=date_col, ghi_col=target_col)
        result['train_end_year'] = train_end
        records.append(result)

        if progress:
            print(f"V0={result['V0']:.3f}  Vt={result['Vt']:.3f}  V_hist={result['V_hist']:.3f}")

    df_out = pd.DataFrame(records, columns=[
        'eval_year', 'train_end_year', 'V0', 'Vt', 'V_hist', 'n_days'
    ])
    return df_out
