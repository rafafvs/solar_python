"""Phase 2: Euler-Maruyama simulation of the OU-mixture solar radiation SDE.

Implements the discrete-time approximation (Eq. 3) of:

    dY_t = θ(Ȳ(t) - Y_t) dt + σ̄(t) η_{B_t} dt + σ̄(t) σ_{B_t} dW_t

where B_t ~ Bernoulli(p_{month(t)}) and η | B_t ~ N(μ_{B_t}, σ²_{B_t}).

The Euler-Maruyama step with sub-day increment dt is:
    Y_{t+dt} = Y_t + θ(Ȳ(t) - Y_t)·dt + σ̄(t)·noise_t·√dt
    noise_t  ~ N(μ_{B_t, month(t)}, σ²_{B_t, month(t)})

After simulation in Y-space, the result is pushed back to GHI via
    R_t = C_t·(1 - α - β·exp(-exp(Y_t)))
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .forecastDensity import clearsky_array, seasonal_mean_Y_at, Y_to_R, OMEGA
from .zzz import number_of_day


def simulate_ou_mixture(
    cal,
    t_start: str,
    t_end: str,
    dt: float = 0.05,
    n_paths: int = 1,
    seed: int | None = None,
) -> pd.DataFrame:
    """
    Euler-Maruyama simulation of Y_t and R_t over [t_start, t_end].

    Parameters
    ----------
    cal     : CalibrationWindowResult  — calibration parameters
    t_start : str  — simulation start date (e.g. '2016-01-01')
    t_end   : str  — simulation end   date (e.g. '2022-12-31')
    dt      : float — sub-day time step in day units (default 0.05)
    n_paths : int   — number of independent paths to simulate
    seed    : int | None — RNG seed for reproducibility

    Returns
    -------
    pd.DataFrame with columns
        date    : actual date (only at integer-day boundaries)
        doy     : day-of-year
        Y_sim   : simulated Y_t  (one column per path, 'Y_sim_0' … 'Y_sim_{n_paths-1}')
        R_sim   : simulated R_t  ('R_sim_0' … 'R_sim_{n_paths-1}')
        C_t     : clear-sky value at that day
    """
    rng = np.random.default_rng(seed)

    start = pd.Timestamp(t_start)
    end   = pd.Timestamp(t_end)

    # All calendar days in the range (integer days only — we store results here)
    all_dates = pd.date_range(start, end, freq='D')
    n_days    = len(all_dates)
    steps_per_day = max(1, int(round(1.0 / dt)))

    # Precompute per-day quantities (on the integer-day grid)
    doys    = all_dates.dayofyear.astype(float).values   # shape (n_days,)
    months  = all_dates.month.values                     # shape (n_days,)

    C_arr   = clearsky_array(doys, cal)                  # C_t  shape (n_days,)

    a0, a1, a2 = cal.a
    Ybar_arr = a0 + a1 * np.sin(OMEGA * doys) + a2 * np.cos(OMEGA * doys)  # Ȳ_t

    c0, c1, c2 = cal.c
    sigma_arr = np.sqrt(
        np.maximum(c0 + c1 * np.sin(OMEGA * doys) + c2 * np.cos(OMEGA * doys), 0.0)
    )  # σ̄_t  shape (n_days,)

    # Monthly mixture parameters as lookup arrays (index = month 1..12)
    mm = cal.monthly_mixture.set_index('Month')
    p_month  = mm['p'].values      # shape (12,) — index 0 → Month 1
    mu1_month = mm['mu1'].values
    mu0_month = mm['mu0'].values
    sd1_month = mm['sd1'].values
    sd0_month = mm['sd0'].values

    # ── Simulation ──────────────────────────────────────────────────────────
    # Initial Y_0: use Ȳ at t_start
    Y0 = Ybar_arr[0]
    Y_cur = np.full(n_paths, Y0, dtype=float)   # current Y for each path

    # Store results at integer-day boundaries
    Y_out = np.zeros((n_days, n_paths), dtype=float)
    Y_out[0] = Y_cur

    for day_idx in range(n_days - 1):
        # Day-level quantities (constant within the day)
        Y_bar_d  = Ybar_arr[day_idx]
        sigma_d  = sigma_arr[day_idx]
        month_d  = months[day_idx] - 1   # 0-based index into month arrays
        p_d      = p_month[month_d]
        theta    = cal.theta

        # Sub-day steps
        for _ in range(steps_per_day):
            # Draw regime B for each path  (1 = cloudy, 0 = sunny)
            B = (rng.random(n_paths) < p_d).astype(int)

            # Draw noise η | B ~ N(μ_B, σ²_B)
            mu_noise  = np.where(B == 1, mu1_month[month_d], mu0_month[month_d])
            sd_noise  = np.where(B == 1, sd1_month[month_d], sd0_month[month_d])
            eta       = rng.normal(mu_noise, sd_noise)

            # Euler-Maruyama step
            Y_cur = (
                Y_cur
                + theta * (Y_bar_d - Y_cur) * dt
                + sigma_d * eta * np.sqrt(dt)
            )

        Y_out[day_idx + 1] = Y_cur

    # ── Push back to R (GHI) ──────────────────────────────────────────────
    alpha, beta = cal.alpha, cal.beta
    R_out = C_arr[:, None] * (1.0 - alpha - beta * np.exp(-np.exp(Y_out)))

    # ── Assemble output DataFrame ─────────────────────────────────────────
    df = pd.DataFrame({'date': all_dates, 'doy': doys.astype(int), 'C_t': C_arr})
    for k in range(n_paths):
        df[f'Y_sim_{k}'] = Y_out[:, k]
        df[f'R_sim_{k}'] = R_out[:, k]

    return df
