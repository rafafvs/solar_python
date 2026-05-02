"""
Solar model helper functions — Python translation.

Functions:
    create_monthly_sequence       — monthly date grid between t_now and t_hor
    martingale_method_seasonal    — mean-reversion θ via martingale estimator
    reparam_seasonal_function     — reparametrise seasonal variance to continuous time
    integral_sigma_numeric        — numerical integral of √(σ̄(s) · e^{-θ(T-s)})
    integral_sigma2_formula       — closed-form integral of σ̄²(s) · e^{-2θ(T-s)}

Version: 1.0.0
"""

from __future__ import annotations

import math
import numpy as np
import pandas as pd
from scipy import integrate as sci_integrate


# ---------------------------------------------------------------------------
# number_of_day  (helper used throughout — mirrors R's number_of_day)
# ---------------------------------------------------------------------------

def number_of_day(date) -> int:
    """Return the day-of-year (1–365/366) for a date or date-like string."""
    return pd.Timestamp(date).day_of_year


# ---------------------------------------------------------------------------
# create_monthly_sequence
# ---------------------------------------------------------------------------

def create_monthly_sequence(
    t_now,
    t_hor,
    last_day: bool = False,
) -> pd.DataFrame:
    """
    Build a monthly date grid between ``t_now`` and ``t_hor``.

    Each row covers one calendar month (or partial month at the boundaries).
    The grid is used to integrate piecewise-constant monthly parameters over
    the horizon [t_now, t_hor].

    Parameters
    ----------
    t_now    : date-like, start date
    t_hor    : date-like, horizon (end) date
    last_day : bool
        When True, the last row is split into two rows: one covering
        [t_{T-1}, T-1] (shared variance) and one for the final day
        [T-1, T] (conditional variance).

    Returns
    -------
    pd.DataFrame with columns:
        Year   — calendar year of the sub-interval start
        Month  — calendar month of the sub-interval end
        n      — start day index (days since t_now, with seasonal offset)
        N      — end   day index
        tau    — total time to maturity (days from t_now to t_hor)

    The ``last_day`` flag is stored as ``df.attrs["last_day"]``.

    Examples
    --------
    >>> create_monthly_sequence("2022-01-01", "2022-03-24")
    >>> create_monthly_sequence("2022-01-01", "2022-03-24", last_day=True)
    """
    t_now = pd.Timestamp(t_now)
    t_hor = pd.Timestamp(t_hor)

    # First day of the start and end months
    t_start = t_now.replace(day=1)
    t_end   = t_hor.replace(day=1)

    # Monthly sequence of first-of-month dates
    month_starts = pd.date_range(t_start, t_end, freq="MS")
    n_months = len(month_starts)

    # dates_now: start of each sub-interval
    #   first entry = t_now itself;
    #   subsequent = last day of each month  (first-of-month minus its own day)
    dates_now = [t_now] + [
        m - pd.Timedelta(days=m.day)          # last day of previous month
        for m in month_starts[1:]
    ]

    # dates_hor: end of each sub-interval
    #   all-but-last = last day of that month;
    #   last entry = t_hor itself
    dates_hor = [
        m + pd.offsets.MonthEnd(0)            # last day of the month
        for m in month_starts[:-1]
    ] + [t_hor]

    # Number of days in each sub-interval
    n_of_day = np.array([
        (h - s).days
        for s, h in zip(dates_now, dates_hor)
    ], dtype=float)

    # Build the base DataFrame
    df = pd.DataFrame({
        "Year":  [d.year  for d in dates_now],
        "Month": [d.month for d in dates_hor],
        "n":     n_of_day,
    })

    # Cumulative day indices (N = end day index measured from the seasonal origin)
    start_doy = number_of_day(t_now)               # day-of-year of t_now
    cumulative_lag = np.concatenate([[0], np.cumsum(n_of_day[:-1])])
    df["N"] = df["n"] + start_doy + cumulative_lag
    df["n"] = df["N"] - df["n"]                    # n = start index of sub-interval

    # Time-to-maturity (total days from t_now to t_hor, measured from n[0])
    total_tau = (t_hor - t_now).days
    df["tau"] = total_tau + df["n"].iloc[0]

    # Optional last-day split
    if last_day:
        last_row = df.iloc[[-1]].copy()
        last_row["n"] = last_row["tau"] - 1        # penultimate day → last sub-interval

        # Trim the end of the second-to-last row by 1 day
        df.iloc[-1, df.columns.get_loc("N")] = df.iloc[-1]["N"] - 1

        df = pd.concat([df, last_row], ignore_index=True)

    df.attrs["last_day"] = last_day
    return df


# ---------------------------------------------------------------------------
# martingale_method_seasonal
# ---------------------------------------------------------------------------

def martingale_method_seasonal(
    Yt: pd.Series | np.ndarray,
    Yt_bar: pd.Series | np.ndarray,
    e_mu: float | pd.Series | np.ndarray = 0.0,
) -> float:
    """
    Estimate the mean-reversion parameter θ via the martingale method.

    The estimator is:
        a_n = Σ w_t · dY_t / Σ w_t · dY_{t-1}
        θ   = −log(a_n)

    where the instrument w_t = (Ȳ_{t-1} − Y_{t-1}) / σ̄²_{t-1}
    and σ̄²_t is estimated from the quadratic variation of Y_t.

    Parameters
    ----------
    Yt     : array-like, transformed solar radiation time series
    Yt_bar : array-like, seasonal mean of Y_t
    e_mu   : float or array-like, expected drift (default 0)

    Returns
    -------
    float, estimated θ
    """
    Yt     = np.asarray(Yt,     dtype=float)
    Yt_bar = np.asarray(Yt_bar, dtype=float)
    e_mu   = np.broadcast_to(np.asarray(e_mu, dtype=float), Yt.shape).copy()

    # Quadratic variation: (Y_{t-1} - Y_{t-2})²
    dYt2 = (np.roll(Yt, 1) - np.roll(Yt, 2)) ** 2
    dYt2[:2] = np.nan

    # Fit a simple seasonal model (constant mean) to the quadratic variation
    # Mirrors seasonalModel$new()$fit("dYt2 ~ 1", …) — intercept-only model
    n = np.arange(len(Yt))
    valid = ~np.isnan(dYt2)
    sigma2_bar = np.full(len(Yt), np.nanmean(dYt2))   # constant seasonal variance

    # Instrument
    Yt_bar_L1 = np.roll(Yt_bar, 1)
    Yt_L1     = np.roll(Yt,     1)
    Y_est_L1  = (Yt_bar_L1 - Yt_L1) / sigma2_bar

    # Differences from seasonal mean
    dYt    = Yt - Yt_bar
    dYt_L1 = np.roll(Yt, 1) - np.roll(Yt_bar, 1) - np.roll(e_mu, 1)

    # Build DataFrame and drop NaN rows (mirrors na.omit)
    df = pd.DataFrame({
        "Y_est_L1": Y_est_L1,
        "dYt":      dYt,
        "dYt_L1":   dYt_L1,
    })
    df = df.dropna()

    a_n = (df["Y_est_L1"] * df["dYt"]).sum() / (df["Y_est_L1"] * df["dYt_L1"]).sum()
    return float(-math.log(a_n))


# ---------------------------------------------------------------------------
# reparam_seasonal_function
# ---------------------------------------------------------------------------

def reparam_seasonal_function(
    par: np.ndarray | dict,
    theta: float,
    omega: float = 2 * math.pi / 365,
) -> dict:
    """
    Reparametrise the OLS seasonal variance coefficients (a0, a1, a2) to
    continuous-time parameters (c0, c1, c2) and their integral counterparts
    (γ0, γ1, γ2), for both the short-term and long-term formulations.

    The seasonal variance function is:
        σ̄²(t) = a0 + a1·sin(ω·t) + a2·cos(ω·t)

    Parameters
    ----------
    par   : array-like of length 3, [a0, a1, a2]
    theta : float, mean-reversion parameter
    omega : float, angular frequency (default 2π/365)

    Returns
    -------
    dict with keys:
        alpha, beta, detM           — transformation intermediates
        a_    — original [a0, a1, a2]
        c_    — short-term [c0, c1, c2]
        c_long— long-term  [c0, c1, c2]
        gamma — short-term integral parameters [γ0, γ1, γ2]
        gamma_long — long-term integral parameters
    """
    if hasattr(par, "values"):
        par = par.values
    par = np.asarray(par, dtype=float).ravel()
    a0, a1, a2 = par[0], par[1], par[2]

    # ---- Long-term (steady-state) reparametrisation ----
    c0_long = a0 * 2 * theta
    c1_long = a1 * 2 * theta - omega * a2
    c2_long = a2 * 2 * theta + omega * a1

    gamma0_long =  c0_long / (2 * theta)
    gamma1_long = (c1_long * 2 * theta + c2_long * omega) / (4 * theta**2 + omega**2)
    gamma2_long = (c2_long * 2 * theta - c1_long * omega) / (4 * theta**2 + omega**2)

    # ---- Short-term reparametrisation ----
    alpha = 1 - math.exp(-2 * theta) * math.cos(omega)
    beta  =     math.exp(-2 * theta) * math.sin(omega)
    detM  = alpha**2 + beta**2

    c0 = (2 * theta * a0) / (1 - math.exp(-2 * theta))
    c1 = ((2 * theta * alpha + omega * beta) * a1 + (2 * theta * beta  - omega * alpha) * a2) / detM
    c2 = ((omega * alpha - 2 * theta * beta) * a1 + (omega * beta  + 2 * theta * alpha) * a2) / detM

    gamma0 =  c0 / (2 * theta)
    gamma1 = (c1 * 2 * theta + c2 * omega) / (4 * theta**2 + omega**2)
    gamma2 = (c2 * 2 * theta - c1 * omega) / (4 * theta**2 + omega**2)

    return {
        "alpha":      alpha,
        "beta":       beta,
        "detM":       detM,
        "a_":         pd.Series({"a0": a0, "a1": a1, "a2": a2}),
        "c_":         pd.Series({"c0": c0, "c1": c1, "c2": c2}),
        "c_long":     pd.Series({"c0": c0_long, "c1": c1_long, "c2": c2_long}),
        "gamma":      pd.Series({"gamma0": gamma0, "gamma1": gamma1, "gamma2": gamma2}),
        "gamma_long": pd.Series({"gamma0": gamma0_long, "gamma1": gamma1_long, "gamma2": gamma2_long}),
    }


# ---------------------------------------------------------------------------
# integral_sigma_numeric
# ---------------------------------------------------------------------------

def integral_sigma_numeric(
    theta: float,
    par: np.ndarray | pd.Series,
    omega: float = 2 * math.pi / 365,
):
    """
    Return a callable that numerically integrates
        ∫_t^s √(σ̄(τ) · e^{−2θ(T−τ)}) dτ

    where σ̄(τ) = par[0] + par[1]·sin(ω·τ) + par[2]·cos(ω·τ).

    Parameters
    ----------
    theta : float, mean-reversion parameter
    par   : array-like of length 3, [c0, c1, c2]
    omega : float

    Returns
    -------
    callable f(t, s, T_) → np.ndarray
        t, s, T_ are array-like of equal length (or scalars).
        Returns the integral value for each (t_i, s_i, T_i) triple.

    Examples
    --------
    >>> fn = integral_sigma_numeric(0.1, [0.5, 0.1, 0.05])
    >>> fn([0], [30], [30])
    """
    if hasattr(par, "values"):
        par = par.values
    par = np.asarray(par, dtype=float).ravel()

    def seasonal_function(t: float) -> float:
        return par[0] + par[1] * math.sin(omega * t) + par[2] * math.cos(omega * t)

    def integrand(tau: float, T_: float) -> float:
        return math.sqrt(max(seasonal_function(tau) * math.exp(-2 * theta * (T_ - tau)), 0.0))

    def fn(
        t:  np.ndarray | list | float,
        s:  np.ndarray | list | float,
        T_: np.ndarray | list | float,
    ) -> np.ndarray:
        t  = np.atleast_1d(np.asarray(t,  dtype=float))
        s  = np.atleast_1d(np.asarray(s,  dtype=float))
        T_ = np.atleast_1d(np.asarray(T_, dtype=float))
        result = np.empty(len(T_))
        for i in range(len(T_)):
            val, _ = sci_integrate.quad(integrand, t[i], s[i], args=(T_[i],))
            result[i] = val
        return result

    return fn


# ---------------------------------------------------------------------------
# integral_sigma2_formula
# ---------------------------------------------------------------------------

def integral_sigma2_formula(
    theta: float,
    par: np.ndarray | pd.Series,
    omega: float = 2 * math.pi / 365,
):
    """
    Return a callable that computes the closed-form integral
        ∫_t^s σ̄²(τ) · e^{−2θ(T−τ)} dτ

    using the formula:
        result = par[0]·f0(t,s,T) + par[1]·f1(t,s,T) + par[2]·f2(t,s,T)

    where:
        f0(t,s,T) = e^{−2θ(T−s)} − e^{−2θ(T−t)}
        f1(t,s,T) = e^{−2θ(T−s)}·sin(ω·s) − e^{−2θ(T−t)}·sin(ω·t)
        f2(t,s,T) = e^{−2θ(T−s)}·cos(ω·s) − e^{−2θ(T−t)}·cos(ω·t)

    Parameters
    ----------
    theta : float, mean-reversion parameter
    par   : array-like of length 3, [γ0, γ1, γ2]
    omega : float

    Returns
    -------
    callable f(t, s, T_) → np.ndarray
        All arguments may be scalars or array-like of equal length.

    Examples
    --------
    >>> fn = integral_sigma2_formula(0.1, [0.5, 0.1, 0.05])
    >>> fn([0], [30], [30])
    """
    if hasattr(par, "values"):
        par = par.values
    par = np.asarray(par, dtype=float).ravel()

    def fn(
        t:  np.ndarray | list | float,
        s:  np.ndarray | list | float,
        T_: np.ndarray | list | float,
    ) -> np.ndarray:
        t  = np.asarray(t,  dtype=float)
        s  = np.asarray(s,  dtype=float)
        T_ = np.asarray(T_, dtype=float)

        exp_s = np.exp(-2 * theta * (T_ - s))
        exp_t = np.exp(-2 * theta * (T_ - t))

        f0 = exp_s - exp_t
        f1 = exp_s * np.sin(omega * s) - exp_t * np.sin(omega * t)
        f2 = exp_s * np.cos(omega * s) - exp_t * np.cos(omega * t)

        return par[0] * f0 + par[1] * f1 + par[2] * f2

    return fn