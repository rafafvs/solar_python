"""Phase 2: conditional moments and forecast density on the GHI scale.

Standalone functions that operate directly on CalibrationWindowResult
(from calibration.py) without going through RadiationModel / SolarModel.

Math reference: Romagnoli & Sartini (2025), Appendix B.
  - Eq. B.6/B.7: conditional mean M_{Y_i}(t, T)
  - Eq. B.16:    conditional variance S^2_{Y_i}(t, T)
  - Eq. 13:      mixture density/CDF of R_T

State convention (matches CalibrationWindowResult.monthly_mixture):
  state 1 = cloudy  (lower mean, columns mu1, sd1)
  state 0 = sunny   (higher mean, columns mu0, sd0)
  p       = P(B_T = 1) = probability of cloudy state
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import norm

from .radiationModel_internals import (
    create_monthly_sequence,
    integral_sigma_numeric,
    integral_sigma2_formula,
)
from .zzz import number_of_day

OMEGA = 2.0 * np.pi / 365.0


# ---------------------------------------------------------------------------
# Point-in-time helpers
# ---------------------------------------------------------------------------

def clearsky_at(t, cal) -> float:
    """C_t = clear-sky radiation at day t (date string, Timestamp, or day-of-year)."""
    doy = number_of_day(t)
    vals = cal.clearsky_model.predict(n=float(doy))
    return float(np.atleast_1d(vals)[0])


def clearsky_array(doys: np.ndarray, cal) -> np.ndarray:
    """Vectorised C_t for an array of day-of-year values."""
    vals = cal.clearsky_model.predict(n=doys.astype(float))
    return np.asarray(vals, dtype=float)


def seasonal_mean_Y_at(t, cal) -> float:
    """Ȳ_t = a0 + a1·sin(ωt) + a2·cos(ωt) at day t."""
    doy = float(number_of_day(t))
    a0, a1, a2 = cal.a
    return float(a0 + a1 * np.sin(OMEGA * doy) + a2 * np.cos(OMEGA * doy))


def R_to_Y(R: float, C_t: float, alpha: float, beta: float) -> float:
    """Inverse transform: Y_t = log(log(β) - log(1 - α - R/C_t))."""
    inner = 1.0 - alpha - np.asarray(R, dtype=float) / C_t
    return float(np.log(np.log(beta) - np.log(inner)))


def Y_to_R(Y: float, C_t: float, alpha: float, beta: float) -> float:
    """Forward transform: R_t = C_t·(1 - α - β·exp(-exp(Y)))."""
    return float(C_t * (1.0 - alpha - beta * np.exp(-np.exp(Y))))


# ---------------------------------------------------------------------------
# Conditional moments of Y_T  (Eq. B.6, B.7, B.16)
# ---------------------------------------------------------------------------

def conditional_moments_Y(
    R_t,
    t_now,
    t_hor,
    cal,
    p_T: float | None = None,
) -> dict:
    """
    Conditional mean and variance of Y_T given (Y_t, B_T = i).

    Parameters
    ----------
    R_t   : float  — observed GHI at t_now
    t_now : date-like  — conditioning date
    t_hor : date-like  — forecast horizon
    cal   : CalibrationWindowResult
    p_T   : P(B_T = 1) at the horizon; defaults to the monthly prior

    Returns
    -------
    dict with keys
        mu1  : E[Y_T | Y_t, B_T=1]   (state 1, cloudy)
        mu0  : E[Y_T | Y_t, B_T=0]   (state 0, sunny)
        var1 : Var[Y_T | Y_t, B_T=1]
        var0 : Var[Y_T | Y_t, B_T=0]
        mu_e : OU forecast, no mixture drift
        p_T  : P(B_T = 1) used in this call
    """
    t_now = pd.Timestamp(t_now)
    t_hor = pd.Timestamp(t_hor)

    # Month-boundary sub-interval grid (last row = day T alone)
    df_date = create_monthly_sequence(t_now, t_hor, last_day=True)

    # Merge monthly mixture parameters on Month
    mm = cal.monthly_mixture[['Month', 'p', 'mu1', 'sd1', 'mu0', 'sd0']].copy()
    df = df_date.merge(mm, on='Month', how='left').sort_values('n').reset_index(drop=True)

    # Integral callables  ─  I(t,s,T) and J(t,s,T)
    int_I = integral_sigma_numeric(cal.theta, cal.c)
    int_J = integral_sigma2_formula(cal.theta, cal.gamma)

    n_arr = df['n'].values.astype(float)
    N_arr = df['N'].values.astype(float)
    tau_arr = df['tau'].values.astype(float)

    df = df.copy()
    df['I_sub'] = int_I(n_arr, N_arr, tau_arr)
    df['J_sub'] = int_J(n_arr, N_arr, tau_arr)

    body = df.iloc[:-1]
    last = df.iloc[-1]

    # ── OU mean-reversion forecast (no mixture drift) ──
    C_now     = clearsky_at(t_now, cal)
    # Clamp R_t to the open interior of the valid GHI range.
    # The dataset clearsky may differ from the calibration model's clearsky,
    # so observed GHI can occasionally exceed the model's upper bound.
    _eps   = 1e-6 * C_now
    _R_t   = float(np.clip(float(R_t),
                            C_now * (1.0 - cal.alpha - cal.beta) + _eps,
                            C_now * (1.0 - cal.alpha) - _eps))
    Y_t       = R_to_Y(_R_t, C_now, cal.alpha, cal.beta)
    tau       = (t_hor - t_now).days
    Y_bar_now = seasonal_mean_Y_at(t_now, cal)
    Y_bar_hor = seasonal_mean_Y_at(t_hor, cal)
    Y_forecast = Y_bar_hor + (Y_t - Y_bar_now) * np.exp(-cal.theta * tau)

    # ── Mixture drift integral (Eq. B.6/B.7) ──
    # Body rows: B_τ unknown → use unconditional mean e_μ_B = p·μ1 + (1-p)·μ0
    e_mu_B    = body['mu1'].values * body['p'].values + body['mu0'].values * (1.0 - body['p'].values)
    sum_drift = float((body['I_sub'].values * e_mu_B).sum())

    mu1_Y = Y_forecast + sum_drift + float(last['I_sub'] * last['mu1'])
    mu0_Y = Y_forecast + sum_drift + float(last['I_sub'] * last['mu0'])

    # ── Mixture variance integral (Eq. B.16) ──
    mu_diff  = body['mu1'].values - body['mu0'].values
    e_sd2_B  = (body['sd1'].values**2 * body['p'].values
                + body['sd0'].values**2 * (1.0 - body['p'].values))

    var_drift  = float((body['J_sub'].values * mu_diff**2
                        * body['p'].values * (1.0 - body['p'].values)).sum())
    var_diff   = float((body['J_sub'].values * e_sd2_B).sum())
    common_var = var_drift + var_diff

    J_last = float(last['J_sub'])
    var1_Y = common_var + J_last * float(last['sd1'])**2
    var0_Y = common_var + J_last * float(last['sd0'])**2

    # Default p_T: unconditional monthly prior for month of t_hor
    if p_T is None:
        m    = t_hor.month
        rows = mm.loc[mm['Month'] == m, 'p']
        p_T  = float(rows.iloc[0]) if len(rows) else 0.5

    return {
        'mu1':  mu1_Y,
        'mu0':  mu0_Y,
        'var1': var1_Y,
        'var0': var0_Y,
        'mu_e': Y_forecast,
        'p_T':  p_T,
    }


# ---------------------------------------------------------------------------
# Internal helpers that work with a precomputed moments dict
# ---------------------------------------------------------------------------

def _moments_to_Y_bracket(mom: dict) -> tuple[float, float]:
    """[Y_lo, Y_hi] bracket covering ±8 SD around both component means."""
    sd1 = np.sqrt(max(mom['var1'], 1e-12))
    sd0 = np.sqrt(max(mom['var0'], 1e-12))
    Y_lo = min(mom['mu1'] - 8 * sd1, mom['mu0'] - 8 * sd0)
    Y_hi = max(mom['mu1'] + 8 * sd1, mom['mu0'] + 8 * sd0)
    return Y_lo, Y_hi


def _mixture_cdf_Y(Y: float, mom: dict) -> float:
    """Mixture CDF at Y using precomputed moments (in Y-space)."""
    p    = mom['p_T']
    sd1  = np.sqrt(max(mom['var1'], 1e-12))
    sd0  = np.sqrt(max(mom['var0'], 1e-12))
    return p * norm.cdf(Y, mom['mu1'], sd1) + (1.0 - p) * norm.cdf(Y, mom['mu0'], sd0)


def _mixture_pdf_Y(Y: float, mom: dict) -> float:
    """Mixture PDF at Y using precomputed moments (in Y-space)."""
    p    = mom['p_T']
    sd1  = np.sqrt(max(mom['var1'], 1e-12))
    sd0  = np.sqrt(max(mom['var0'], 1e-12))
    return p * norm.pdf(Y, mom['mu1'], sd1) + (1.0 - p) * norm.pdf(Y, mom['mu0'], sd0)


# ---------------------------------------------------------------------------
# Public density / CDF / quantile / expectation  (on the GHI scale)
# ---------------------------------------------------------------------------

def mixture_pdf_R(R, R_t, t_now, t_hor, cal, p_T: float | None = None) -> float:
    """
    Mixture density of R_T = R, conditioned on R_t at t_now.

    Valid R range: (C_T·(1-α-β), C_T·(1-α)).
    Outside that range the transform Y = h⁻¹(R) is undefined → density = 0.
    """
    mom  = conditional_moments_Y(R_t, t_now, t_hor, cal, p_T)
    C_T  = clearsky_at(t_hor, cal)
    alpha, beta = cal.alpha, cal.beta

    R    = float(R)
    inner = 1.0 - alpha - R / C_T   # = β·exp(-exp(Y)) in the forward map
    # Valid iff 0 < inner < beta  (boundary: inner>=beta → Y→-∞; inner<=0 → Y→+∞)
    if inner <= 0.0 or inner >= beta:
        return 0.0

    Y        = np.log(np.log(beta) - np.log(inner))
    jacobian = 1.0 / (C_T * inner * np.exp(Y))   # |dY/dR|
    return _mixture_pdf_Y(Y, mom) * jacobian


def mixture_cdf_R(R, R_t, t_now, t_hor, cal, p_T: float | None = None) -> float:
    """Mixture CDF of R_T at value R, conditioned on R_t.

    CDF = 0 for R <= R_min = C_T·(1-α-β),
    CDF = 1 for R >= R_max = C_T·(1-α).
    """
    mom  = conditional_moments_Y(R_t, t_now, t_hor, cal, p_T)
    C_T  = clearsky_at(t_hor, cal)
    alpha, beta = cal.alpha, cal.beta

    R    = float(R)
    inner = 1.0 - alpha - R / C_T
    if inner <= 0.0:    # R above upper bound → CDF = 1
        return 1.0
    if inner >= beta:   # R below lower bound → CDF = 0
        return 0.0

    Y = np.log(np.log(beta) - np.log(inner))
    return _mixture_cdf_Y(Y, mom)


def mixture_quantile_R(q: float, R_t, t_now, t_hor, cal, p_T: float | None = None) -> float:
    """q-th quantile of R_T via brentq on the mixture CDF."""
    mom  = conditional_moments_Y(R_t, t_now, t_hor, cal, p_T)
    C_T  = clearsky_at(t_hor, cal)
    alpha, beta = cal.alpha, cal.beta

    Y_lo, Y_hi = _moments_to_Y_bracket(mom)
    R_lo = Y_to_R(Y_lo, C_T, alpha, beta)
    R_hi = Y_to_R(Y_hi, C_T, alpha, beta)

    def _cdf(R_):
        inner = 1.0 - alpha - R_ / C_T
        if inner <= 0.0:
            return 1.0
        if inner >= 1.0 - alpha:
            return 0.0
        Y_ = np.log(np.log(beta) - np.log(inner))
        return _mixture_cdf_Y(Y_, mom) - q

    return brentq(_cdf, R_lo, R_hi, xtol=1e-6)


def mixture_expected_R(R_t, t_now, t_hor, cal, p_T: float | None = None) -> float:
    """E[R_T] via numerical integration of R · f(R)."""
    mom  = conditional_moments_Y(R_t, t_now, t_hor, cal, p_T)
    C_T  = clearsky_at(t_hor, cal)
    alpha, beta = cal.alpha, cal.beta

    Y_lo, Y_hi = _moments_to_Y_bracket(mom)
    R_lo = Y_to_R(Y_lo, C_T, alpha, beta)
    R_hi = Y_to_R(Y_hi, C_T, alpha, beta)

    def integrand(R_):
        inner = 1.0 - alpha - R_ / C_T
        if inner <= 0.0 or inner >= beta:
            return 0.0
        Y_       = np.log(np.log(beta) - np.log(inner))
        jacobian = 1.0 / (C_T * inner * np.exp(Y_))
        return R_ * _mixture_pdf_Y(Y_, mom) * jacobian

    result, _ = quad(integrand, R_lo, R_hi, limit=200)
    return result
