"""Phase 3: SoRad (Solar Radiation Derivative) pricing engine.

Implements Eq. 19 from Romagnoli & Sartini (2025):

    V_{i,t}/B(t,T) = (K - C_T·(1-α)) · Φ((K_Y - M_i)/S_i)
                     + C_T·β · G_{i,0}(K_Y; M_i, S_i)

where G_{i,0} = ∫_{-∞}^{K_Y} exp(-exp(y)) · φ(y; M_i, S_i) dy   (Eq. 20)

Total price (mixture over states):
    V_t = p_T · V_{1,t} + (1 - p_T) · V_{0,t}

Under λ^R = 0, Q ≈ P, so V^P = V^Q*.

State convention (matches CalibrationWindowResult.monthly_mixture):
    state 1 = cloudy  (mu1, sd1)
    state 0 = sunny   (mu0, sd0)
    p_T     = P(B_T = 1)
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad
from scipy.stats import norm

from .forecastDensity import clearsky_at, conditional_moments_Y, seasonal_mean_Y_at

OMEGA = 2.0 * np.pi / 365.0


# ---------------------------------------------------------------------------
# Strike transform
# ---------------------------------------------------------------------------

def strike_to_Y(K: float, t_hor, cal) -> float:
    """K_Y = log(log(β) - log(1 - α - K/C_T)).

    Raises ValueError if K is outside the valid GHI range (K ≥ C_T·(1-α)).
    """
    C_T = clearsky_at(t_hor, cal)
    alpha, beta = cal.alpha, cal.beta

    inner = 1.0 - alpha - float(K) / C_T
    if inner <= 0.0:
        raise ValueError(
            f"K={K} >= C_T·(1-α)={C_T*(1-alpha):.4f}: strike above upper GHI bound."
        )
    if inner >= beta:
        raise ValueError(
            f"K={K} <= C_T·(1-α-β)={C_T*(1-alpha-beta):.4f}: strike below lower GHI bound."
        )
    return float(np.log(np.log(beta) - np.log(inner)))


def daily_strike(t_hor, cal) -> float:
    """K_n = C_n · (1 - α - β · exp(-exp(Ȳ_n))).

    The strike that makes the SoRad ATM relative to the seasonal expectation.
    """
    C_T = clearsky_at(t_hor, cal)
    Y_bar = seasonal_mean_Y_at(t_hor, cal)
    alpha, beta = cal.alpha, cal.beta
    return float(C_T * (1.0 - alpha - beta * np.exp(-np.exp(Y_bar))))


# ---------------------------------------------------------------------------
# G integral  (Eq. 20, q=0)
# ---------------------------------------------------------------------------

def G_integral(K_Y: float, mu_i: float, sd_i: float) -> float:
    """∫_{-∞}^{K_Y} exp(-exp(y)) · φ(y; mu_i, sd_i) dy via scipy.integrate.quad."""
    lo = mu_i - 10.0 * sd_i

    def integrand(y: float) -> float:
        return np.exp(-np.exp(y)) * norm.pdf(y, mu_i, sd_i)

    result, _ = quad(integrand, -np.inf, K_Y, limit=200)
    return float(result)


# ---------------------------------------------------------------------------
# Per-state SoRad price
# ---------------------------------------------------------------------------

def _sorad_price_state(K: float, K_Y: float, C_T: float, alpha: float, beta: float,
                       mu_i: float, sd_i: float) -> float:
    """V_{i,t}/B(t,T) for a single state i (Eq. 19)."""
    F_i = norm.cdf((K_Y - mu_i) / sd_i)
    G_i = G_integral(K_Y, mu_i, sd_i)
    return float((K - C_T * (1.0 - alpha)) * F_i + C_T * beta * G_i)


# ---------------------------------------------------------------------------
# Public SoRad price
# ---------------------------------------------------------------------------

def sorad_price(
    K: float,
    R_t,
    t_now,
    t_hor,
    cal,
    p_T: float | None = None,
) -> float:
    """Mixture SoRad price V_t = p_T · V_{1,t} + (1-p_T) · V_{0,t}.

    Parameters
    ----------
    K     : float  — strike GHI (kWh/m²)
    R_t   : float  — observed GHI at t_now (conditioning value)
    t_now : date-like  — conditioning date
    t_hor : date-like  — forecast horizon (delivery date)
    cal   : CalibrationWindowResult
    p_T   : P(B_T = 1); defaults to the monthly prior from cal

    Returns
    -------
    float — SoRad price (same units as GHI × tick convention)
    """
    mom  = conditional_moments_Y(R_t, t_now, t_hor, cal, p_T)
    C_T  = clearsky_at(t_hor, cal)
    alpha, beta = cal.alpha, cal.beta
    p    = mom['p_T']

    K_Y  = strike_to_Y(K, t_hor, cal)

    sd1 = float(np.sqrt(max(mom['var1'], 1e-12)))
    sd0 = float(np.sqrt(max(mom['var0'], 1e-12)))

    V1 = _sorad_price_state(K, K_Y, C_T, alpha, beta, mom['mu1'], sd1)
    V0 = _sorad_price_state(K, K_Y, C_T, alpha, beta, mom['mu0'], sd0)

    return float(p * V1 + (1.0 - p) * V0)
