"""
ARMA helper functions — Python translation of the R source.

Functions included:
    ARMA_vector_b          — state-space selection vector
    ARMA_companion_matrix  — companion / transition matrix
    ARMA_expectation       — h-step conditional mean
    ARMA_variance          — h-step forecast variance
    ARMA_covariance        — h-step forecast covariance
    ARMA_next_step         — iterate the state forward h steps
    ARMA_filter            — one-step-ahead filter (replaces the C backend)
    ARMA_forecast          — full state trajectory
    AR_variance            — long-run AR variance (closed-form for p≤4)
    AR_params_to_zeta      — phi → unconstrained zeta  (step-down Durbin-Levinson + atanh)
    ARMA_params_to_zeta    — ARMA wrapper for AR_params_to_zeta
    AR_params_to_phi       — zeta → constrained phi   (tanh + step-up Durbin-Levinson + Jacobian)
    ARMA_params_to_phi     — ARMA wrapper for AR_params_to_phi
    ARMA_loss_logLik       — Gaussian log-likelihood (constrained parametrisation)
    ARMA_loss_CSS          — conditional sum-of-squares loss
    ARMA_fit               — full estimation pipeline

Version: 1.0.0
"""

from __future__ import annotations

import warnings
import numpy as np
import pandas as pd
from scipy import optimize
from scipy.linalg import solve   # solve(A, b) → A^{-1} b


# ---------------------------------------------------------------------------
# ARMA_vector_b
# ---------------------------------------------------------------------------

def ARMA_vector_b(ar_order: int, ma_order: int) -> np.ndarray:
    """
    Build the selection (input) vector b for the ARMA state-space model.

    The state vector is laid out as
        [y_t, y_{t-1}, …, y_{t-p+1},  eps_t, eps_{t-1}, …, eps_{t-q+1}]
    and b marks which positions the new innovation eps_{t+1} enters.

    Parameters
    ----------
    ar_order : int  (p)
    ma_order : int  (q)

    Returns
    -------
    np.ndarray of shape (p+q, 1)

    Examples
    --------
    >>> ARMA_vector_b(2, 2)   # [1, 0, 1, 0]
    >>> ARMA_vector_b(2, 0)   # [1, 0]
    >>> ARMA_vector_b(0, 2)   # [1, 0]
    """
    # AR block: e_p = [1, 0, ..., 0]  (length p), empty if p == 0
    if ar_order == 0:
        e_p = []
    elif ar_order == 1:
        e_p = [1.0]
    else:
        e_p = [1.0] + [0.0] * (ar_order - 1)

    # MA block: e_q = [1, 0, ..., 0]  (length q), empty if q == 0
    if ma_order == 0:
        e_q = []
    elif ma_order == 1:
        e_q = [1.0]
    else:
        e_q = [1.0] + [0.0] * (ma_order - 1)

    return np.array(e_p + e_q, dtype=float).reshape(-1, 1)


# ---------------------------------------------------------------------------
# ARMA_companion_matrix
# ---------------------------------------------------------------------------

def ARMA_companion_matrix(
    phi: np.ndarray | None = None,
    theta: np.ndarray | None = None,
) -> np.ndarray:
    """
    Build the companion (state-transition) matrix A for an ARMA(p, q) model.

    Matrix structure (p+q) × (p+q):
        Row 0   : [phi_1, …, phi_p, theta_1, …, theta_q]
        Rows 1…p-1 : AR shift block
        Row p   : all-zeros boundary row (only when p≥1 AND q≥1)
        Rows p+1…p+q-1 : MA shift block

    Parameters
    ----------
    phi   : array-like of length p, AR parameters (optional / can be None)
    theta : array-like of length q, MA parameters (optional / can be None)

    Returns
    -------
    np.ndarray of shape (p+q, p+q)
    Attributes stored on the array:
        .arOrder  (int)
        .maOrder  (int)

    Examples
    --------
    >>> ARMA_companion_matrix(phi=[0.4])
    >>> ARMA_companion_matrix(theta=[0.4])
    >>> ARMA_companion_matrix(phi=[0.4, 0.3, 0.1])
    >>> ARMA_companion_matrix(phi=[0.4, 0.2], theta=[0.3])
    """
    phi   = np.asarray(phi,   dtype=float) if phi   is not None else np.array([])
    theta = np.asarray(theta, dtype=float) if theta is not None else np.array([])
    p, q = len(phi), len(theta)

    # --- AR shift block (rows 1 … p-1) ---
    L_p: np.ndarray | None = None
    if p > 0:
        # [I_{p-1} | 0_col | 0_{q cols}]   shape: (p-1) x (p+q)
        L_p = np.zeros((p - 1, p + q))
        if p > 1:
            L_p[:, :p - 1] = np.eye(p - 1)

    # --- MA shift block (rows p+1 … p+q-1) ---
    L_q: np.ndarray | None = None
    if q > 0:
        # [0_{p cols} | I_{q-1} | 0_col]   shape: (q-1) x (p+q)
        L_q = np.zeros((q - 1, p + q))
        if q > 1:
            L_q[:, p : p + q - 1] = np.eye(q - 1)

    # --- all-zeros boundary row (only when both p≥1 and q≥1) ---
    zeros: np.ndarray | None = None
    if p >= 1 and q >= 1:
        zeros = np.zeros((1, p + q))

    # --- stack rows ---
    top_row = np.concatenate([phi, theta]).reshape(1, -1)
    blocks  = [top_row]
    if L_p  is not None: blocks.append(L_p)
    if zeros is not None: blocks.append(zeros)
    if L_q  is not None: blocks.append(L_q)

    A = np.vstack(blocks)

    # Store orders as array attributes (mirrors R's attr(A, "arOrder"))
    A.arOrder = p   # type: ignore[attr-defined]
    A.maOrder = q   # type: ignore[attr-defined]
    return A


# ---------------------------------------------------------------------------
# ARMA_expectation
# ---------------------------------------------------------------------------

def ARMA_expectation(
    h: int,
    X0: np.ndarray,
    A: np.ndarray,
    b: np.ndarray,
    intercept: float = 0.0,
) -> np.ndarray:
    """
    Compute the conditional mean of an ARMA process h steps ahead.

    Formula:
        E[y_{t+k} | X0] = [(I - A^k)(I - A)^{-1} c + A^k X0][0]

    where c = [intercept, 0, …, 0]^T.

    Parameters
    ----------
    h         : int, number of steps ahead
    X0        : array of length p+q, current state
    A         : (p+q)×(p+q) companion matrix
    b         : (p+q, 1) or (p+q,) selection vector (unused here, kept for API parity)
    intercept : float

    Returns
    -------
    np.ndarray of shape (h,) with names "t+1", "t+2", …
    """
    X0  = np.asarray(X0, dtype=float).ravel()
    A   = np.asarray(A,  dtype=float)
    pq  = A.shape[0]

    # Intercept vector c = [intercept, 0, …, 0]
    c_ = np.zeros(pq)
    c_[0] = intercept

    # Long-run mean contribution: (I - A)^{-1} c
    I = np.eye(pq)
    try:
        I_A_c_inv = np.linalg.solve(I - A, c_)
    except np.linalg.LinAlgError:
        I_A_c_inv = np.zeros(pq)

    expectations = np.zeros(h)
    A_pow = I.copy()
    for k in range(1, h + 1):
        A_pow = A_pow @ A              # A^k
        mean_k = (I - A_pow) @ I_A_c_inv + A_pow @ X0
        expectations[k - 1] = mean_k[0]

    return pd.Series(expectations, index=[f"t+{k}" for k in range(1, h + 1)])


# ---------------------------------------------------------------------------
# ARMA_variance
# ---------------------------------------------------------------------------

def ARMA_variance(
    h: int,
    A: np.ndarray,
    b: np.ndarray,
    sigma2: float = 1.0,
) -> np.ndarray:
    """
    Compute h-step-ahead forecast variances of an ARMA process.

    var[t+k] = sigma2 * [sum_{j=0}^{k-1} A^j b b^T (A^j)^T]_{0,0}

    Parameters
    ----------
    h      : int, steps ahead
    A      : (p+q)×(p+q) companion matrix
    b      : (p+q, 1) or (p+q,) selection vector
    sigma2 : float, innovation variance (NOT std dev)

    Returns
    -------
    pd.Series of shape (h,) indexed "t+1", …, "t+h"
    """
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float).ravel()
    n = A.shape[0]

    bb     = np.outer(b, b)            # b b^T
    I      = np.eye(n)
    cumsum = I.copy()                  # j=0 term: A^0 bb (A^0)^T = I (scaled by b)
    A_pow  = I.copy()

    variances = np.zeros(h)
    variances[0] = cumsum[0, 0] * sigma2   # t+1

    for i in range(1, h):              # i = 1 … h-1  (corresponds to step t+i+1)
        A_pow  = A @ A_pow             # A^i
        cumsum = cumsum + A_pow @ bb @ A_pow.T
        variances[i] = cumsum[0, 0] * sigma2

    return pd.Series(variances, index=[f"t+{k}" for k in range(1, h + 1)])


# ---------------------------------------------------------------------------
# ARMA_covariance
# ---------------------------------------------------------------------------

def ARMA_covariance(
    h: int,
    k: int,
    A: np.ndarray,
    b: np.ndarray,
    sigma2: float = 1.0,
) -> np.ndarray:
    """
    Compute the forecast covariance Cov(y_{t+h}, y_{t+k} | F_t).

    cov = sigma2 * [sum_{j=0}^{min(h,k)-1} A^{h-1-j} b b^T (A^{k-1-j})^T]

    Parameters
    ----------
    h, k   : int, step indices (h and k steps ahead)
    A      : companion matrix
    b      : selection vector
    sigma2 : float, innovation variance

    Returns
    -------
    np.ndarray, full covariance matrix (use [0, 0] for scalar covariance)
    """
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float).ravel()
    bb = np.outer(b, b)

    n   = A.shape[0]
    cov = np.zeros((n, n))

    for j in range(min(h, k)):
        A_hj = np.linalg.matrix_power(A, h - 1 - j)
        A_kj = np.linalg.matrix_power(A, k - 1 - j)
        cov  = cov + A_hj @ bb @ A_kj.T

    return sigma2 * cov


# ---------------------------------------------------------------------------
# ARMA_next_step
# ---------------------------------------------------------------------------

def ARMA_next_step(
    h: int,
    X0: np.ndarray,
    A: np.ndarray,
    b: np.ndarray,
    intercept: float = 0.0,
    eps: np.ndarray | float = 0.0,
) -> np.ndarray:
    """
    Iterate the ARMA state vector h steps forward.

    X_{t+1} = A X_t + b eps_{t+1},  then X_{t+1}[0] += intercept

    Parameters
    ----------
    h         : int, steps ahead
    X0        : array of length p+q, initial state
    A         : companion matrix
    b         : selection vector (p+q,1) or (p+q,)
    intercept : float
    eps       : scalar 0 (pure forecast) or array of length h (simulation)

    Returns
    -------
    np.ndarray of shape (p+q,), the state vector at time t+h
    """
    X0  = np.asarray(X0,  dtype=float).ravel()
    A   = np.asarray(A,   dtype=float)
    b   = np.asarray(b,   dtype=float).ravel()
    eps = np.asarray(eps, dtype=float).ravel()

    if eps.size == 1 and float(eps) == 0.0:
        eps = np.zeros(h)
    elif len(eps) != h:
        raise ValueError(f"Length of `eps` ({len(eps)}) must equal `h` ({h}).")
    if len(X0) != A.shape[1]:
        raise ValueError(f"Length of `X0` ({len(X0)}) must equal ncol(A) ({A.shape[1]}).")

    x_t = X0.copy()
    for step in range(h):
        x_t    = A @ x_t + b * eps[step]
        x_t[0] += intercept

    return x_t


# ---------------------------------------------------------------------------
# ARMA_filter
# ---------------------------------------------------------------------------

def ARMA_filter(
    x: np.ndarray,
    A: np.ndarray,
    b: np.ndarray,
    intercept: float = 0.0,
) -> np.ndarray:
    """
    One-step-ahead filter: compute ŷ_t for each t.

    Replaces the C backend (.Call("ARMA_filter_c", …)) from the R source.

    State transition:
        X_{t+1} = A X_t + b * eps_t,   eps_t = y_t - ŷ_t
        X_{t+1}[0] += intercept

    Parameters
    ----------
    x         : array of length T, observed time series
    A         : (p+q)×(p+q) companion matrix
    b         : (p+q, 1) or (p+q,) selection vector
    intercept : float

    Returns
    -------
    np.ndarray of shape (T,), one-step-ahead fitted values ŷ_t
    """
    x = np.asarray(x, dtype=float).ravel()
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float).ravel()
    T = len(x)
    n = A.shape[0]

    state  = np.zeros(n)
    y_hat  = np.zeros(T)

    for t in range(T):
        # one-step prediction: ŷ_t = A[0] @ state + intercept
        pred   = float(A[0] @ state) + intercept if n > 0 else intercept
        y_hat[t] = pred
        # update state with realised residual
        eps    = x[t] - pred
        if n > 0:
            state = A @ state + b * eps

    return y_hat


# ---------------------------------------------------------------------------
# ARMA_forecast
# ---------------------------------------------------------------------------

def ARMA_forecast(
    h: int,
    X0: np.ndarray,
    A: np.ndarray,
    b: np.ndarray,
    intercept: float = 0.0,
) -> pd.DataFrame:
    """
    Compute the full h-step state trajectory (pure forecast, eps=0).

    Replaces the C backend (.Call("ARMA_forecast_c", …)).

    Returns
    -------
    pd.DataFrame with columns:
        'step'     : 1 … h
        'forecast' : ŷ_{t+k} (first element of the state vector)
        'weights'  : (only on last row) nested DataFrame of all steps
    """
    X0 = np.asarray(X0, dtype=float).ravel()
    A  = np.asarray(A,  dtype=float)
    b  = np.asarray(b,  dtype=float).ravel()
    n  = A.shape[0]

    steps     = np.arange(1, h + 1)
    forecasts = np.zeros(h)
    state     = X0.copy()

    for step in range(h):
        state        = A @ state          # eps = 0
        state[0]    += intercept
        forecasts[step] = state[0]

    df = pd.DataFrame({"step": steps, "forecast": forecasts})

    # Attach the full trajectory as a nested column on the last row (mirrors R)
    last_row         = df.tail(1).copy()
    last_row["weights"] = [df]
    return last_row.reset_index(drop=True)


# ---------------------------------------------------------------------------
# AR_variance  (closed-form for p ≤ 4, companion matrix otherwise)
# ---------------------------------------------------------------------------

def AR_variance(phi: np.ndarray, sigma2: float = 1.0) -> float:
    """
    Long-run (unconditional) variance of a stationary AR(p) process.

    Closed-form formulas are used for p ≤ 4.
    For p ≥ 5 the companion-matrix approach via ARMA_variance is used.

    Parameters
    ----------
    phi    : array-like of length p, AR coefficients
    sigma2 : float, innovation variance

    Returns
    -------
    float
    """
    phi = np.asarray(phi, dtype=float).ravel()
    p   = len(phi)

    if p == 1:
        return sigma2 / (1 - phi[0] ** 2)

    elif p == 2:
        denom = (1 - phi[1]) * (1 - phi[0]**2 - phi[1]**2) - 2 * phi[0]**2 * phi[1]
        return sigma2 * (1 - phi[1]) / denom

    elif p == 3:
        a, b_, c_ = phi[0], phi[1], phi[2]
        pt1 = (a + b_ * c_) / (1 - b_ - c_**2 - a * c_)
        pt2 = (a + c_) * pt1 + b_
        pt3 = a * pt2 + b_ * pt1 + c_
        pt0 = 1.0 / (1 - a * pt1 - b_ * pt2 - c_ * pt3)
        return pt0 * sigma2

    elif p == 4:
        phi1, phi2, phi3, phi4 = phi
        p1 = (phi1 + phi3) / (1 - phi4)
        p0 = phi2  / (1 - phi4)
        psi1_num = p1 * phi3 + phi1 * p0 * phi4 + phi1 + phi3 * phi4
        psi1_den = 1 - p1 * (phi3 + phi1 * phi4) - phi2 * (1 + phi4) - phi4**2
        psi1 = psi1_num / psi1_den
        psi2 = p1 * psi1 + p0
        psi3 = phi1 * psi2 + phi2 * psi1 + phi4 * psi1 + phi3
        psi4 = phi1 * psi3 + phi2 * psi2 + phi3 * psi1 + phi4
        return sigma2 / (1 - phi1 * psi1 - phi2 * psi2 - phi3 * psi3 - phi4 * psi4)

    else:
        # General case via companion matrix
        A = ARMA_companion_matrix(phi=phi, theta=None)
        b = ARMA_vector_b(ar_order=p, ma_order=0)
        return float(ARMA_variance(h=1, A=A, b=b, sigma2=sigma2).iloc[0])


# ---------------------------------------------------------------------------
# AR_params_to_zeta  (step-down Durbin-Levinson + atanh)
# ---------------------------------------------------------------------------

def AR_params_to_zeta(phi: np.ndarray) -> np.ndarray:
    """
    Map constrained AR coefficients φ to unconstrained ζ = atanh(κ),
    where κ are the partial-autocorrelation (reflection) coefficients.

    Uses the step-down Levinson-Durbin recursion:
        κ_k = φ^(k)[k]
        φ^(k-1)[j] = (φ^(k)[j] + κ_k * φ^(k)[k-j]) / (1 - κ_k²)

    Parameters
    ----------
    phi : array-like of length p

    Returns
    -------
    np.ndarray of length p, unconstrained parameters ζ
    """
    phi = np.asarray(phi, dtype=float).ravel()
    p   = len(phi)
    if p == 0:
        return np.array([])

    phik  = [None] * p
    kappa = np.zeros(p)
    phik[p - 1] = phi.copy()

    for k in range(p, 0, -1):           # k = p, p-1, …, 1  (1-indexed conceptually)
        kappa_k = phik[k - 1][k - 1]    # last element of φ^(k)
        # clip for numerical safety
        if abs(kappa_k) >= 1.0:
            kappa_k = np.sign(kappa_k) * (1.0 - 1e-12)
        kappa[k - 1] = kappa_k

        if k > 1:
            prev     = phik[k - 1][: k - 1]      # φ^(k)[1..k-1]
            rev_prev = prev[::-1]                  # φ^(k)[k-1..1]
            denom    = 1.0 - kappa_k ** 2
            phik[k - 2] = (prev + kappa_k * rev_prev) / denom

    zeta = np.arctanh(kappa)
    return zeta


def ARMA_params_to_zeta(
    phi: np.ndarray | None = None,
    theta: np.ndarray | None = None,
) -> dict:
    """
    Map ARMA(p, q) constrained parameters (φ, θ) to unconstrained (ζ_φ, ζ_θ).

    Parameters
    ----------
    phi   : array-like of length p, AR parameters
    theta : array-like of length q, MA parameters

    Returns
    -------
    dict with keys 'zeta_phi' and 'zeta_theta'
    """
    zeta_phi   = np.array([])
    zeta_theta = np.array([])

    if phi is not None:
        phi = np.asarray(phi, dtype=float).ravel()
        if len(phi) > 0 and phi[0] != 0.0:
            zeta_phi = AR_params_to_zeta(phi)
            zeta_phi = pd.Series(
                zeta_phi, index=[f"zeta_phi_{i+1}" for i in range(len(zeta_phi))]
            )

    if theta is not None:
        theta = np.asarray(theta, dtype=float).ravel()
        if len(theta) > 0 and theta[0] != 0.0:
            zeta_theta = AR_params_to_zeta(theta)
            zeta_theta = pd.Series(
                zeta_theta, index=[f"zeta_theta_{i+1}" for i in range(len(zeta_theta))]
            )

    return {"zeta_phi": zeta_phi, "zeta_theta": zeta_theta}


# ---------------------------------------------------------------------------
# AR_params_to_phi  (tanh + step-up Durbin-Levinson + Jacobian)
# ---------------------------------------------------------------------------

def AR_params_to_phi(zeta: np.ndarray) -> dict:
    """
    Map unconstrained ζ back to constrained AR coefficients φ, and compute
    the Jacobian ∂φ/∂ζ.

    Uses κ = tanh(ζ) and the step-up Levinson-Durbin recursion:
        φ^(k)[j] = φ^(k-1)[j] - κ_k * φ^(k-1)[k-j]   for j = 1 … k-1
        φ^(k)[k] = κ_k

    Parameters
    ----------
    zeta : array-like of length p

    Returns
    -------
    dict with keys:
        'phi'         — np.ndarray, constrained AR coefficients
        'kappa'       — np.ndarray, reflection coefficients
        'J_phi_kappa' — np.ndarray (p×p), ∂φ/∂κ
        'J_phi_zeta'  — np.ndarray (p×p), ∂φ/∂ζ  (chain rule via sech²)
    """
    zeta  = np.asarray(zeta, dtype=float).ravel()
    p     = len(zeta)
    kappa = np.tanh(zeta)

    phis = [None] * p
    G    = [None] * p       # G[k] = (k+1)×(k+1) Jacobian ∂φ^(k+1)/∂κ_{1..k+1}

    # --- k = 1 (0-indexed: k=0) ---
    phis[0] = np.array([kappa[0]])
    G[0]    = np.array([[1.0]])

    # --- k = 2 … p ---
    for k in range(1, p):              # k is 0-indexed; conceptually K = k+1
        phi_prev = phis[k - 1]         # φ^(K-1), length K-1
        K        = k + 1               # 1-indexed order

        phi_new = np.zeros(K)
        for j in range(1, K):          # j = 1 … K-1
            phi_new[j - 1] = phi_prev[j - 1] - kappa[K - 1] * phi_prev[K - j - 1]
        phi_new[K - 1] = kappa[K - 1]
        phis[k] = phi_new

        # Propagate the Jacobian G^(K) from G^(K-1)
        G_prev = G[k - 1]              # (K-1) × (K-1)
        Gk     = np.zeros((K, K))

        # Columns m = 0 … K-2 (propagate from G_prev)
        for m in range(K - 1):
            for j in range(1, K):      # row j-1 in 0-indexed
                Gk[j - 1, m] = G_prev[j - 1, m] - kappa[K - 1] * G_prev[K - j - 1, m]

        # Column m = K-1 (new κ_K)
        for j in range(1, K):
            Gk[j - 1, K - 1] = -phi_prev[K - j - 1]
        Gk[K - 1, K - 1] = 1.0

        G[k] = Gk

    phi          = phis[p - 1]
    J_phi_kappa  = G[p - 1]

    # Chain rule: ∂φ/∂ζ = (∂φ/∂κ) * diag(sech²(ζ))
    sech2       = 1.0 - kappa ** 2    # sech²(ζ) = 1 - tanh²(ζ)
    J_phi_zeta  = J_phi_kappa * sech2[np.newaxis, :]   # broadcast

    return {
        "phi":          phi,
        "kappa":        kappa,
        "J_phi_kappa":  J_phi_kappa,
        "J_phi_zeta":   J_phi_zeta,
    }


def ARMA_params_to_phi(
    zeta_phi: np.ndarray | None = None,
    zeta_theta: np.ndarray | None = None,
) -> dict:
    """
    Map unconstrained ARMA parameters (ζ_φ, ζ_θ) to constrained (φ, θ)
    and compute the block-diagonal Jacobian J = diag(∂φ/∂ζ_φ, ∂θ/∂ζ_θ).

    Parameters
    ----------
    zeta_phi   : array-like of length p
    zeta_theta : array-like of length q

    Returns
    -------
    dict with keys: 'phi', 'theta', 'J'

    Examples
    --------
    >>> ARMA_params_to_phi(zeta_phi=[0.7, 0.465, 0.1], zeta_theta=[0.549, 1.099])
    """
    p = q = 0
    phi = theta = np.array([])
    J_phi_block = J_theta_block = None

    if zeta_phi is not None and len(zeta_phi) > 0:
        zeta_phi = np.asarray(zeta_phi, dtype=float).ravel()
        p        = len(zeta_phi)
        res_phi  = AR_params_to_phi(zeta_phi)
        phi      = pd.Series(res_phi["phi"], index=[f"phi_{i+1}"   for i in range(p)])
        J_phi_block = res_phi["J_phi_zeta"]

    if zeta_theta is not None and len(zeta_theta) > 0:
        zeta_theta  = np.asarray(zeta_theta, dtype=float).ravel()
        q           = len(zeta_theta)
        res_theta   = AR_params_to_phi(zeta_theta)
        theta       = pd.Series(res_theta["phi"], index=[f"theta_{i+1}" for i in range(q)])
        J_theta_block = res_theta["J_phi_zeta"]

    # Block-diagonal Jacobian
    J = np.zeros((p + q, p + q))
    if p > 0 and J_phi_block is not None:
        J[:p, :p] = J_phi_block
    if q > 0 and J_theta_block is not None:
        J[p:p + q, p:p + q] = J_theta_block

    return {"phi": phi, "theta": theta, "J": J}


# ---------------------------------------------------------------------------
# Loss functions (constrained parametrisation)
# ---------------------------------------------------------------------------

def ARMA_loss_logLik(
    params: np.ndarray,
    p: int = 0,
    q: int = 0,
    y: np.ndarray = None,
    per_obs: bool = False,
) -> np.ndarray | float:
    """
    Gaussian log-likelihood in the unconstrained (ζ) parametrisation.

    Parameters
    ----------
    params  : array-like, [ζ_φ₁,…,ζ_φₚ, ζ_θ₁,…,ζ_θ_q]
              Optionally prepended with intercept if params['intercept'] is not NaN.
    p, q    : int, AR/MA orders
    y       : array-like, observed time series
    per_obs : bool — if True, return per-observation log-likelihoods (array);
              if False, return the scalar sum.

    Returns
    -------
    float or np.ndarray
    """
    params = np.asarray(params, dtype=float).ravel()
    y      = np.asarray(y,      dtype=float).ravel()

    # --- extract intercept if provided as a named element ---
    intercept = 0.0
    # (In the pure-numpy API, intercept is the last element when p+q < len(params))

    # --- AR unconstrained parameters ---
    zeta_phi   = params[:p]   if p > 0 else None
    zeta_theta = params[p:p+q] if q > 0 else None

    # --- convert to constrained ---
    res   = ARMA_params_to_phi(zeta_phi, zeta_theta)
    A     = ARMA_companion_matrix(
        phi=res["phi"].values if len(res["phi"]) > 0 else None,
        theta=res["theta"].values if len(res["theta"]) > 0 else None,
    )
    b     = ARMA_vector_b(ar_order=p, ma_order=q)

    # --- filter ---
    y_hat = ARMA_filter(y, A, b, intercept=intercept)

    # --- residuals (exclude warm-up) ---
    burn      = max(p, q) if max(p, q) > 0 else 0
    eps_hat   = (y - y_hat)[burn:]
    sigma_hat = float(np.sqrt(np.mean(eps_hat ** 2)))
    if sigma_hat == 0.0:
        sigma_hat = 1e-10

    z_hat  = eps_hat / sigma_hat
    loglik = -0.5 * np.log(2 * np.pi) - np.log(sigma_hat) - 0.5 * z_hat ** 2

    return loglik if per_obs else float(np.sum(loglik))


def ARMA_loss_CSS(
    params: np.ndarray,
    p: int = 0,
    q: int = 0,
    y: np.ndarray = None,
) -> float:
    """
    Conditional Sum-of-Squares loss in the unconstrained (ζ) parametrisation.

    Parameters
    ----------
    params : array-like, [ζ_φ₁,…,ζ_φₚ, ζ_θ₁,…,ζ_θ_q]
    p, q   : int, AR/MA orders
    y      : array-like, observed time series

    Returns
    -------
    float, SSE (excluding warm-up observations)
    """
    params = np.asarray(params, dtype=float).ravel()
    y      = np.asarray(y,      dtype=float).ravel()

    intercept  = 0.0
    zeta_phi   = params[:p]    if p > 0 else None
    zeta_theta = params[p:p+q] if q > 0 else None

    res   = ARMA_params_to_phi(zeta_phi, zeta_theta)
    A     = ARMA_companion_matrix(
        phi=res["phi"].values   if len(res["phi"])   > 0 else None,
        theta=res["theta"].values if len(res["theta"]) > 0 else None,
    )
    b     = ARMA_vector_b(ar_order=p, ma_order=q)

    y_hat   = ARMA_filter(y, A, b, intercept=intercept)
    burn    = max(p, q) if max(p, q) > 0 else 0
    eps_hat = (y - y_hat)[burn:]

    return float(np.sum(eps_hat ** 2))


# ---------------------------------------------------------------------------
# ARMA_fit
# ---------------------------------------------------------------------------

def ARMA_fit(
    y: np.ndarray,
    ar_order: int = 1,
    ma_order: int = 0,
    method: str = "CSS",
) -> dict:
    """
    Fit an ARMA(p, q) model by optimising over unconstrained parameters.

    Mirrors the R function ARMA_fit().

    Parameters
    ----------
    y        : array-like, time series
    ar_order : int (p)
    ma_order : int (q)
    method   : 'CSS' or 'ML'

    Returns
    -------
    dict with keys:
        'phi'              — pd.Series, constrained AR parameters
        'theta'            — pd.Series, constrained MA parameters
        'J'                — Jacobian matrix ∂(φ,θ)/∂ζ
        'loss'             — scalar, optimised loss value
        'H'                — Hessian of log-likelihood at ζ*
        'vcov'             — dict of variance-covariance matrices
        'std_errors'       — dict with keys 'zeta', 'phi' (classical std errors)
        'std_errors_rob'   — dict with keys 'zeta', 'phi' (sandwich std errors)

    Notes
    -----
    Requires ``numdifftools`` for numerical Hessian and Jacobian:
        pip install numdifftools
    Falls back to a warning if numdifftools is not available.
    """
    y  = np.asarray(y,  dtype=float).ravel()
    pq = ar_order + ma_order
    method = method.upper()
    if method not in ("CSS", "ML"):
        raise ValueError("`method` must be 'CSS' or 'ML'.")

    # --- initial unconstrained parameters (small random values) ---
    rng       = np.random.default_rng(seed=42)
    zeta_init = rng.uniform(-0.1, 0.1, size=pq)
    zeta_init_names = (
        [f"zeta_phi_{i+1}"   for i in range(ar_order)] +
        [f"zeta_theta_{i+1}" for i in range(ma_order)]
    )

    # --- loss closure ---
    if method == "CSS":
        loss_fn = lambda params: ARMA_loss_CSS(params, p=ar_order, q=ma_order, y=y)
    else:
        loss_fn = lambda params: -ARMA_loss_logLik(params, p=ar_order, q=ma_order, y=y)

    # --- optimise ---
    opt    = optimize.minimize(loss_fn, zeta_init, method="Nelder-Mead")
    zeta_opt = opt.x

    # --- convert to constrained parameters ---
    zeta_phi_opt   = zeta_opt[:ar_order]  if ar_order > 0 else None
    zeta_theta_opt = zeta_opt[ar_order:]  if ma_order > 0 else None
    res            = ARMA_params_to_phi(zeta_phi_opt, zeta_theta_opt)

    # --- numerical Hessian and Jacobian of log-likelihood ---
    try:
        import numdifftools as nd

        loglik_fn = lambda params: ARMA_loss_logLik(params, p=ar_order, q=ma_order, y=y)
        loglik_per_obs_fn = lambda params: ARMA_loss_logLik(
            params, p=ar_order, q=ma_order, y=y, per_obs=True
        )

        H  = nd.Hessian(loglik_fn)(zeta_opt)        # (pq × pq) Hessian of sum(loglik)
        S  = nd.Jacobian(loglik_per_obs_fn)(zeta_opt) # (n_obs × pq) per-obs scores

    except ImportError:
        warnings.warn(
            "numdifftools is not installed; Hessian and sandwich SE not computed. "
            "Install with: pip install numdifftools"
        )
        H, S = None, None

    # --- variance-covariance matrices ---
    vcov = vcov_rob = {}
    std_errors     = {"zeta": None, "phi": None}
    std_errors_rob = {"zeta": None, "phi": None}

    if H is not None:
        # Classical: V_star = (-H)^{-1}
        try:
            V_star = np.linalg.solve(-H, np.eye(pq))
        except np.linalg.LinAlgError:
            V_star = np.full((pq, pq), np.nan)

        J        = res["J"]
        V_orig   = J @ V_star @ J.T           # delta method to constrained space

        std_errors["zeta"] = pd.Series(
            np.sqrt(np.abs(np.diag(V_star))), index=zeta_init_names
        )
        phi_theta_names = (
            [f"phi_{i+1}"   for i in range(ar_order)] +
            [f"theta_{i+1}" for i in range(ma_order)]
        )
        std_errors["phi"] = pd.Series(
            np.sqrt(np.abs(np.diag(V_orig))), index=phi_theta_names
        )

        # Sandwich: V_rob = V_star B V_star
        if S is not None:
            B         = S.T @ S                # (pq × pq) outer-product of scores
            V_rob_star = V_star @ B @ V_star
            V_rob_orig = J @ V_rob_star @ J.T

            std_errors_rob["zeta"] = pd.Series(
                np.sqrt(np.abs(np.diag(V_rob_star))), index=zeta_init_names
            )
            std_errors_rob["phi"] = pd.Series(
                np.sqrt(np.abs(np.diag(V_rob_orig))), index=phi_theta_names
            )
            vcov_rob = {"V_star": V_rob_star, "V_orig": V_rob_orig}
        else:
            V_rob_star = V_rob_orig = None

        vcov = {
            "V_star":     V_star,
            "V_orig":     V_orig,
            "V_rob_star": V_rob_star,
            "V_rob_orig": V_rob_orig,
        }

    res["loss"]          = pd.Series({method: opt.fun})
    res["H"]             = H
    res["vcov"]          = vcov
    res["vcov_rob"]      = vcov_rob
    res["std_errors"]    = std_errors
    res["std_errors_rob"]= std_errors_rob

    return res