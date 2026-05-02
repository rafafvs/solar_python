"""
Solar model QMLE estimation — Python translation.

Functions:
    solar_model_params_to_zeta   — constrained  → unconstrained parameter vector
    solar_model_params_to_phi    — unconstrained → constrained + Jacobian
    solar_model_quasi_loglik     — quasi log-likelihood (delegates to C extension)
    solar_model_QMLE             — full QMLE estimation loop
"""

from __future__ import annotations

import warnings
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.optimize import approx_fprime


# ---------------------------------------------------------------------------
# solar_model_params_to_zeta
# ---------------------------------------------------------------------------

def solar_model_params_to_zeta(model) -> dict:
    """
    Map constrained model parameters to an unconstrained vector ζ
    suitable for gradient-based optimisation.

    Steps
    -----
    1. Seasonal mean coefficients (passed through as-is).
    2. ARMA (φ, θ) → unconstrained via ``ARMA_params_to_zeta``.
    3. Seasonal variance coefficients b → unconstrained via
       ``seasonal_model_params_to_zeta``.
    4. GARCH (ω, α, β) → unconstrained via ``sGARCH_params_to_zeta``
       (only when ``model.spec["garch_variance"]`` is True).

    Returns
    -------
    dict with keys:
        params      — np.ndarray, unconstrained parameter vector
        orig_names  — list[str], original constrained parameter names
        GARCH_order — tuple(int, int), (archOrder, garchOrder)
    """
    params = np.array([], dtype=float)
    params_names: list[str] = []

    # 1) Seasonal mean
    a = np.asarray(list(model.coefficients["seasonal_model_Yt"].values()), dtype=float)
    params = np.concatenate([params, a])
    params_names += list(model.coefficients["seasonal_model_Yt"].keys())

    # 2) ARMA
    phi   = model.ARMA.phi
    theta = model.ARMA.theta
    if phi is None or len(phi) == 0:
        phi = pd.Series({"phi_1": 0.0})
    if theta is None or len(theta) == 0:
        theta = pd.Series({"theta_1": 0.0})

    zeta_arma_dict = ARMA_params_to_zeta(phi.values, theta.values)
    zeta_flat      = np.concatenate([
        np.atleast_1d(zeta_arma_dict["zeta_phi"]),
        np.atleast_1d(zeta_arma_dict["zeta_theta"]),
    ])
    params = np.concatenate([params, zeta_flat])
    params_names += list(phi.index) + list(theta.index)

    # 3) Seasonal variance
    b      = np.asarray(list(model.coefficients["seasonal_variance"].values()), dtype=float)
    b_star = seasonal_model_params_to_zeta(b)
    b_star_names = [f"b_star_{i+1}" for i in range(len(b_star))]
    params = np.concatenate([params, b_star])
    params_names += list(model.coefficients["seasonal_variance"].keys())

    # 4) GARCH (optional)
    if model.spec.get("garch_variance", False):
        arch_order  = model.GARCH.arch_order
        garch_order = model.GARCH.garch_order
        if arch_order != 0 and garch_order != 0:
            coefs = np.concatenate([[model.GARCH.omega],
                                    np.atleast_1d(model.GARCH.alpha),
                                    np.atleast_1d(model.GARCH.beta)])
        elif arch_order != 0:
            coefs = np.concatenate([[model.GARCH.omega],
                                    np.atleast_1d(model.GARCH.alpha)])
        else:
            coefs = np.array([model.GARCH.omega])

        coefs_star = sGARCH_params_to_zeta(coefs, arch_order, garch_order)
        coefs_star = coefs_star[:-1]    # drop last element (mirror R's head(..., -1))
        params = np.concatenate([params, coefs_star])
        params_names += (
            ["omega"]
            + [f"alpha{i+1}" for i in range(arch_order)]
            + [f"beta{i+1}"  for i in range(garch_order)]
        )

    return {
        "params":      params,
        "orig_names":  params_names,
        "GARCH_order": (
            (model.GARCH.arch_order, model.GARCH.garch_order)
            if model.spec.get("garch_variance", False)
            else (0, 0)
        ),
    }


# ---------------------------------------------------------------------------
# solar_model_params_to_phi
# ---------------------------------------------------------------------------

def solar_model_params_to_phi(
    params: np.ndarray,
    orig_names: list[str],
    GARCH_order: tuple[int, int],
) -> dict:
    """
    Map an unconstrained vector ζ back to constrained model parameters,
    and build the delta-method Jacobian ∂θ/∂ζ.

    Parameters
    ----------
    params      : np.ndarray, unconstrained parameter vector
    orig_names  : list[str], constrained parameter names (rows of J)
    GARCH_order : tuple(int, int), (archOrder, garchOrder)

    Returns
    -------
    dict with keys:
        params      — dict of constrained sub-vectors {a, phi, theta, b, omega, alpha, beta}
        theta_star  — np.ndarray, the original unconstrained vector (ζ)
        theta       — np.ndarray, the constrained vector (θ)
        J           — np.ndarray, Jacobian matrix ∂θ/∂ζ  shape (d+1, d)
    """
    params       = np.asarray(params, dtype=float)
    params_names = np.asarray(
        [f"param_{i}" for i in range(len(params))]
        if not hasattr(params, "index") else list(params.index)
    )
    # If params was passed as a plain ndarray, build placeholder names
    # (the caller should pass a pd.Series or annotated array for correct naming)
    d = len(params)
    row_names = list(orig_names)
    col_names = [f"param_{i}" for i in range(d)]   # overridden below

    # Build Jacobian scaffold (d+1 × d)
    J = np.eye(d + 1, d)

    coefs = np.array([], dtype=float)

    # Helper: boolean mask over params by name prefix
    def starts(prefix: str) -> np.ndarray:
        return np.array([str(n).startswith(prefix) for n in params_names], dtype=bool)

    def regex_match(pattern: str) -> np.ndarray:
        import re
        return np.array([bool(re.search(pattern, str(n))) for n in params_names], dtype=bool)

    # ---- 1) Seasonal mean ----
    idx_a = starts("a_")
    a     = params[idx_a]
    coefs = np.concatenate([coefs, a])

    # ---- 2) ARMA ----
    idx_zeta_phi   = starts("zeta_phi")
    idx_zeta_theta = starts("zeta_theta")
    zeta_phi   = params[idx_zeta_phi]
    zeta_theta = params[idx_zeta_theta]

    J_zeta = ARMA_params_to_phi(zeta_phi, zeta_theta)
    phi    = J_zeta["phi"]
    theta  = J_zeta["theta"]
    coefs  = np.concatenate([coefs, phi, theta])

    if len(phi)   == 0: phi   = np.array([0.0])
    if len(theta) == 0: theta = np.array([0.0])

    # Insert ARMA block into Jacobian
    idx_row_arma = np.array([bool(__import__('re').search(r"phi|theta", n)) for n in row_names])
    idx_col_arma = idx_zeta_phi | idx_zeta_theta
    J[np.ix_(idx_row_arma, idx_col_arma)] = J_zeta["J"]

    # ---- 3) Seasonal variance ----
    idx_b_star = starts("b_star_")
    b_star     = params[idx_b_star]
    b          = seasonal_model_params_to_phi(b_star)
    b          = np.asarray(b, dtype=float)
    # Standard names: c_0, c_sin_1_365, c_cos_1_365
    coefs = np.concatenate([coefs, b])

    idx_row_b = np.array([str(n).startswith("c_") for n in row_names])
    J[np.ix_(idx_row_b, idx_b_star)] = seasonal_model_params_to_zeta_jacobian(b_star)

    # ---- 4) GARCH ----
    import re as _re
    idx_garch_star = np.array([bool(_re.search(r"eta0|kappa[0-9]", str(n))) for n in params_names])
    omega = np.array([1.0])
    alpha = np.array([0.0])
    beta  = np.array([0.0])

    if idx_garch_star.any():
        coef_star = np.concatenate([params[idx_garch_star], [0.0]])
        par       = sGARCH_params_to_phi(coef_star, GARCH_order[0], GARCH_order[1])
        par       = np.asarray(par, dtype=float)
        coefs     = np.concatenate([coefs, par])

        par_names = (
            ["omega"]
            + [f"alpha{i+1}" for i in range(GARCH_order[0])]
            + [f"beta{i+1}"  for i in range(GARCH_order[1])]
        )
        omega = par[[n.startswith("omega") for n in par_names]]
        alpha = par[[n.startswith("alpha") for n in par_names]]
        beta  = par[[n.startswith("beta")  for n in par_names]]

        J_garch = sGARCH_params_to_zeta_jacobian(coef_star, GARCH_order[0], GARCH_order[1])
        idx_row_g = np.array([bool(_re.search(r"alpha|beta|omega", n)) for n in row_names])
        idx_col_g = np.array([bool(_re.search(r"eta0|kappa[0-9]", str(n))) for n in params_names])
        J[np.ix_(idx_row_g, idx_col_g)] = J_garch

    return {
        "params": {
            "a":     a,
            "phi":   phi,
            "theta": theta,
            "b":     b,
            "omega": omega,
            "alpha": alpha,
            "beta":  beta,
        },
        "theta_star": params,
        "theta":      coefs,
        "J":          J,
    }


# ---------------------------------------------------------------------------
# solar_model_quasi_loglik
# ---------------------------------------------------------------------------

def solar_model_quasi_loglik(
    params: np.ndarray,
    Yt: np.ndarray,
    w: np.ndarray,
    t: np.ndarray,
    orig_names: list[str],
    GARCH_order: tuple[int, int],
    neg_loglik: bool = False,
    per_obs: bool = False,
) -> float | np.ndarray:
    """
    Evaluate the quasi log-likelihood at unconstrained parameters ζ.

    This function mirrors R's ``solarModel_quasi_loglik``, which calls a
    compiled C routine. In Python the equivalent is a Cython/C extension
    ``solarrpy._cext.solar_model_quasi_loglik_c`` (you must build this
    separately), or you can substitute a pure-Python implementation.

    Parameters
    ----------
    params      : np.ndarray, unconstrained parameter vector ζ
    Yt          : np.ndarray, transformed solar radiation time series
    w           : np.ndarray, observation weights
    t           : np.ndarray, time index (day-of-year)
    orig_names  : list[str], constrained parameter names
    GARCH_order : tuple(int, int)
    neg_loglik  : bool, return −ℓ when True
    per_obs     : bool, return per-observation contributions when True

    Returns
    -------
    float or np.ndarray
    """
    conv = solar_model_params_to_phi(params, orig_names, GARCH_order)
    par  = conv["params"]

    # Try the compiled C extension first (solarrpy._cext must be built)
    try:
        from solarrpy._cext import solar_model_quasi_loglik_c
        return solar_model_quasi_loglik_c(
            Yt.astype(float),
            w.astype(float),
            t.astype(float),
            par["a"].astype(float),
            par["phi"][par["phi"] != 0].astype(float),
            par["theta"][par["theta"] != 0].astype(float),
            par["b"].astype(float),
            float(par["omega"]),
            par["alpha"][par["alpha"] != 0].astype(float),
            par["beta"][par["beta"] != 0].astype(float),
            bool(neg_loglik),
            bool(per_obs),
        )
    except ImportError:
        raise NotImplementedError(
            "solar_model_quasi_loglik requires the compiled C extension "
            "'solarrpy._cext'. Build it with 'python setup.py build_ext --inplace', "
            "or replace this function body with a pure-Python implementation."
        )


# ---------------------------------------------------------------------------
# solar_model_QMLE
# ---------------------------------------------------------------------------

def solar_model_QMLE(
    model,
    max_restarts: int = 1,
    seed: int = 1,
    quiet: bool = True,
):
    """
    Quasi-Maximum Likelihood Estimation (QMLE) for the solar model.

    Steps
    -----
    1. Extract training data and initialise unconstrained parameters.
    2. Minimise −ℓ(ζ) via ``scipy.optimize.minimize`` (L-BFGS-B).
    3. Optionally restart from random perturbations of ζ.
    4. Recover constrained parameters and compute the sandwich
       variance–covariance matrix via numerical Hessian and score.
    5. Clone the model, update all sub-models, re-filter, re-fit
       the mixture, and update moments and log-likelihood.

    Parameters
    ----------
    model        : SolarModel
    max_restarts : int, number of random restarts (1 = no restarts)
    seed         : int, random seed for restarts
    quiet        : bool, suppress progress printing

    Returns
    -------
    SolarModel — updated clone with QMLE estimates and robust std. errors
    """
    # ---- Training data ----
    data = model.data[model.data["isTrain"]].copy()
    Yt = data["Yt"].values
    t  = data["n"].values
    w  = data["weights"].values

    # ---- Initial unconstrained parameters ----
    init      = solar_model_params_to_zeta(model)
    zeta0     = init["params"]
    orig_names  = init["orig_names"]
    GARCH_order = init["GARCH_order"]

    # ---- Objective (−ℓ) ----
    def neg_loglik_fn(zeta: np.ndarray) -> float:
        return float(solar_model_quasi_loglik(
            zeta, Yt=Yt, w=w, t=t,
            orig_names=orig_names, GARCH_order=GARCH_order,
            neg_loglik=True, per_obs=False,
        ))

    init_loglik = neg_loglik_fn(zeta0)

    opt = minimize(neg_loglik_fn, zeta0, method="L-BFGS-B")
    if not quiet:
        print(f"Log-lik improved by: {init_loglik - opt.fun:.6f}")

    best_zeta    = opt.x
    best_negloglik = opt.fun

    # ---- Random restarts ----
    if max_restarts > 1:
        rng = np.random.default_rng(seed)
        for restart in range(1, max_restarts):
            if not quiet:
                print(f"Restarting: {restart}/{max_restarts - 1}")
            rand_zeta = zeta0 * rng.uniform(size=len(zeta0))
            opt_r = minimize(neg_loglik_fn, rand_zeta, method="L-BFGS-B")
            if not quiet:
                print(f"  New −ℓ: {opt_r.fun:.6f}  Old: {best_negloglik:.6f}")
            if opt_r.fun < best_negloglik:
                if not quiet:
                    print(f"  Log-lik improved by: {abs(opt_r.fun) - abs(best_negloglik):.6f}")
                best_negloglik = opt_r.fun
                best_zeta      = opt_r.x

    # ---- Constrained parameters + Jacobian ----
    par_qmle = solar_model_params_to_phi(best_zeta, orig_names, GARCH_order)
    J_qmle   = par_qmle["J"]           # shape (d+1, d)
    theta_qml = par_qmle["theta"]      # constrained parameter vector

    # ---- Numerical Hessian (of +ℓ, i.e. −neg_loglik) ----
    def loglik_fn(zeta: np.ndarray) -> float:
        return -neg_loglik_fn(zeta)

    def loglik_per_obs(zeta: np.ndarray) -> np.ndarray:
        return np.atleast_1d(solar_model_quasi_loglik(
            zeta, Yt=Yt, w=w, t=t,
            orig_names=orig_names, GARCH_order=GARCH_order,
            neg_loglik=False, per_obs=True,
        ))

    p    = len(best_zeta)
    eps  = np.sqrt(np.finfo(float).eps)

    # Hessian of ℓ (negative for information matrix I = −H)
    H = _numerical_hessian(loglik_fn, best_zeta, eps=eps)

    # Score matrix S[i, :] = ∇_ζ ℓ_i(ζ*)  via finite differences
    S = _numerical_jacobian(loglik_per_obs, best_zeta, eps=eps)  # shape (n, p)

    # Outer-product meat: B = Σ_i s_i s_i'
    B = S.T @ S                          # shape (p, p)

    # Sandwich variance
    I        = -H
    V_star   = np.linalg.solve(I, np.eye(p))       # I^{-1}
    V_rob_star = V_star @ B @ V_star                # sandwich in ζ-space
    V_rob_orig = J_qmle @ V_rob_star @ J_qmle.T    # delta-method to θ-space

    # Robust standard errors (drop the extra row from d+1 × d+1)
    std_errors = np.sqrt(np.diag(V_rob_orig)[: len(theta_qml)])

    # ---- Update model ----
    model_upd = model.clone(deep=True)

    theta_dict = dict(zip(orig_names, theta_qml))
    # Exclude omega from the update (mirrors R's `theta_qml[!(names == "omega")]`)
    theta_no_omega = {k: v for k, v in theta_dict.items() if k != "omega"}
    model_upd.update(theta_no_omega)

    se_dict = dict(zip(orig_names, std_errors))
    for sub in ("seasonal_model_Yt", "ARMA", "seasonal_variance", "GARCH"):
        try:
            getattr(model_upd, sub).update_std_errors(se_dict)
        except Exception:
            pass

    model_upd.filter()
    model_upd.fit_NM_model()
    model_upd.update_moments()
    model_upd.update_logLik()

    return model_upd


# ---------------------------------------------------------------------------
# Numerical differentiation helpers
# ---------------------------------------------------------------------------

def _numerical_hessian(f, x0: np.ndarray, eps: float = 1e-5) -> np.ndarray:
    """
    Central-difference Hessian of scalar-valued f at x0.
    Mirrors ``numDeriv::hessian``.
    """
    n = len(x0)
    H = np.zeros((n, n))
    f0 = f(x0)
    for i in range(n):
        for j in range(i, n):
            ei = np.zeros(n); ei[i] = eps
            ej = np.zeros(n); ej[j] = eps
            if i == j:
                H[i, i] = (f(x0 + 2*ei) - 2*f(x0 + ei) + f0) / eps**2
            else:
                H[i, j] = H[j, i] = (
                    f(x0 + ei + ej) - f(x0 + ei - ej)
                    - f(x0 - ei + ej) + f(x0 - ei - ej)
                ) / (4 * eps**2)
    return H


def _numerical_jacobian(f, x0: np.ndarray, eps: float = 1e-5) -> np.ndarray:
    """
    Forward-difference Jacobian of vector-valued f at x0.
    Mirrors ``numDeriv::jacobian``. Returns shape (m, n).
    """
    f0 = np.atleast_1d(f(x0))
    n  = len(x0)
    m  = len(f0)
    J  = np.zeros((m, n))
    for j in range(n):
        dx      = np.zeros(n); dx[j] = eps
        J[:, j] = (np.atleast_1d(f(x0 + dx)) - f0) / eps
    return J