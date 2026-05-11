"""
Phase 1 calibration pipeline — Appendix A.1 of Romagnoli & Sartini (2025).

Implements the strict six-step sequential estimation pipeline on a single
training window and emits the Tables A1, A2, A3 records.  The function
``calibrate_window`` returns a typed :class:`CalibrationWindowResult` that
holds every intermediate diagnostic; ``calibrate_all_windows`` runs the
loop over the 10 published cutoffs (2013-2022 by default) and aggregates
the per-window results into the three table DataFrames.

Pipeline (one window):
    1. Constrained-OLS clear sky model C_t = δ0 + δ1·H_t + δ2·sin(ωt) + δ3·cos(ωt)
       subject to C_t > R_t for all t.
    2. Data-driven α, β with ε = 1e-3·min_t(1 - R/C); transform Y_t (Eq. A.2).
    3. OLS seasonal mean Ȳ_t = a0 + a1·sin(ωt) + a2·cos(ωt).
    4. Two-stage θ: (a) OLS on (ΔY)^2 for a θ-free variance; (b) Bibby-
       Sørensen martingale estimator on the standardised increments.
    5. Residual ε_t = (Y_t - Ȳ_t) - (Y_{t-1} - Ȳ_{t-1})·exp(-θ); OLS on ε_t^2
       gives (b0, b1, b2).  Proposition A.2 inversion ⇒ continuous-time
       (c0, c1, c2) and (γ0, γ1, γ2).
    6. J(t-1, t, t) and I(t-1, t, t) from Propositions A.1 / Eq. B.5.
       Standardise ε̃_t = ε_t / √J.  Per-calendar-month two-component EM
       mixture (sklearn) sorted by ascending mean (cloudy=lower, sunny=
       higher).  Correction factor √(I/J) applied to the means before the
       analytic mixture moments are evaluated for Table A3.

The pipeline does **not** call ``SolarModel.fit()`` because that
orchestrator carries deferred SoREd / GARCH machinery; instead it composes
the stable Phase-0 building blocks directly.

Notes on convergence and fidelity to the paper
----------------------------------------------
* ``ω = 2π/365`` is enforced everywhere (constant in
  :mod:`solarrpy.radiationModel_internals`).
* The clear-sky constrained fit is implemented by
  :func:`solarrpy.seasonalClearsky.clearsky_optimizer` (Nelder-Mead with a
  large per-violation penalty).  Standard errors for δ are taken from the
  pre-constrained OLS fit and are *not* re-derived under the constraint;
  this matches the existing reference notebook and the plan's accepted
  R-vs-Python tolerance for Table A1.
* The mean-correction factor in Step 6 follows the plan literally:
  ``mu_continuous = mu_discrete * sqrt(I/J)`` (the alternative ``I/√J``
  used in the legacy notebook is preserved in the result diagnostics for
  inspection but is not what is reported in Table A3).
"""

from __future__ import annotations

import calendar
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from scipy.optimize import minimize as _scipy_minimize

from .seasonalClearsky import (
    SeasonalClearsky,
    control_seasonalClearsky,
    clearsky_outliers,
)
from .solarTransform import SolarTransform
from .seasonalModel import SeasonalModel
from .radiationModel_internals import (
    martingale_method_seasonal,
    reparam_seasonal_function,
    integral_sigma2_formula,
    integral_sigma_numeric,
)
from .gaussianMixture import GaussianMixtureModel
from .gaussianMixture_internals import gm_moments

OMEGA = 2.0 * np.pi / 365.0  # paper-mandated angular frequency
DEFAULT_TRAIN_END_YEARS: tuple[int, ...] = tuple(range(2013, 2023))


# ---------------------------------------------------------------------------
# Result containers
# ---------------------------------------------------------------------------


@dataclass
class CalibrationWindowResult:
    """Typed bundle of all intermediate quantities for one training window."""

    train_end_year: int
    n_obs: int

    # --- Step 1: clear-sky δ
    delta: np.ndarray            # length-4 [δ0, δ1, δ2, δ3]
    delta_se: np.ndarray         # length-4 standard errors (pre-constraint OLS)
    clearsky_model: object       # the fitted SeasonalClearsky instance

    # --- Step 2: bound transform
    alpha: float
    beta: float
    epsilon: float
    transform: SolarTransform

    # --- Step 3: seasonal mean of Y
    a: np.ndarray                # length-3 [a0, a1, a2]
    a_se: np.ndarray             # length-3 standard errors
    seasonal_mean_model: SeasonalModel

    # --- Step 4: martingale θ
    theta: float                 # > 0 expected (mean-reverting)
    b_star: np.ndarray           # length-3 stage-4a OLS [b*0, b*1, b*2]
    b_star_se: np.ndarray
    seasonal_b_star_model: SeasonalModel

    # --- Step 5: variance (b → c, γ)
    b: np.ndarray                # length-3 [b0, b1, b2] from ε^2 OLS
    b_se: np.ndarray
    c: np.ndarray                # length-3 continuous-time short-term [c0, c1, c2]
    gamma: np.ndarray            # length-3 long-term [γ0, γ1, γ2]
    seasonal_var_model: SeasonalModel

    # --- Step 6: monthly mixture
    monthly_mixture: pd.DataFrame
    """One row per calendar month, columns:
        Month, Month_abbr, n_obs, p, mu1, sd1, N1, mu0, sd0, N0,
        mean, variance, skewness, kurtosis,
        I_avg, J_avg, mean_corr_paper, mean_corr_legacy
    """

    # --- diagnostics: full enriched training frame
    df_train: pd.DataFrame = field(repr=False)


@dataclass
class CalibrationAllWindowsResult:
    """Aggregated outputs across all training windows."""

    table_a1: pd.DataFrame  # 10 rows, schema matches results/TableA1.csv + train_years/n_obs
    table_a2: pd.DataFrame  # 10 rows: train_end_year, theta, b*, b, c, γ
    table_a3: Dict[int, pd.DataFrame]  # year -> 12-row monthly mixture frame
    windows: List[CalibrationWindowResult] = field(repr=False)


# ---------------------------------------------------------------------------
# Single-window pipeline
# ---------------------------------------------------------------------------


def _step1_clear_sky(
    df_train: pd.DataFrame,
    *,
    target_col: str,
    clearsky_col: str,
    date_col: str,
    lat: float,
    altitude: Optional[float],
    optimiser: str,
) -> tuple[SeasonalClearsky, np.ndarray, np.ndarray]:
    """Constrained OLS clear-sky fit.  Returns (model, delta, delta_se).

    Three ``optimiser`` paths are supported:

    * ``'delta_optimiser'`` (**default** — Table A1 parity): unconstrained
      OLS on CAMS clearsky regressed on ``(1, H_t, \sin, \cos)`` followed by
      the scalar ``δ`` search in :func:`clearsky_delta_optimizer` that scales
      ``δ·C_t`` above ``R_t`` everywhere.  This matches ``vignettespy/
      notebook-appendixA`` and aligns the published ``(δ, β)`` block and
      seasonal ``(a₀, a₁, a₂)`` with the Bologna values in the thesis
      ``paper/TableA1.png`` far more closely than penalty or SLSQP paths
      (verified by mean absolute error over 2013–2022).

    * ``'constrained'``: penalty Nelder-Mead in :func:`clearsky_optimizer`;
      can drift δ away from the CAMS OLS surface — Table A1 discrepancies
      grow materially.

    * ``'slsqp'``: Appendix A.1.1–literal inequality-constrained OLS via
      :func:`scipy.optimize.minimize` (SLSQP), using the same H0 vector as
      ``predict``.  δ shape diverges from the thesis table because the
      binding constraint is much tighter than the legacy scalar scaling.

    Standard errors are always pulled from the *unconstrained* OLS
    ``_ols_bse`` snapshot — under SLSQP the constraint is active on a few
    days only, but a proper constrained variance-covariance matrix would
    require active-set sensitivity analysis which the paper does not use.
    Reporting OLS SEs matches Table A1's convention.
    """

    control = control_seasonalClearsky(
        orders=1, order_H0=1, periods=365,
        include_intercept=True, include_trend=False,
        delta0=1.4, lower=0, upper=3, by=0.001, ntol=0, quiet=True,
    )
    model = SeasonalClearsky(control=control)

    # IMPORTANT — H0 sourcing.  ``SeasonalClearsky.fit`` recomputes H0 via
    # ``Hon(n, lat)`` but ``SeasonalClearsky.predict`` honours any
    # pre-existing ``H0`` column in ``newdata``.  The published Table A1
    # values were generated by the existing ``vignettespy/notebook-appendixA``
    # which calls ``fit(...)`` *without* an explicit H0 (so fit uses
    # Hon-recomputed H0) and then calls ``predict(newdata=df_train)`` (so
    # predict uses the dataset's CSV H0).  These two H0 vectors disagree
    # by ~1% on the Bologna dataset.  Faithfully reproducing the published
    # δ requires preserving that asymmetry, so we *do not* override H0
    # during fit, even though it is technically inconsistent.
    #
    # The SLSQP path is plan-literal and uses a self-consistent H0 (the
    # CSV column when available, else Hon).  Its δ values therefore differ
    # from Table A1 in shape — but it satisfies the inequality
    # ``C_t ≥ R_t`` exactly under the same H0 used by ``predict``.

    # Always run the underlying OLS first so ``_ols_params`` and ``_ols_bse``
    # snapshots are populated.
    underlying = "delta_optimiser" if optimiser == "delta_optimiser" else "constrained"
    model.fit(
        x=df_train[target_col],
        date=df_train[date_col],
        lat=lat,
        clearsky=df_train[clearsky_col],
        alt=altitude,
        optimiser=underlying,
    )

    delta_ols = np.asarray(model._ols_params.values[:4], dtype=float)
    delta_se = np.asarray(model._ols_bse.values[:4], dtype=float)

    if optimiser == "slsqp":
        # Reconstruct the design matrix used by ``predict``: it pulls H0
        # from ``newdata['H0']`` if present, else from ``Hon(n, lat)``.
        if "H0" in df_train.columns:
            H0 = np.asarray(df_train["H0"].values, dtype=float)
        else:
            from .seasonalSolarFunctions import SeasonalSolarFunctions as _SSF
            ssf = _SSF()
            H0 = np.asarray(ssf.Hon(df_train["n"].values.astype(float), lat), dtype=float)
        n = df_train["n"].values.astype(float)
        X = np.column_stack([
            np.ones(len(df_train)),
            H0,
            np.sin(OMEGA * n),
            np.cos(OMEGA * n),
        ])
        y_cams = df_train[clearsky_col].values.astype(float)
        R_t = df_train[target_col].values.astype(float)

        def _obj(p):
            r = y_cams - X @ p
            return float(r @ r)

        def _grad(p):
            return -2.0 * X.T @ (y_cams - X @ p)

        cons = ({
            "type": "ineq",
            "fun": lambda p: X @ p - R_t,
            "jac": lambda p: X,
        },)

        result = _scipy_minimize(
            _obj, x0=delta_ols, jac=_grad, method="SLSQP",
            constraints=cons,
            options={"maxiter": 500, "ftol": 1e-10, "disp": False},
        )
        if not result.success:
            import warnings as _w
            _w.warn(
                f"SLSQP did not converge ({result.message}); using unconstrained OLS δ."
            )
            delta_final = delta_ols
        else:
            delta_final = np.asarray(result.x, dtype=float)

        # Push the SLSQP coefficients back so ``cs_model.predict(...)``
        # returns C_t under the constrained δ.
        param_index = list(model._model.params.index)
        new_params = {k: float(v) for k, v in zip(param_index[:4], delta_final)}
        for k in param_index[4:]:
            new_params[k] = float(model._model.params[k])
        model.update(new_params)
        delta = delta_final
    else:
        delta = np.asarray(model._model.params.values[:4], dtype=float)

    return model, delta, delta_se


def _step2_alpha_beta(
    df_train: pd.DataFrame,
    *,
    Ct_col: str,
    target_col_clean: str,
    epsilon: float = 1e-3,
    clamp_floor_rel: float = 1e-3,
) -> tuple[SolarTransform, dict]:
    """Compute α, β and instantiate a SolarTransform with link='invgumbel'.

    Numerical safeguard
    -------------------
    The paper's data-driven recipe sets
    ``ε = 1e-3 · min_t(1 - R/C)``,
    ``α = min - ε``, ``β = max - min + 2ε``,
    which gives ``Xt_prime = (Xt - α) / β = 1 - ε/β`` at the most-cloudy
    training day.  Under the *constrained* clear-sky fit ``min_t(1 - R/C)``
    is essentially zero (the constraint allows C to touch R from above),
    so ε ≈ 1e-7 and ``Yt = log(-log(1 - 1e-7)) ≈ -16``.

    The pathology only affects the upper bound (the most-cloudy training
    day pinned exactly at ``α + β - ε``), so we widen ``β`` just enough to
    enforce ``max Xt_prime ≤ 1 - clamp_floor_rel``, keeping ``α``
    unchanged.  This bounds the worst-case ``|Yt|`` at
    ``-log(clamp_floor_rel) ≈ 7`` for ``clamp_floor_rel = 1e-3``.  The
    paper's α formula is preserved verbatim; β picks up an additive
    correction of order 1e-3 that vanishes in the limit of a non-binding
    clear-sky constraint.
    """

    transform = SolarTransform(link="invgumbel")
    Xt_raw = np.asarray(transform.X(df_train[target_col_clean], df_train[Ct_col]))
    bounds = transform.fit(Xt_raw, epsilon=epsilon)

    eps_old = float(transform.epsilon)
    beta_old = float(transform.beta)
    eps_floor = beta_old * float(clamp_floor_rel)
    if eps_old < eps_floor:
        # Asymmetric fix: only widen β; α stays positive (BoundTransform.update
        # validates α >= 0 and would silently reject a negative value).
        transform.epsilon = eps_floor
        new_beta = beta_old + (eps_floor - eps_old)
        transform.update(alpha=transform.alpha, beta=new_beta)
        bounds["epsilon"] = eps_floor
        bounds["beta"] = new_beta

    return transform, bounds


def _step3_seasonal_mean(
    df_train: pd.DataFrame,
) -> tuple[SeasonalModel, np.ndarray, np.ndarray]:
    """OLS Ȳ_t = a0 + a1·sin(ωn) + a2·cos(ωn)."""

    model = SeasonalModel(orders=[1], periods=[365])
    model.fit(data=df_train, target_col="Yt", time_col="n", include_intercept=True)
    a = np.asarray(model._model.params.values[:3], dtype=float)
    a_se = np.asarray(model._model.bse.values[:3], dtype=float)
    return model, a, a_se


def _step4_theta(
    Yt: np.ndarray, Yt_bar: np.ndarray
) -> tuple[float, np.ndarray, np.ndarray, SeasonalModel]:
    """Two-stage θ recovery.

    4a — fit (ΔY)^2 = b*0 + b*1·sin(ωn) + b*2·cos(ωn) by OLS.
    4b — Bibby-Sørensen martingale estimator (delegated to
         ``martingale_method_seasonal``, which performs the stage-4a OLS
         internally so we re-derive the (b*) coefficients here for Table A2).
    """
    Yt = np.asarray(Yt)
    Yt_bar = np.asarray(Yt_bar)
    n = len(Yt)

    # Stage 4a — explicit (ΔY)^2 OLS so we can return (b*, se_b*).
    dYt2 = np.full(n, np.nan)
    dYt2[2:] = (Yt[1:-1] - Yt[:-2]) ** 2
    fit_df = pd.DataFrame({"dYt2": dYt2, "n": np.arange(1, n + 1)}).dropna()
    sm_b_star = SeasonalModel(orders=[1], periods=[365])
    sm_b_star.fit(data=fit_df, target_col="dYt2", time_col="n", include_intercept=True)
    b_star = np.asarray(sm_b_star._model.params.values[:3], dtype=float)
    b_star_se = np.asarray(sm_b_star._model.bse.values[:3], dtype=float)

    # Stage 4b — martingale θ.
    theta = float(martingale_method_seasonal(Yt, Yt_bar))
    return theta, b_star, b_star_se, sm_b_star


def _step5_variance_reparam(
    df_train: pd.DataFrame, theta: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, SeasonalModel, dict]:
    """OLS on ε^2 + Proposition A.2 inversion.

    Returns ``(b, b_se, c, gamma, model, reparam_dict)`` where
    ``ε_t = (Y_t - Ȳ_t) - (Y_{t-1} - Ȳ_{t-1})·exp(-θ)``.
    """
    Yt = df_train["Yt"].values
    Yt_bar = df_train["Yt_bar"].values
    Yt_tilde = Yt - Yt_bar  # length N

    # ε_t requires Y_{t-1}; first entry is NaN.
    shifted = np.empty_like(Yt_tilde)
    shifted[0] = np.nan
    shifted[1:] = Yt_tilde[:-1]
    eps = Yt_tilde - shifted * np.exp(-theta)
    df_train["eps"] = eps

    fit_df = df_train.dropna(subset=["eps"]).copy()
    fit_df["eps2"] = fit_df["eps"] ** 2

    sm_var = SeasonalModel(orders=[1], periods=[365])
    sm_var.fit(data=fit_df, target_col="eps2", time_col="n", include_intercept=True)

    b = np.asarray(sm_var._model.params.values[:3], dtype=float)
    b_se = np.asarray(sm_var._model.bse.values[:3], dtype=float)
    reparam = reparam_seasonal_function([b[0], b[1], b[2]], theta, omega=OMEGA)
    c = np.asarray(reparam["c_"], dtype=float)
    gamma = np.asarray(reparam["gamma"], dtype=float)
    return b, b_se, c, gamma, sm_var, reparam


def _step6_monthly_mixture(
    df_train: pd.DataFrame,
    theta: float,
    reparam: dict,
    *,
    em_max_iter: int = 5000,
    em_max_restarts: int = 50,
) -> pd.DataFrame:
    """Per-calendar-month EM with √(I/J) mean correction.

    Adds ``J``, ``I``, ``eps_tilde``, ``corr_paper`` (= √(I/J)), and
    ``corr_legacy`` (= I/√J) columns to ``df_train`` *in place* so the
    caller can inspect the day-level diagnostics.
    """

    days = df_train["n"].values.astype(float)

    # J(t-1, t, t) — long-term variance integral, closed form, uses γ.
    calc_J = integral_sigma2_formula(theta, reparam["gamma"], omega=OMEGA)
    J_vals = calc_J(days - 1.0, days, days)

    # I(t-1, t, t) — short-term standard-deviation integral, numerical, uses c.
    calc_I = integral_sigma_numeric(theta, reparam["c_"], omega=OMEGA)
    I_vals = calc_I(days - 1.0, days, days)

    df_train["J"] = J_vals
    df_train["I"] = I_vals

    # Standardise ε.  J should be > 0 by construction; clamp to defend
    # against unlikely numerical underflow at the seasonal trough.
    safe_J = np.maximum(J_vals, 1e-12)
    df_train["eps_tilde"] = df_train["eps"] / np.sqrt(safe_J)

    # Two correction factors retained for diagnostics; the published Table
    # A3 uses corr_paper per the user's "literal interpretation" choice.
    df_train["corr_paper"] = np.sqrt(I_vals / safe_J)
    df_train["corr_legacy"] = I_vals / np.sqrt(safe_J)

    rows: list[dict] = []
    for month in range(1, 13):
        mask = df_train["Month"] == month
        sub = df_train.loc[mask].dropna(subset=["eps_tilde"])
        if len(sub) < 2:
            continue

        x = sub["eps_tilde"].values
        I_avg = float(sub["I"].mean())
        J_avg = float(sub["J"].mean())
        corr_paper = float(np.sqrt(I_avg / max(J_avg, 1e-12)))
        corr_legacy = float(I_avg / np.sqrt(max(J_avg, 1e-12)))

        gmm = GaussianMixtureModel(
            components=2, maxit=em_max_iter, maxrestarts=em_max_restarts
        )
        gmm.fit(x)

        mu_raw = np.asarray(gmm.means, dtype=float)
        sd = np.asarray(gmm.sd, dtype=float)
        p = np.asarray(gmm.p, dtype=float)

        # Components are sorted ascending by mean inside gm_fit
        # (Phase-0 confirmed).  cloudy = state 1 = index 0; sunny =
        # state 0 = index 1.
        classifications = gmm.fitted["classification"].values
        N1_cloudy = int((classifications == 1).sum())
        N0_sunny = int((classifications == 2).sum())

        mu_corrected = mu_raw * corr_paper
        moms = gm_moments(mu_corrected, sd, p).iloc[0]

        rows.append({
            "Month": month,
            "Month_abbr": calendar.month_abbr[month],
            "n_obs": int(len(x)),
            "p": float(p[0]),
            "mu1": float(mu_corrected[0]),
            "sd1": float(sd[0]),
            "N1": N1_cloudy,
            "mu0": float(mu_corrected[1]),
            "sd0": float(sd[1]),
            "N0": N0_sunny,
            "mean": float(moms["mean"]),
            "variance": float(moms["variance"]),
            "skewness": float(moms["skewness"]),
            "kurtosis": float(moms["kurtosis"]),
            "I_avg": I_avg,
            "J_avg": J_avg,
            "mean_corr_paper": corr_paper,
            "mean_corr_legacy": corr_legacy,
        })

    return pd.DataFrame(rows)


def calibrate_window(
    df: pd.DataFrame,
    train_end_year: int,
    *,
    coords: Dict[str, float],
    target_col: str = "GHI",
    clearsky_col: str = "clearsky",
    date_col: str = "date",
    optimiser: str = "delta_optimiser",
) -> CalibrationWindowResult:
    """Run the six-step Appendix A.1 pipeline on training data ≤ ``train_end_year``.

    Parameters
    ----------
    df : pd.DataFrame
        Daily records.  Must contain columns ``date``, ``n`` (1-based day
        of year), ``Month``, ``GHI``, ``clearsky``.
    train_end_year : int
        Inclusive end-of-year for the training window (2013..2022 for the
        published reproduction).
    coords : dict
        ``{'lat': ..., 'lon': ..., 'alt': ...}``.  Latitude is mandatory.
    optimiser : {'delta_optimiser', 'slsqp', 'constrained'}
        Clear-sky fitting strategy.

        * ``'delta_optimiser'`` (default, reference-matching): unconstrained
          OLS-vs-CAMS followed by a scalar shift δ that lifts ``δ·C_t``
          above ``R_t`` everywhere.  Reproduces the published Table A1
          values (``results/TableA1.csv``) on Bologna and matches the thesis
          figure ``paper/TableA1.png`` most closely overall—especially for
          ``δ₀ … δ₃`` and ``β``—compared with ``constrained`` (penalty) or
          ``slsqp``.  Uses the legacy notebook fitting path.
        * ``'slsqp'`` (plan-literal): true constrained OLS via
          :func:`scipy.optimize.minimize` (SLSQP) with inequality
          constraints ``X·δ ≥ R_t``, warm-started from the unconstrained
          OLS.  Self-consistent H0 (uses the CSV column when present),
          so ``predict(...)`` exactly satisfies the constraint.  δ shape
          differs from Table A1 because the constraint binds tightly.
        * ``'constrained'`` (legacy): penalty Nelder-Mead inside
          :func:`clearsky_optimizer`; kept for back-compat only.
    """

    # --- 0. Training slice
    df = df.copy()
    df[date_col] = pd.to_datetime(df[date_col])
    mask = df[date_col].dt.year <= train_end_year
    df_train = df.loc[mask].reset_index(drop=True).copy()
    n_obs = len(df_train)

    # --- 1. Clear sky
    cs_model, delta, delta_se = _step1_clear_sky(
        df_train,
        target_col=target_col,
        clearsky_col=clearsky_col,
        date_col=date_col,
        lat=float(coords["lat"]),
        altitude=None,
        optimiser=optimiser,
    )
    columns_to_drop = [col for col in df_train.columns if 'H0' in col]
    df_train["Ct"] = cs_model.predict(newdata=df_train.drop(columns=columns_to_drop))

    # Mandatory scrubbing of any residual outliers (CAMS occasionally
    # records GHI marginally above the constrained envelope; impute them).
    scrubbed = clearsky_outliers(
        df_train[target_col].values,
        df_train["Ct"].values,
        df_train[date_col],
        quiet=True,
    )
    df_train["GHI_clean"] = scrubbed["x"]

    # --- 2. Alpha / beta + Y_t
    transform, bounds = _step2_alpha_beta(
        df_train,
        Ct_col="Ct",
        target_col_clean="GHI_clean",
        epsilon=1e-3,
    )
    alpha = float(bounds["alpha"])
    beta = float(bounds["beta"])
    epsilon = float(bounds["epsilon"])
    df_train["Yt"] = transform.RY(df_train["GHI_clean"], df_train["Ct"])

    # --- 3. Seasonal mean of Y
    sm_yt, a, a_se = _step3_seasonal_mean(df_train)
    df_train["Yt_bar"] = sm_yt.predict(data=df_train, time_col="n")

    # --- 4. θ via martingale
    theta, b_star, b_star_se, sm_b_star = _step4_theta(
        df_train["Yt"].values, df_train["Yt_bar"].values
    )

    # --- 5. variance OLS + reparam
    b, b_se, c, gamma, sm_var, reparam = _step5_variance_reparam(df_train, theta)

    # --- 6. monthly EM mixture
    monthly_mixture = _step6_monthly_mixture(df_train, theta, reparam)

    return CalibrationWindowResult(
        train_end_year=train_end_year,
        n_obs=n_obs,
        delta=delta,
        delta_se=delta_se,
        clearsky_model=cs_model,
        alpha=alpha,
        beta=beta,
        epsilon=epsilon,
        transform=transform,
        a=a,
        a_se=a_se,
        seasonal_mean_model=sm_yt,
        theta=theta,
        b_star=b_star,
        b_star_se=b_star_se,
        seasonal_b_star_model=sm_b_star,
        b=b,
        b_se=b_se,
        c=c,
        gamma=gamma,
        seasonal_var_model=sm_var,
        monthly_mixture=monthly_mixture,
        df_train=df_train,
    )


# ---------------------------------------------------------------------------
# All-windows loop and table assembly
# ---------------------------------------------------------------------------


def _table_a1_row(res: CalibrationWindowResult) -> dict:
    return {
        "train_years": f"2005-{res.train_end_year}",
        "train_end_year": res.train_end_year,
        "n_obs": res.n_obs,
        "alpha": res.alpha,
        "beta": res.beta,
        "delta0": res.delta[0],
        "delta1": res.delta[1],
        "delta2": res.delta[2],
        "delta3": res.delta[3],
        "delta0_err": res.delta_se[0],
        "delta1_err": res.delta_se[1],
        "delta2_err": res.delta_se[2],
        "delta3_err": res.delta_se[3],
        "a0": res.a[0],
        "a1": res.a[1],
        "a2": res.a[2],
        "a0_err": res.a_se[0],
        "a1_err": res.a_se[1],
        "a2_err": res.a_se[2],
    }


def _table_a2_row(res: CalibrationWindowResult) -> dict:
    return {
        "train_end_year": res.train_end_year,
        "theta": res.theta,
        "b_star_0": res.b_star[0],
        "b_star_1": res.b_star[1],
        "b_star_2": res.b_star[2],
        "b_star_0_err": res.b_star_se[0],
        "b_star_1_err": res.b_star_se[1],
        "b_star_2_err": res.b_star_se[2],
        "b0": res.b[0],
        "b1": res.b[1],
        "b2": res.b[2],
        "b0_err": res.b_se[0],
        "b1_err": res.b_se[1],
        "b2_err": res.b_se[2],
        "c0": res.c[0],
        "c1": res.c[1],
        "c2": res.c[2],
        "gamma0": res.gamma[0],
        "gamma1": res.gamma[1],
        "gamma2": res.gamma[2],
    }


def calibrate_all_windows(
    df: pd.DataFrame,
    *,
    coords: Dict[str, float],
    end_years: Sequence[int] = DEFAULT_TRAIN_END_YEARS,
    target_col: str = "GHI",
    clearsky_col: str = "clearsky",
    date_col: str = "date",
    optimiser: str = "delta_optimiser",
    progress: bool = True,
) -> CalibrationAllWindowsResult:
    """Run :func:`calibrate_window` for each end year and assemble the tables."""

    a1_rows: list[dict] = []
    a2_rows: list[dict] = []
    a3_tables: Dict[int, pd.DataFrame] = {}
    windows: list[CalibrationWindowResult] = []

    for y in end_years:
        if progress:
            print(f"[calibrate_all_windows] training period 2005-{y} ...", flush=True)
        res = calibrate_window(
            df, y,
            coords=coords,
            target_col=target_col,
            clearsky_col=clearsky_col,
            date_col=date_col,
            optimiser=optimiser,
        )
        a1_rows.append(_table_a1_row(res))
        a2_rows.append(_table_a2_row(res))
        a3_tables[y] = res.monthly_mixture.copy()
        windows.append(res)

    table_a1 = pd.DataFrame(a1_rows)
    table_a2 = pd.DataFrame(a2_rows)

    return CalibrationAllWindowsResult(
        table_a1=table_a1,
        table_a2=table_a2,
        table_a3=a3_tables,
        windows=windows,
    )


__all__ = [
    "OMEGA",
    "DEFAULT_TRAIN_END_YEARS",
    "CalibrationWindowResult",
    "CalibrationAllWindowsResult",
    "calibrate_window",
    "calibrate_all_windows",
]
