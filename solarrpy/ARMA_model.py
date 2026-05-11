"""
ARMA(p, q) Model — Python translation of the R6 class ARMA_modelR6.

Dependencies: numpy, pandas, statsmodels
Version: 1.0.2

Numerical-parity note (R vs Python)
-----------------------------------
The R reference uses ``stats::arima(..., method = "CSS")`` (conditional sum of
squares). This Python port delegates fitting to
``statsmodels.tsa.arima.model.ARIMA``, which defaults to **exact state-space
maximum likelihood**. The two estimators produce coefficient estimates that
differ on the order of ``1e-4`` for a ~6000-day Bologna training window and
yield slightly different ``sigma2``. This drift is documented and accepted
because the SoRad calibration pipeline does not parity-gate against R output:
``theta`` is recovered via the Bibby-Sørensen martingale estimator, not from
``ARIMA.fit`` coefficients, so the drift does not propagate into the published
Tables A1, A2, A3 in any first-order way.
"""

from __future__ import annotations

import warnings
import numpy as np
import pandas as pd
from statsmodels.tsa.arima.model import ARIMA
from .ARMA_model_internals import (
    arma_expectation,
    arma_variance,
    arma_next_step,
    arma_filter,
    arma_companion_matrix,
    arma_vector_b,
)

# ---------------------------------------------------------------------------
# ARMAModel
# ---------------------------------------------------------------------------


class ARMAModel:
    """
    ARMA(p, q) model.

    Python translation of the R6 class ``ARMA_modelR6``.
    Fitting delegates to ``statsmodels.tsa.arima.ARIMA``.

    Parameters
    ----------
    ar_order : int
        Autoregressive order. Default 1.
    ma_order : int
        Moving-Average order. Default 1.
    include_intercept : bool
        Include a constant term. Default False.

    State-attribute contract (Phase 0 stabilization)
    -----------------------------------------------
    The constructor and ``fit`` both leave the following attributes as
    ``pd.Series`` (possibly empty) so that the ``coefficients``, ``tidy``,
    and downstream consumers in ``SolarModel`` can index/concatenate them
    without type-drift failures:

    - ``self._intercept`` -- ``pd.Series`` with at most one entry "intercept"
    - ``self._phi``       -- ``pd.Series`` indexed "phi_1", "phi_2", ...
    - ``self._theta``     -- ``pd.Series`` indexed "theta_1", "theta_2", ...
    - ``self._std_errors``-- ``pd.Series`` over the union of the above indices
    - ``self._A``         -- companion matrix (numpy array)
    - ``self._sigma2``    -- innovation **variance** (not std. dev.)

    See Also
    --------
    statsmodels.tsa.arima.model.ARIMA : underlying fitting engine.
    """

    _version = "1.0.2"

    # ------------------------------------------------------------------
    # Constructor
    # ------------------------------------------------------------------

    def __init__(
        self,
        ar_order: int = 1,
        ma_order: int = 1,
        include_intercept: bool = False,
    ) -> None:
        self._include_intercept = include_intercept
        self._ar_order = ar_order
        self._ma_order = ma_order
        self._b = arma_vector_b(ar_order, ma_order)

        # Placeholders populated by fit() — kept as Series so concat works
        # before fit and so update() can iterate (name -> value) pairs.
        self._model = None
        self._intercept = pd.Series(dtype=float)
        self._phi = pd.Series(dtype=float)
        self._theta = pd.Series(dtype=float)
        self._A = arma_companion_matrix(np.array([]), np.array([]))
        self._sigma2 = 1.0
        self._std_errors = pd.Series(dtype=float)

    # ------------------------------------------------------------------
    # Public methods
    # ------------------------------------------------------------------

    def fit(self, x: np.ndarray) -> "ARMAModel":
        """
        Fit the ARMA model. Wraps statsmodels ARIMA to mimic R's stats::arima.

        Notes
        -----
        ``ARIMAResults.params`` in modern statsmodels is a positional
        ``numpy.ndarray`` (not a ``pandas.Series``); the names live in
        ``param_names``. ``sigma2`` is also a parameter (last entry).
        """
        trend = "c" if self._include_intercept else "n"

        sm_model = ARIMA(
            np.asarray(x), order=(self._ar_order, 0, self._ma_order), trend=trend
        )
        sm_results = sm_model.fit()
        self._model = sm_results

        # statsmodels >= 0.13 returns ndarrays for params/bse and uses the
        # separate ``param_names`` list for symbolic access.
        names = list(sm_results.param_names)
        params = pd.Series(np.asarray(sm_results.params), index=names, dtype=float)
        std_errors = pd.Series(np.asarray(sm_results.bse), index=names, dtype=float)

        # ---------- Intercept ----------
        if self._include_intercept:
            intercept_key = "const" if "const" in names else names[0]
            self._intercept = pd.Series(
                {"intercept": float(params[intercept_key])}, dtype=float
            )
            ic_se = pd.Series(
                {"intercept": float(std_errors[intercept_key])}, dtype=float
            )
        else:
            self._intercept = pd.Series(dtype=float)
            ic_se = pd.Series(dtype=float)

        # ---------- AR phi ----------
        if self._ar_order > 0:
            ar_keys = [k for k in names if k.startswith("ar.")]
            phi_vals = [float(params[k]) for k in ar_keys]
            phi_se_vals = [float(std_errors[k]) for k in ar_keys]
            phi_index = [f"phi_{i}" for i in range(1, len(phi_vals) + 1)]
            self._phi = pd.Series(phi_vals, index=phi_index, dtype=float)
            phi_se = pd.Series(phi_se_vals, index=phi_index, dtype=float)
        else:
            self._phi = pd.Series(dtype=float)
            phi_se = pd.Series(dtype=float)

        # ---------- MA theta ----------
        if self._ma_order > 0:
            ma_keys = [k for k in names if k.startswith("ma.")]
            theta_vals = [float(params[k]) for k in ma_keys]
            theta_se_vals = [float(std_errors[k]) for k in ma_keys]
            theta_index = [f"theta_{i}" for i in range(1, len(theta_vals) + 1)]
            self._theta = pd.Series(theta_vals, index=theta_index, dtype=float)
            theta_se = pd.Series(theta_se_vals, index=theta_index, dtype=float)
        else:
            self._theta = pd.Series(dtype=float)
            theta_se = pd.Series(dtype=float)

        # Concatenated standard errors, ordered intercept -> phi -> theta to
        # match the `coefficients` property.
        self._std_errors = pd.concat([ic_se, phi_se, theta_se])

        # Companion matrix and innovation variance.
        self._A = arma_companion_matrix(self._phi.values, self._theta.values)
        # Pull sigma2 from the parameter vector (statsmodels appends it last).
        if "sigma2" in names:
            self._sigma2 = float(params["sigma2"])
        else:
            # Fallback for older statsmodels that exposes the attribute.
            self._sigma2 = float(getattr(sm_results, "sigma2", 1.0))

        return self

    def filter(self, x: np.ndarray) -> dict[str, np.ndarray]:
        """Filter the time series, returning fitted values and residuals."""
        intercept_val = (
            float(self._intercept["intercept"])
            if "intercept" in self._intercept.index
            else 0.0
        )
        return arma_filter(np.asarray(x), self._A, self._b, intercept_val)

    def next_step(
        self, x: np.ndarray, n_ahead: int = 1, eps: float = 0.0
    ) -> np.ndarray:
        """Forecast ``n_ahead`` steps from the current state vector."""
        intercept_val = (
            float(self._intercept["intercept"])
            if "intercept" in self._intercept.index
            else 0.0
        )
        return arma_next_step(
            n_ahead, np.asarray(x), self._A, self._b, intercept_val, eps
        )

    def expectation(self, h: int = 1, X0: np.ndarray | None = None) -> np.ndarray:
        """``h``-step-ahead expected value."""
        if X0 is None:
            X0 = np.zeros(sum(self.order.values()))
        intercept_val = (
            float(self._intercept["intercept"])
            if "intercept" in self._intercept.index
            else 0.0
        )
        return arma_expectation(h, np.asarray(X0), self._A, self._b, intercept_val)

    def variance(self, h: int = 1, sigma2: float = 1.0) -> np.ndarray:
        """``h``-step-ahead forecast variance.

        Parameters
        ----------
        h : int
            Forecast horizon.
        sigma2 : float
            Innovation **variance** (not std. dev.). When called from inside a
            fitted ``ARMAModel``, pass ``self.sigma2``.
        """
        return arma_variance(h, self._A, self._b, sigma2)

    def update(self, coefficients) -> None:
        """Update model coefficients in-place.

        ``coefficients`` may be a ``pd.Series`` or any mapping with the
        index/keys ``"intercept"``, ``"phi_1"`` ... and ``"theta_1"`` ....
        Unknown names are ignored with a warning. Setting any coefficient
        invalidates that name's standard error (set to NaN).
        """
        if coefficients is None:
            return

        if hasattr(coefficients, "items"):
            items = list(coefficients.items())
        else:
            # Treat array-like as positional update over current coefficient names.
            items = list(zip(self.coefficients.index, list(coefficients)))

        if len(items) != len(self.coefficients):
            warnings.warn(
                "ARMAModel.update(): length of new `coefficients` does not "
                "match the current coefficients."
            )

        for name, val in items:
            if name in self._intercept.index:
                self._intercept[name] = float(val)
            elif name in self._phi.index:
                self._phi[name] = float(val)
            elif name in self._theta.index:
                self._theta[name] = float(val)
            else:
                warnings.warn(f"ARMAModel.update(): unknown coefficient '{name}' ignored.")
                continue
            if name in self._std_errors.index:
                self._std_errors[name] = np.nan

        self._A = arma_companion_matrix(self._phi.values, self._theta.values)

    def update_std_errors(self, std_errors: pd.Series | None = None) -> None:
        """Update parameter standard errors in-place."""
        if std_errors is None or len(std_errors) == 0:
            return
        if len(std_errors) != len(self._std_errors):
            warnings.warn("ARMAModel.update_std_errors(): length mismatch.")
        for name, val in std_errors.items():
            if name in self._std_errors.index:
                self._std_errors[name] = float(val)

    def update_sigma2(self, sigma2: float | None = None) -> None:
        """Update the innovation variance."""
        if sigma2 is not None:
            self._sigma2 = float(sigma2)

    # ------------------------------------------------------------------
    # __str__ / __repr__
    # ------------------------------------------------------------------

    def __str__(self) -> str:
        GREEN = "\033[1;32m"
        RED = "\033[1;31m"
        RESET = "\033[0m"

        def col_bool(flag: bool) -> str:
            return f"{GREEN if flag else RED}{flag}{RESET}"

        def fmt_param(val: float, se: float) -> str:
            se_str = f"{se:.3g}" if se is not None and not np.isnan(se) else "NA"
            return f"{val:.4g} ({se_str})"

        def fmt_block(series: pd.Series) -> str:
            parts = []
            for name in series.index:
                se = (
                    float(self._std_errors[name])
                    if name in self._std_errors.index
                    else np.nan
                )
                parts.append(f"{name}={fmt_param(float(series[name]), se)}")
            return ", ".join(parts) if parts else "(none)"

        ic_se = (
            float(self._std_errors["intercept"])
            if "intercept" in self._std_errors.index
            else np.nan
        )
        ic_val = (
            float(self._intercept["intercept"])
            if "intercept" in self._intercept.index
            else 0.0
        )
        fmt_ic = fmt_param(ic_val, ic_se)

        model_name = f"ARMA({self._ar_order}, {self._ma_order})"
        sep = "-" * 50
        return "\n".join(
            [
                f"--------------------- {model_name} ---------------------",
                f"Include Intercept : {col_bool(self._include_intercept)}",
                f"AR                : {col_bool(self._ar_order != 0)}",
                f"MA                : {col_bool(self._ma_order != 0)}",
                f"Version           : {self._version}",
                sep,
                f"Intercept         : {fmt_ic}",
                f"AR parameters     : {fmt_block(self._phi)}",
                f"MA parameters     : {fmt_block(self._theta)}",
            ]
        )

    def __repr__(self) -> str:
        return self.__str__()

    # ------------------------------------------------------------------
    # Properties (read-only — mirror R active bindings)
    # ------------------------------------------------------------------

    @property
    def model(self):
        """Fitted ``statsmodels`` ARIMAResults object."""
        return self._model

    @property
    def ar_order(self) -> int:
        """Autoregressive order."""
        return self._ar_order

    @property
    def ma_order(self) -> int:
        """Moving-Average order."""
        return self._ma_order

    @property
    def order(self) -> dict[str, int]:
        """``{'AR': p, 'MA': q}``."""
        return {"AR": self._ar_order, "MA": self._ma_order}

    @property
    def intercept(self) -> pd.Series:
        """Named Series ``{'intercept': value}`` (empty if not included)."""
        return self._intercept

    @property
    def phi(self) -> pd.Series:
        """AR parameters — named Series ``(phi_1, phi_2, ...)``."""
        return self._phi

    @property
    def theta(self) -> pd.Series:
        """MA parameters — named Series ``(theta_1, theta_2, ...)``."""
        return self._theta

    @property
    def coefficients(self) -> pd.Series:
        """All parameters concatenated: intercept -> AR -> MA."""
        return pd.concat([self._intercept, self._phi, self._theta])

    @property
    def std_errors(self) -> pd.Series:
        """Standard errors matching ``coefficients``."""
        return self._std_errors

    @property
    def sigma2(self) -> float:
        """Innovation **variance** (not std. dev.)."""
        return self._sigma2

    @property
    def A(self) -> np.ndarray:
        """Companion (state-transition) matrix."""
        return self._A

    @property
    def b(self) -> np.ndarray:
        """Unit selection vector for the innovation."""
        return self._b

    @property
    def tidy(self) -> pd.DataFrame:
        """DataFrame with columns: ``term``, ``estimate``, ``std_error``."""
        coefs = self.coefficients
        # Align std_errors index to coefs index, NaN where missing.
        ses = self._std_errors.reindex(coefs.index).values
        return pd.DataFrame(
            {
                "term": coefs.index,
                "estimate": coefs.values,
                "std_error": ses,
            }
        )
