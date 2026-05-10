"""
ARMA(p, q) Model — Python translation of the R6 class ARMA_modelR6.

Dependencies: numpy, pandas, statsmodels
Version: 1.0.1
"""

from __future__ import annotations

import re
import warnings
import numpy as np
import pandas as pd
from statsmodels.tsa.arima.model import ARIMA
from .ARMA_model_internals import arma_expectation, arma_variance, arma_next_step, arma_filter, arma_companion_matrix, arma_vector_b

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

    Notes
    -----
    Version 1.0.1

    See Also
    --------
    statsmodels.tsa.arima.model.ARIMA : underlying fitting engine.
    """

    _version = "1.0.1"

    # ------------------------------------------------------------------
    # Constructor
    # ------------------------------------------------------------------

    def __init__(self, ar_order: int = 1, ma_order: int = 1, include_intercept: bool = False,) -> None:
        self._include_intercept = include_intercept
        self._ar_order = ar_order
        self._ma_order = ma_order
        self._b = arma_vector_b(ar_order, ma_order)

        # Placeholders populated by fit()
        self._model = None
        self._intercept = 0.0
        self._phi = pd.Series(dtype=float)
        self._theta = pd.Series(dtype=float)
        self._A = arma_companion_matrix(np.array([]), np.array([]))
        self._sigma2 = 1.0
        self._std_errors = None

        self._coefficients = None

    # ------------------------------------------------------------------
    # Public methods
    # ------------------------------------------------------------------

    def fit(self, x: np.ndarray):
        """
        Fit the ARMA model. Wraps statsmodels ARIMA to mimic R's stats::arima.
        """
        # statsmodels uses order=(p, d, q) and trend='c' for intercept, 'n' for none.
        trend = 'c' if self._include_intercept else 'n'
        
        # To strictly mimic R's arima(..., method="CSS"), we might need conditional likelihood, 
        # but statsmodels defaults to exact state-space MLE which is generally superior.
        sm_model = ARIMA(x, order=(self._ar_order, 0, self._ma_order), trend=trend)
        sm_results = sm_model.fit()
        self._model = sm_results
        
        # Extract parameters and standard errors
        params = sm_results.params
        std_errors = sm_results.bse
        
        # Parse Intercept
        if self._include_intercept:
            # statsmodels typically names this 'const'
            intercept_key = 'const' if 'const' in params.index else params.index[0]
            self._intercept = params[intercept_key]
            self._coefficients['intercept'] = self._intercept
            self._std_errors['intercept'] = std_errors[intercept_key]
        else:
            self._intercept = 0.0

        # Parse AR coefficients
        phi_list, phi_se_list = [], []
        if self._ar_order > 0:
            ar_keys = [k for k in params.index if 'ar.' in k]
            for i, k in enumerate(ar_keys, 1):
                val = params[k]
                se = std_errors[k]
                name = f"phi_{i}"
                phi_list.append(val)
                self._coefficients[name] = val
                self._std_errors[name] = se
        self._phi = np.array(phi_list)

        # Parse MA coefficients
        theta_list, theta_se_list = [], []
        if self._ma_order > 0:
            ma_keys = [k for k in params.index if 'ma.' in k]
            for i, k in enumerate(ma_keys, 1):
                val = params[k]
                se = std_errors[k]
                name = f"theta_{i}"
                theta_list.append(val)
                self._coefficients[name] = val
                self._std_errors[name] = se
        self._theta = np.array(theta_list)

        # Update Companion Matrix and Variance
        self._a_matrix = arma_companion_matrix(self._phi, self._theta)
        self._sigma2 = sm_results.sigma2

    def filter(self, x: np.ndarray) -> dict[str, np.ndarray]:
        """
        Filter the time series, returning fitted values and residuals.
        """
        return arma_filter(np.asarray(x), self._A, self._b, self._intercept)

    def next_step(self, x: np.ndarray, n_ahead: int = 1, eps: float = 0.0) -> np.ndarray:
        """
        Forecast `n_ahead` steps from the current state vector.
        """
        return arma_next_step(
            n_ahead, np.asarray(x), self._A, self._b, self._intercept, eps
        )

    def expectation(self, h: int = 1, X0: np.ndarray | None = None) -> np.ndarray:
        """
        `h`-step-ahead expected value.
        """
        if X0 is None:
            X0 = np.zeros(sum(self.order.values()))
        return arma_expectation(h, np.asarray(X0), self._A, self._b, self._intercept)

    def variance(self, h: int = 1, sigma2: float = 1.0) -> np.ndarray:
        """
        h-step-ahead forecast variance.

        Parameters
        ----------
        h : int
        sigma2 : float
            Residual std. deviation.
        """
        return arma_variance(h, self.A, self.b, sigma2)

    def update(self, coefficients) -> None:
        """
        Update model coefficients in-place.
        """
        if coefficients is None:
            return

        if len(coefficients) != len(self.coefficients):
            warnings.warn(
                "ARMAModel.update(): length of new `coefficients` does not "
                "match the current coefficients."
            )

        for name, val in coefficients.items():
            if name in self.coefficients.index:
                self.coefficients[name] = val
                self._std_errors[name] = np.nan
        
        self._A = arma_companion_matrix(self._phi.values, self._theta.values)

    def update_std_errors(self, std_errors: pd.Series | None = None) -> None:
        """
        Update parameter standard errors in-place.

        Parameters
        ----------
        std_errors : pd.Series
        """
        if (std_errors is None or len(std_errors) == 0 or self._std_errors is None):
            return

        if len(std_errors) != len(self._std_errors):
            warnings.warn("ARMAModel.update_std_errors(): length mismatch.")

        for name, val in std_errors.items():
            if name in self._std_errors.index:
                self._std_errors[name] = val

    def update_sigma2(self, sigma2: float | None = None) -> None:
        """Update the residual std. deviation."""
        if sigma2 is not None:
            self._sigma2 = sigma2

    def __str__(self) -> str:
        GREEN = "\033[1;32m"
        RED   = "\033[1;31m"
        RESET = "\033[0m"

        def col_bool(flag: bool) -> str:
            return f"{GREEN if flag else RED}{flag}{RESET}"

        def fmt_param(val: float, se: float) -> str:
            se_str = f"{se:.3g}" if not np.isnan(se) else "NA"
            return f"{val:.4g} ({se_str})"

        ic_se = (
            float(self._std_errors["intercept"])
            if self._std_errors is not None and "intercept" in self._std_errors.index
            else np.nan
        )
        fmt_ic = fmt_param(float(self._intercept["intercept"]), ic_se)

        def fmt_block(series: pd.Series) -> str:
            parts = []
            for name in series.index:
                se = (
                    float(self._std_errors[name])
                    if self._std_errors is not None and name in self._std_errors.index
                    else np.nan
                )
                parts.append(f"{name}={fmt_param(float(series[name]), se)}")
            return ", ".join(parts)

        model_name = f"ARMA({self._ar_order}, {self._ma_order})"
        sep = "-" * 50
        return "\n".join([
            f"--------------------- {model_name} ---------------------",
            f"Include Intercept : {col_bool(self._include_intercept)}",
            f"AR                : {col_bool(self._ar_order != 0)}",
            f"MA                : {col_bool(self._ma_order != 0)}",
            f"Version           : {self._version}",
            sep,
            f"Intercept         : {fmt_ic}",
            f"AR parameters     : {fmt_block(self._phi)}",
            f"MA parameters     : {fmt_block(self._theta)}",
        ])

    def __repr__(self) -> str:
        return self.__str__()

    # ------------------------------------------------------------------
    # Properties  (read-only — mirror R active bindings)
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
        """Named Series ``{'intercept': value}``."""
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
        """All parameters concatenated: intercept → AR → MA."""
        return pd.concat([self._intercept, self._phi, self._theta])

    @property
    def std_errors(self) -> pd.Series | None:
        """Standard errors matching ``coefficients``."""
        return self._std_errors

    @property
    def sigma2(self) -> float:
        """Residual std. deviation (square root of the innovation variance)."""
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

        return pd.DataFrame({
            "term": self.coefficients.index,
            "estimate": self.coefficients.values,
            "std_error": self._std_errors.values,
        })