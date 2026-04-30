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


# ---------------------------------------------------------------------------
# Helper functions  (replace with your own implementations if needed)
# ---------------------------------------------------------------------------

def ARMA_vector_b(ar_order: int, ma_order: int) -> np.ndarray:
    """
    Unit selection vector of length (ar_order + ma_order).
    Marks the position where the innovation enters the state vector.
    Returns an empty array when both orders are 0, or when there is no
    MA component (ma_order == 0).
    """
    size = ar_order + ma_order
    b = np.zeros(size)
    # The residual enters only if there is an MA part
    if ma_order > 0 and ar_order < size:
        b[ar_order] = 1.0
    return b


def ARMA_companion_matrix(phi: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """
    Build the companion / state-transition matrix A for the ARMA model.

    State vector layout: [y_{t-1}, ..., y_{t-p}, eps_{t-1}, ..., eps_{t-q}]
    """
    p = len(phi)
    q = len(theta)
    n = p + q
    if n == 0:
        return np.empty((0, 0))

    A = np.zeros((n, n))

    # AR block (top-left p×p companion)
    if p > 0:
        A[0, :p] = phi
        if p > 1:
            A[1:p, :p - 1] = np.eye(p - 1)

    # MA block (top rows, columns p : p+q)
    if q > 0:
        A[0, p : p + q] = theta
        if q > 1:
            A[p + 1 : p + q, p : p + q - 1] = np.eye(q - 1)

    return A


def ARMA_filter(
    x: np.ndarray,
    A: np.ndarray,
    b: np.ndarray,
    intercept: float,
) -> dict[str, np.ndarray]:
    """
    Recursive state-space filter.
    Returns dict with keys 'fitted' (one-step predictions) and 'residuals'.
    """
    T = len(x)
    n = len(b)
    state = np.zeros(n)
    fitted = np.zeros(T)
    residuals = np.zeros(T)

    for t in range(T):
        y_hat = float(A[0] @ state) + intercept if n > 0 else intercept
        eps = x[t] - y_hat
        fitted[t] = y_hat
        residuals[t] = eps
        if n > 0:
            state = A @ state + b * eps

    return {"fitted": fitted, "residuals": residuals}


def ARMA_next_step(
    n_ahead: int,
    x: np.ndarray,
    A: np.ndarray,
    b: np.ndarray,
    intercept: float,
    eps: float = 0.0,
) -> np.ndarray:
    """
    Forecast n_ahead steps given current state vector x.
    eps is the most recent realised residual (used for step 1 only).
    """
    n = len(b)
    if n == 0:
        return np.full(n_ahead, intercept)

    state = x.copy().astype(float)
    forecasts = np.zeros(n_ahead)

    for h in range(n_ahead):
        current_eps = eps if h == 0 else 0.0
        state = A @ state + b * current_eps
        forecasts[h] = float(A[0] @ state) + intercept

    return forecasts


def ARMA_expectation(
    h: int,
    X0: np.ndarray,
    A: np.ndarray,
    b: np.ndarray,
    intercept: float,
) -> np.ndarray:
    """h-step-ahead expected values starting from state X0."""
    n = len(b)
    if n == 0:
        return np.full(h, intercept)

    expectations = np.zeros(h)
    A_power = np.eye(n)

    for step in range(h):
        A_power = A_power @ A
        expectations[step] = float(A[0] @ A_power @ X0) + intercept

    return expectations


def ARMA_variance(
    h: int,
    A: np.ndarray,
    b: np.ndarray,
    sigma2: float,
) -> np.ndarray:
    """h-step-ahead forecast variances via the VMA representation."""
    n = len(b)
    if n == 0:
        return np.full(h, sigma2 ** 2)

    variances = np.zeros(h)
    psi = b.copy().astype(float)

    for step in range(h):
        variances[step] = (sigma2 ** 2) * float(psi @ psi)
        psi = A @ psi

    return variances


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

    def __init__(
        self,
        ar_order: int = 1,
        ma_order: int = 1,
        include_intercept: bool = False,
    ) -> None:
        self._include_intercept = include_intercept
        self._ar_order = ar_order
        self._ma_order = ma_order
        self._b = ARMA_vector_b(ar_order, ma_order)

        # Placeholders populated by fit()
        self._model = None
        self._intercept = pd.Series({"intercept": 0.0})
        self._phi = pd.Series(dtype=float)
        self._theta = pd.Series(dtype=float)
        self._A = ARMA_companion_matrix(np.array([]), np.array([]))
        self._sigma2: float = 1.0
        self._std_errors: pd.Series | None = None

    # ------------------------------------------------------------------
    # Public methods
    # ------------------------------------------------------------------

    def fit(self, x: np.ndarray) -> None:
        """
        Fit the model using ``statsmodels.tsa.arima.ARIMA``.

        Mirrors R's ``arima(x, order=c(p,0,q), include.mean=..., method='CSS')``.

        Parameters
        ----------
        x : array-like
            Time series.
        """
        x = np.asarray(x, dtype=float)
        trend = "c" if self._include_intercept else "n"

        result = ARIMA(
            x,
            order=(self._ar_order, 0, self._ma_order),
            trend=trend,
        ).fit(method="innovations_mle")

        # statsmodels returns plain numpy arrays; use param_names for lookup
        param_names: list[str] = result.param_names
        params: np.ndarray = result.params
        bse: np.ndarray = result.bse

        # ---- intercept --------------------------------------------------
        intercept_val = 0.0
        intercept_se = np.nan
        ic_idx = [i for i, n in enumerate(param_names) if n in ("const", "intercept")]
        if self._include_intercept and ic_idx:
            idx = ic_idx[0]
            intercept_val = float(params[idx])
            intercept_se = float(bse[idx])
        self._intercept = pd.Series({"intercept": intercept_val})

        # ---- AR coefficients -------------------------------------------
        ar_idx = [i for i, n in enumerate(param_names) if n.startswith("ar.")]
        phi_names = [f"phi_{i+1}" for i in range(len(ar_idx))]
        self._phi = pd.Series(params[ar_idx], index=phi_names)
        phi_se = pd.Series(bse[ar_idx], index=phi_names)

        # ---- MA coefficients -------------------------------------------
        ma_idx = [i for i, n in enumerate(param_names) if n.startswith("ma.")]
        theta_names = [f"theta_{i+1}" for i in range(len(ma_idx))]
        self._theta = pd.Series(params[ma_idx], index=theta_names)
        theta_se = pd.Series(bse[ma_idx], index=theta_names)

        # ---- store -------------------------------------------------------
        self._model = result
        self._A = ARMA_companion_matrix(self._phi.values, self._theta.values)

        # R's sigma2 field stores the std dev (sqrt of variance)
        sigma2_idx = param_names.index("sigma2")
        self._sigma2 = float(np.sqrt(params[sigma2_idx]))

        self._std_errors = pd.concat([
            pd.Series({"intercept": intercept_se}),
            phi_se,
            theta_se,
        ])

    def filter(self, x: np.ndarray) -> dict[str, np.ndarray]:
        """
        Filter the time series, returning fitted values and residuals.

        Parameters
        ----------
        x : array-like

        Returns
        -------
        dict with keys ``'fitted'`` and ``'residuals'``.
        """
        return ARMA_filter(np.asarray(x), self.A, self.b, self.intercept_value)

    def next_step(
        self,
        x: np.ndarray,
        n_ahead: int = 1,
        eps: float = 0.0,
    ) -> np.ndarray:
        """
        Forecast ``n_ahead`` steps from the current state vector.

        Parameters
        ----------
        x : array-like
            State vector of length ``p + q``.
        n_ahead : int
        eps : float
            Most recent realised residual (used at step 1 only).
        """
        return ARMA_next_step(
            n_ahead, np.asarray(x), self.A, self.b, self.intercept_value, eps
        )

    def expectation(self, h: int = 1, X0: np.ndarray | None = None) -> np.ndarray:
        """
        h-step-ahead expected value.

        Parameters
        ----------
        h : int
        X0 : array-like, optional
            Initial state of length ``p + q``. Defaults to zeros.
        """
        if X0 is None:
            X0 = np.zeros(sum(self.order.values()))
        return ARMA_expectation(h, np.asarray(X0), self.A, self.b, self.intercept_value)

    def variance(self, h: int = 1, sigma2: float = 1.0) -> np.ndarray:
        """
        h-step-ahead forecast variance.

        Parameters
        ----------
        h : int
        sigma2 : float
            Residual std. deviation.
        """
        return ARMA_variance(h, self.A, self.b, sigma2)

    def update(self, coefficients: pd.Series | None = None) -> None:
        """
        Update model coefficients in-place.

        Parameters
        ----------
        coefficients : pd.Series
            Named series. Unrecognised names are silently ignored.
            Matching the R behaviour, the std. errors for updated
            parameters are set to NaN.
        """
        if coefficients is None:
            return

        new_coefs = self.coefficients.copy()
        old_names = set(new_coefs.index)

        if len(coefficients) != len(new_coefs):
            warnings.warn(
                "ARMAModel.update(): length of new `coefficients` does not "
                "match the current coefficients."
            )

        for name, val in coefficients.items():
            if name in old_names:
                new_coefs[name] = val
                if self._std_errors is not None and name in self._std_errors.index:
                    self._std_errors[name] = np.nan

        if self._include_intercept:
            self._intercept = pd.Series({"intercept": new_coefs["intercept"]})

        if self._ar_order > 0:
            phi_names = [n for n in new_coefs.index if re.search(r"phi", n)]
            self._phi = new_coefs[phi_names]

        if self._ma_order > 0:
            theta_names = [n for n in new_coefs.index if re.search(r"theta", n)]
            self._theta = new_coefs[theta_names]

        self._A = ARMA_companion_matrix(self._phi.values, self._theta.values)

    def update_std_errors(self, std_errors: pd.Series | None = None) -> None:
        """
        Update parameter standard errors in-place.

        Parameters
        ----------
        std_errors : pd.Series
        """
        if (
            std_errors is None
            or len(std_errors) == 0
            or self._std_errors is None
        ):
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

    def __repr__(self) -> str:
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
    def intercept_value(self) -> float:
        """Intercept as a plain float (0.0 when not included)."""
        return float(self._intercept["intercept"])

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
        coefs = self.coefficients
        ses = (
            self._std_errors.reindex(coefs.index)
            if self._std_errors is not None
            else pd.Series([np.nan] * len(coefs), index=coefs.index)
        )
        return pd.DataFrame({
            "term": coefs.index,
            "estimate": coefs.values,
            "std_error": ses.values,
        })