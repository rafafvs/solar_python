"""
Solar model diagnostic tests — Python translation.

Functions:
    solar_model_test_distribution  — per-month KS test on residuals
    solar_model_test_autocorr      — autocorrelation tests (BG / BP / LB)
    solar_model_tests              — combined autocorr + distribution tests
    solar_model_test_PIT           — Probability Integral Transform
    solar_model_test_LPD           — Log Predictive Density
    solar_model_test_forecast      — point-forecast metrics (RMSE, MAE, …)
    solar_model_test_pricing       — option pricing accuracy metrics
"""

from __future__ import annotations

import warnings
import numpy as np
import pandas as pd
from scipy import stats as sci_stats


# ---------------------------------------------------------------------------
# solar_model_test_distribution
# ---------------------------------------------------------------------------

def solar_model_test_distribution(
    model,
    H0: str = "gm",
    ci: float = 0.05,
    min_quantile: float = 0.025,
    max_quantile: float = 0.985,
    type: str = "train",
) -> pd.DataFrame:
    """
    Per-month Kolmogorov-Smirnov test on the standardised residuals
    of a fitted SolarModel.

    Parameters
    ----------
    model        : SolarModel
    H0           : str, 'gm' (Gaussian mixture) or 'norm' (normality)
    ci           : float, significance level
    min_quantile : float, lower quantile trim for the KS test
    max_quantile : float, upper quantile trim for the KS test
    type         : str, 'train' | 'test' | 'full'

    Returns
    -------
    pd.DataFrame with columns:
        Month, test_H0, data, statistic, p_value, H0, result
    """
    if H0 not in ("gm", "norm"):
        raise ValueError("`H0` must be 'gm' or 'norm'.")
    if type not in ("train", "test", "full"):
        raise ValueError("`type` must be 'train', 'test', or 'full'.")

    data = model.data.copy()
    if type == "train":
        data = data[(data["isTrain"]) & (data["weights"] != 0)]
    elif type == "test":
        data = data[~data["isTrain"]]

    rows = []
    for month in range(1, 13):
        x = data[data["Month"] == month]["u_tilde"].values

        if H0 == "norm":
            cdf_Yt = lambda v: sci_stats.norm.cdf(v, loc=np.mean(x), scale=np.std(x, ddof=1))
        else:
            gm = model.NM_model.model[month - 1]   # 0-indexed list
            cdf_Yt = lambda v: pmixnorm(v, gm["means"], gm["sd"], gm["p"])

        row = ks_test(x, cdf_Yt, ci=ci, min_quantile=min_quantile, max_quantile=max_quantile)
        rows.append(row)

    tests = pd.DataFrame(rows)
    tests.insert(0, "Month",   range(1, 13))
    tests.insert(1, "test_H0", H0)
    tests.insert(2, "data",    type)
    tests["result"] = tests["H0"].apply(lambda h: "Passed" if h == "Non-Rejected" else "Not-passed")
    return tests


# ---------------------------------------------------------------------------
# solar_model_test_autocorr
# ---------------------------------------------------------------------------

def solar_model_test_autocorr(
    model,
    lag_max: int = 3,
    ci: float = 0.05,
    method: str = "bg",
    type: str = "train",
) -> pd.DataFrame:
    """
    Autocorrelation test on all residual series of a fitted SolarModel.

    Parameters
    ----------
    model   : SolarModel
    lag_max : int, maximum lag
    ci      : float, significance level
    method  : str, 'bg' (Breusch-Godfrey) | 'bp' (Box-Pierce) | 'lb' (Ljung-Box)
    type    : str, 'train' | 'test' | 'full'

    Returns
    -------
    pd.DataFrame with one row per tested series.
    """
    if method not in ("bg", "bp", "lb"):
        raise ValueError("`method` must be 'bg', 'bp', or 'lb'.")
    if type not in ("train", "test", "full"):
        raise ValueError("`type` must be 'train', 'test', or 'full'.")

    data = model.data.copy()
    if type == "train":
        data = data[(data["isTrain"]) & (data["weights"] != 0)]
    elif type == "test":
        data = data[~data["isTrain"]]

    # Standardise mixture residuals: u_tilde = (u - μ) / √v
    nm_moments = model.NM_model.moments[["Month", "mean", "variance"]]
    data = data.merge(nm_moments, on="Month", how="left")
    data["u"]       = data["u_tilde"]
    data["u_tilde"] = (data["u"] - data["mean"]) / np.sqrt(data["variance"])

    # Squared series
    for col in ("Yt_tilde", "eps", "eps_tilde", "u", "u_tilde"):
        data[f"{col}2"] = data[col] ** 2

    # Expected outcome per target (mirrors R's `expected` vector)
    targets = ["Yt_tilde", "eps", "eps_tilde", "u", "u_tilde",
               "Yt_tilde2", "eps2", "eps_tilde2", "u_2", "u_tilde2"]
    expected = {
        "Yt_tilde":  "Rejected",
        "eps":       "Not-rejected",
        "eps_tilde": "Not-rejected",
        "u":         "Not-rejected",
        "u_tilde":   "Not-rejected",
        "Yt_tilde2": "Rejected",
        "eps2":      "Rejected",
        "eps_tilde2":"Rejected",
        "u_2":       "Not-rejected",
        "u_tilde2":  "Not-rejected",
    }

    # ---- Inner test helpers ----

    def _bg_test(series: np.ndarray, lag: int) -> dict:
        """Breusch-Godfrey via auxiliary OLS regression of residuals on lags."""
        from statsmodels.regression.linear_model import OLS
        from statsmodels.tools import add_constant
        n = len(series)
        # Build lag matrix
        X = np.column_stack([series[lag - k : n - k] for k in range(1, lag + 1)])
        e = series[lag:]
        X = add_constant(X)
        res = OLS(e, X).fit()
        stat = n * res.rsquared
        pval = 1 - sci_stats.chi2.cdf(stat, df=lag)
        return {"statistic": stat, "p_value": pval, "lags": lag, "H0_method": "Breusch-Godfrey"}

    def _box_test(series: np.ndarray, lag: int, box_type: str) -> dict:
        """Box-Pierce or Ljung-Box via statsmodels."""
        from statsmodels.stats.diagnostic import acorr_ljungbox
        lb = box_type == "lb"
        result = acorr_ljungbox(series, lags=[lag], boxpierce=True, return_df=True)
        if lb:
            stat = float(result["lb_stat"].iloc[0])
            pval = float(result["lb_pvalue"].iloc[0])
            method_name = "Ljung-Box"
        else:
            stat = float(result["bp_stat"].iloc[0])
            pval = float(result["bp_pvalue"].iloc[0])
            method_name = "Box-Pierce"
        return {"statistic": stat, "p_value": pval, "lags": lag, "H0_method": method_name}

    rows = []
    for target in targets:
        if target not in data.columns:
            continue
        series = data[target].dropna().values
        exp    = expected[target]

        if method == "bg":
            res = _bg_test(series, lag_max)
        elif method == "bp":
            res = _box_test(series, lag_max, "bp")
        else:
            res = _box_test(series, lag_max, "lb")

        h0_label = "Not-rejected" if res["p_value"] > ci else "Rejected"
        rows.append({
            "target":    target,
            "statistic": round(res["statistic"], 5),
            "p_value":   round(res["p_value"],   5),
            "H0":        h0_label,
            "lags":      res["lags"],
            "method":    res["H0_method"],
            "result":    "passed" if h0_label == exp else "Not-passed",
        })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# solar_model_tests
# ---------------------------------------------------------------------------

def solar_model_tests(
    model,
    lags: list[int] = None,
    ci: float = 0.05,
    min_quantile: float = 0.025,
    max_quantile: float = 0.985,
    method: str = "bg",
    type: str = "train",
) -> dict:
    """
    Combined autocorrelation + distribution test battery.

    Parameters
    ----------
    model        : SolarModel
    lags         : list[int], lags for autocorrelation tests (default [7])
    ci           : float, significance level
    min_quantile : float, lower quantile trim for KS tests
    max_quantile : float, upper quantile trim for KS tests
    method       : str, autocorrelation method ('bg' | 'bp' | 'lb')
    type         : str, 'train' | 'test' | 'full'

    Returns
    -------
    dict with keys:
        autocorr  — dict keyed by 'lag_{k}', each a pd.DataFrame
        normality — pd.DataFrame (KS vs normal)
        mixture   — pd.DataFrame (KS vs Gaussian mixture)
    """
    if lags is None:
        lags = [7]

    autocorr_tests = {
        f"lag_{lag}": solar_model_test_autocorr(
            model, lag_max=lag, ci=ci, method=method, type=type
        )
        for lag in lags
    }

    normality_test = solar_model_test_distribution(
        model, H0="norm", ci=ci,
        min_quantile=min_quantile, max_quantile=max_quantile, type=type,
    )
    mixture_test = solar_model_test_distribution(
        model, H0="gm", ci=ci,
        min_quantile=min_quantile, max_quantile=max_quantile, type=type,
    )

    return {
        "autocorr":  autocorr_tests,
        "normality": normality_test,
        "mixture":   mixture_test,
    }


# ---------------------------------------------------------------------------
# solar_model_test_PIT
# ---------------------------------------------------------------------------

def solar_model_test_PIT(
    model,
    ci: float = 0.05,
    type: str = "train",
) -> dict:
    """
    Probability Integral Transform (PIT) uniformity test.

    Transforms observed GHI values through the estimated conditional CDF
    and tests whether the resulting grades u are Uniform(0, 1).

    Parameters
    ----------
    model : SolarModel
    ci    : float, significance level for the KS uniformity test
    type  : str, 'train' | 'test' | 'full'

    Returns
    -------
    dict with keys:
        data — pd.DataFrame with columns (date, Year, Month, Day, u)
        test — pd.DataFrame with KS test result
    """
    if type not in ("train", "test", "full"):
        raise ValueError("`type` must be 'train', 'test', or 'full'.")

    moments = model.moments["conditional"].copy()
    Rt      = model.data["GHI"].values

    if type == "train":
        mask    = (model.data["isTrain"]) & (model.data["weights"] != 0)
        moments = moments[mask.values]
        Rt      = Rt[mask.values]
    elif type == "test":
        mask    = ~model.data["isTrain"]
        moments = moments[mask.values]
        Rt      = Rt[mask.values]

    link = model.spec["transform"]["link"]
    u    = np.empty(len(moments))

    for i, (_, row) in enumerate(moments.iterrows()):
        def cdf_Y(x, _row=row):
            return pmixnorm(
                x,
                mean=[_row["M_Y1"], _row["M_Y0"]],
                sd=[_row["S_Y1"],   _row["S_Y0"]],
                alpha=[_row["p1"],  1 - _row["p1"]],
            )
        u[i] = p_solar_ghi(Rt[i], row["Ct"], row["alpha"], row["beta"], cdf_Y, link=link)

    moments = moments.copy()
    moments["u"] = u

    test_df = pd.DataFrame([
        {"link": link}
        | ks_test(u, lambda x: sci_stats.uniform.cdf(x), ci=ci)
    ])

    return {
        "data": moments[["date", "Year", "Month", "Day", "u"]].reset_index(drop=True),
        "test": test_df,
    }


# ---------------------------------------------------------------------------
# solar_model_test_LPD
# ---------------------------------------------------------------------------

def solar_model_test_LPD(
    model,
    type: str = "train",
) -> pd.DataFrame:
    """
    Mean Log Predictive Density (LPD) of a fitted SolarModel.

    Parameters
    ----------
    model : SolarModel
    type  : str, 'train' | 'test' | 'full'

    Returns
    -------
    pd.DataFrame with columns: type, link, LPD
    """
    if type not in ("train", "test", "full"):
        raise ValueError("`type` must be 'train', 'test', or 'full'.")

    moments = model.moments["conditional"].copy()

    if type == "train":
        mask    = (model.data["isTrain"]) & (model.data["weights"] != 0)
        moments = moments[mask.values]
    elif type == "test":
        mask    = ~model.data["isTrain"]
        moments = moments[mask.values]

    LPD = model.log_lik(moments, target="GHI")
    LPD = LPD[~np.isinf(LPD)]

    return pd.DataFrame([{
        "type": type,
        "link": model.spec["transform"]["link"],
        "LPD":  float(np.nanmean(LPD)),
    }])


# ---------------------------------------------------------------------------
# solar_model_test_forecast
# ---------------------------------------------------------------------------

def solar_model_test_forecast(
    model,
    ci: float = 0.1,
    type: str = "train",
) -> pd.DataFrame:
    """
    Point-forecast and interval-coverage metrics for a fitted SolarModel.

    Metrics: SSE, RMSE, MAE, MAPE, upper/lower VaR violation rates.

    Parameters
    ----------
    model : SolarModel
    ci    : float, nominal coverage (two-sided, so total = 2·ci)
    type  : str, 'train' | 'test' | 'full'

    Returns
    -------
    pd.DataFrame with one row of summary statistics.
    """
    if type not in ("train", "test", "full"):
        raise ValueError("`type` must be 'train', 'test', or 'full'.")

    moments = model.moments["conditional"].copy()

    if type == "train":
        mask    = (model.data["isTrain"]) & (model.data["weights"] != 0)
        moments = moments[mask.values]
    elif type == "test":
        mask    = ~model.data["isTrain"]
        moments = moments[mask.values]

    forecast = solar_model_forecast(model, moments, ci=ci)

    errors = forecast["Rt"].values - forecast["e_Rt"].values

    SSE  = float(np.sum(errors ** 2))
    RMSE = float(np.sqrt(np.mean(errors ** 2)))
    MAE  = float(np.mean(np.abs(errors)))
    MAPE = float(np.mean(np.abs(errors) / forecast["Rt"].values * 100))

    viol_hi = float(np.mean(forecast["Rt"].values > forecast["ci_mix_hi"].values))
    viol_lo = float(np.mean(forecast["Rt"].values < forecast["ci_mix_lo"].values))

    return pd.DataFrame([{
        "type":       type,
        "SSE":        SSE,
        "RMSE":       RMSE,
        "MAE":        MAE,
        "MAPE":       MAPE,
        "ci":         ci * 2,
        "viol_ci_hi": viol_hi,
        "viol_ci_lo": viol_lo,
    }])


# ---------------------------------------------------------------------------
# solar_model_test_pricing
# ---------------------------------------------------------------------------

def solar_model_test_pricing(
    model,
    type: str = "train",
    control=None,
) -> dict:
    """
    Option-pricing accuracy metrics for a fitted SolarModel.

    Computes model-implied vs. realised PUT and CALL statistics at both
    individual-day and annual-index levels.

    Parameters
    ----------
    model   : SolarModel
    type    : str, 'train' | 'test' | 'full'
    control : control object from ``control_solar_option()``; uses defaults
              when None.

    Returns
    -------
    dict with keys: put, call, put_idx, call_idx
        Each value is a pd.DataFrame with one row of summary statistics.
    """
    if type not in ("train", "test", "full"):
        raise ValueError("`type` must be 'train', 'test', or 'full'.")

    if control is None:
        control = control_solar_option()

    moments = model.moments["conditional"].copy()

    if type == "train":
        mask    = (model.data["isTrain"]) & (model.data["weights"] != 0)
        moments = moments[mask.values]
    elif type == "test":
        mask    = ~model.data["isTrain"]
        moments = moments[mask.values]

    put_prices  = solar_option_model(model, moments, put=True,  control_options=control)
    call_prices = solar_option_model(model, moments, put=False, control_options=control)

    def put_call_stats(price: np.ndarray, exercise: np.ndarray, gamma: np.ndarray) -> dict:
        errors = price - gamma
        ex_mask = exercise != 0
        return {
            "SSE":     float(np.sum(errors ** 2)),
            "premium": float(np.sum(price)),
            "payoff":  float(np.sum(gamma)),
            "diff":    float(np.sum(price) - np.sum(gamma)),
            "cPt":     float(np.mean(price[ex_mask]))   if ex_mask.any() else np.nan,
            "cGamma":  float(np.mean(gamma[ex_mask]))   if ex_mask.any() else np.nan,
            "cdiff":   float(np.mean(gamma[ex_mask]) - np.mean(price[ex_mask])) if ex_mask.any() else np.nan,
        }

    put_stats  = put_call_stats(
        put_prices["payoff"]["premium"].values,
        put_prices["payoff"]["exercise"].values,
        put_prices["payoff"]["payoff"].values,
    )
    call_stats = put_call_stats(
        call_prices["payoff"]["premium"].values,
        call_prices["payoff"]["exercise"].values,
        call_prices["payoff"]["payoff"].values,
    )

    # ---- Annual index stats ----
    Pt_hist = solar_option_historical(model, put=True,  control_options=control)["payoff_year"]["premium"].values
    Pt      = put_prices["payoff_year"]["premium"].values
    Gamma_p = put_prices["payoff_year"]["payoff"].values

    put_idx = pd.DataFrame({
        "Pt_hist":          Pt_hist,
        "Pt":               Pt,
        "Gamma":            Gamma_p,
        "diff_Gamma_Pt":    Gamma_p - Pt,
        "diff_Gamma_Pt_hist": Gamma_p - Pt_hist,
        "diff_Pt_hist_Pt":  Pt_hist - Pt,
    })

    Ct_hist = solar_option_historical(model, put=False, control_options=control)["payoff_year"]["premium"].values
    Ct      = call_prices["payoff_year"]["premium"].values
    Gamma_c = call_prices["payoff_year"]["payoff"].values

    call_idx = pd.DataFrame({
        "Ct_hist":          Ct_hist,
        "Ct":               Ct,
        "Gamma":            Gamma_c,
        "diff_Gamma_Ct":    Gamma_c - Ct,
        "diff_Gamma_Ct_hist": Gamma_c - Ct_hist,
        "diff_Ct_hist_Ct":  Ct_hist - Ct,
    })

    K = control.get("K", np.nan)
    n = len(moments)

    return {
        "put":      pd.DataFrame([{"type": type, "K": K, "n": n} | put_stats]),
        "call":     pd.DataFrame([{"type": type, "K": K, "n": n} | call_stats]),
        "put_idx":  pd.concat([pd.DataFrame([{"type": type, "K": K, "n": n}]), put_idx],  axis=1),
        "call_idx": pd.concat([pd.DataFrame([{"type": type, "K": K, "n": n}]), call_idx], axis=1),
    }