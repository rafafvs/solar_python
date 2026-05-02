"""
Solar model prediction functions — Python translation.

Functions:
    solar_model_covariance      — conditional covariance between Y_t and Y_T
    solar_model_predict         — single-date forecast with PDF/CDF/CI
    solar_model_forecast        — iterate forecast over multiple dates
    solar_model_match_params    — match a flat parameter vector into a nested dict
"""

from __future__ import annotations

import warnings
import numpy as np
import pandas as pd
from scipy import integrate as sci_integrate
from scipy.stats import multivariate_normal


# ---------------------------------------------------------------------------
# solar_model_covariance
# ---------------------------------------------------------------------------

def solar_model_covariance(
    t_now,
    mom_t: pd.Series,
    mom_T: pd.Series,
    GARCH,
    NM_model,
    theta: float = 0.0,
    tol: float = 0.01,
) -> pd.Series:
    """
    Compute the conditional covariance between Y_t and Y_T and attach
    joint/marginal mixture PDFs to the horizon moments row.

    Parameters
    ----------
    t_now   : date-like, current date (unused directly; carried for context)
    mom_t   : pd.Series, conditional moments at date t  (single row)
    mom_T   : pd.Series, conditional moments at date T  (single row)
    GARCH   : fitted GARCH object
    NM_model: fitted Normal-Mixture model
    theta   : float, mean-reversion parameter
    tol     : float, cubature tolerance (kept for API parity; cubature not used)

    Returns
    -------
    pd.Series — mom_T augmented with mixture means, variances, covariances,
                and joint/marginal PDF callables.
    """
    mom_T = mom_T.copy()

    t_cond = mom_t["date"]
    t_hor  = mom_T["date"]
    print(f"Covariance {t_cond} - {t_hor}")

    mom_T_1 = solar_moments_path(mom_T, GARCH, NM_model, theta=theta, B=1, t_cond=t_cond)
    mom_T_0 = solar_moments_path(mom_T, GARCH, NM_model, theta=theta, B=0, t_cond=t_cond)

    # ---- Conditional means of Y_T ----
    mom_T["M_YT_11"] = mom_T_1["M_Y1"]
    mom_T["M_YT_10"] = mom_T_1["M_Y0"]
    mom_T["M_YT_01"] = mom_T_0["M_Y1"]
    mom_T["M_YT_00"] = mom_T_0["M_Y0"]

    # ---- Conditional variances of Y_T ----
    mom_T["S_YT_11"] = mom_T_1["S_Y1"]
    mom_T["S_YT_10"] = mom_T_1["S_Y0"]
    mom_T["S_YT_01"] = mom_T_0["S_Y1"]
    mom_T["S_YT_00"] = mom_T_0["S_Y0"]

    # ---- Conditional means / variances of Y_t (independent of B_T) ----
    mom_T["M_Yt_11"] = mom_T["M_Yt_10"] = mom_t["M_Y1"]
    mom_T["M_Yt_01"] = mom_T["M_Yt_00"] = mom_t["M_Y0"]
    mom_T["S_Yt_11"] = mom_T["S_Yt_10"] = mom_t["S_Y1"]
    mom_T["S_Yt_01"] = mom_T["S_Yt_00"] = mom_t["S_Y0"]

    # ---- Covariance terms ----
    data_psi_t   = mom_t["psi_j"][0].copy()
    data_psi_T_1 = mom_T_1["psi_j"][0].copy()
    data_psi_T_1 = data_psi_T_1[data_psi_T_1["date"] <= mom_t["date"]]
    data_psi_T_0 = mom_T_0["psi_j"][0].copy()
    data_psi_T_0 = data_psi_T_0[data_psi_T_0["date"] <= mom_t["date"]]

    t = len(data_psi_t)

    # psi_t_1: weight using component-1 variance at the last row, mixture elsewhere
    variance_col = data_psi_t["variance"].values.copy()
    variance_col_1 = np.concatenate([
        variance_col[: t - 1],
        [data_psi_t["sd1"].values[t - 1] ** 2],
    ])
    variance_col_0 = np.concatenate([
        variance_col[: t - 1],
        [data_psi_t["sd2"].values[t - 1] ** 2],
    ])
    data_psi_t = data_psi_t.copy()
    data_psi_t["psi_t_1"] = data_psi_t["psi_j"] * data_psi_t["sigma_bar"] ** 2 * variance_col_1
    data_psi_t["psi_t_0"] = data_psi_t["psi_j"] * data_psi_t["sigma_bar"] ** 2 * variance_col_0

    cv_given_1 = float((data_psi_t["psi_t_1"].values * data_psi_T_1["psi_j"].values).sum())
    cv_given_0 = float((data_psi_t["psi_t_0"].values * data_psi_T_0["psi_j"].values).sum())

    # ---- Bivariate normal PDFs for each (B_t, B_T) combination ----
    def make_bvn_pdf(mu: list[float], cov_matrix: np.ndarray):
        rv = multivariate_normal(mean=mu, cov=cov_matrix, allow_singular=True)
        return rv.pdf

    mu_11    = [mom_T["M_YT_11"], mom_T["M_Yt_11"]]
    sigma_11 = np.array([[mom_T["S_YT_11"] ** 2, cv_given_1],
                         [cv_given_1,             mom_T["S_Yt_11"] ** 2]])
    pdf_11   = make_bvn_pdf(mu_11, sigma_11)

    mu_01    = [mom_T["M_YT_01"], mom_T["M_Yt_01"]]
    sigma_01 = np.array([[mom_T["S_YT_01"] ** 2, cv_given_0],
                         [cv_given_0,             mom_T["S_Yt_01"] ** 2]])
    pdf_01   = make_bvn_pdf(mu_01, sigma_01)

    mu_10    = [mom_T["M_YT_10"], mom_T["M_Yt_10"]]
    sigma_10 = np.array([[mom_T["S_YT_10"] ** 2, cv_given_1],
                         [cv_given_1,             mom_T["S_Yt_10"] ** 2]])
    pdf_10   = make_bvn_pdf(mu_10, sigma_10)

    mu_00    = [mom_T["M_YT_00"], mom_T["M_Yt_00"]]
    sigma_00 = np.array([[mom_T["S_YT_00"] ** 2, cv_given_0],
                         [cv_given_0,             mom_T["S_Yt_00"] ** 2]])
    pdf_00   = make_bvn_pdf(mu_00, sigma_00)

    # ---- Marginal mixture PDFs ----
    def pdf_Yt(x):
        return dmixnorm(
            x,
            means=[mom_t["M_Y1"], mom_t["M_Y0"]],
            sds=[mom_t["S_Y1"], mom_t["S_Y0"]],
            probs=[mom_t["p1"], 1 - mom_t["p1"]],
        )

    def pdf_YT(x):
        return dmixnorm(
            x,
            means=[mom_T["M_Y1"], mom_T["M_Y0"]],
            sds=[mom_T["S_Y1"], mom_T["S_Y0"]],
            probs=[mom_T["p1"], 1 - mom_T["p1"]],
        )

    p1_t = mom_t["p1"]
    p1_T = mom_T["p1"]

    def joint_pdf_YtT(x):
        """Joint mixture PDF of (Y_T, Y_t). x must be shape (..., 2)."""
        return (
              p1_T *       p1_t  * pdf_11(x)
            + p1_T * (1 - p1_t) * pdf_01(x)
            + (1 - p1_T) * p1_t  * pdf_10(x)
            + (1 - p1_T) * (1 - p1_t) * pdf_00(x)
        )

    mom_T["pdf"]    = [{"Yt": pdf_Yt, "YT": pdf_YT, "joint": joint_pdf_YtT}]
    mom_T["t_cond"] = mom_t["date"]

    return mom_T


# ---------------------------------------------------------------------------
# solar_model_predict
# ---------------------------------------------------------------------------

def solar_model_predict(
    model,
    moments: pd.Series,
    lambda_: float = 0.0,
    ci: float = 0.01,
) -> dict:
    """
    Produce a single-date solar radiation forecast.

    Parameters
    ----------
    model   : SolarModel with attributes ``transform``, ``data``, ``spec``
    moments : pd.Series, one row from model.moments.conditional
    lambda_ : float, Sugeno parameter (0 = no distortion)
    ci      : float, confidence level for prediction intervals

    Returns
    -------
    dict with keys: ``grid`` (pd.DataFrame), ``df_n`` (pd.Series), ``ci``
    and class attribute ``solarModelForecast``.
    """
    df_n = moments.copy()

    comb = pd.DataFrame({
        "mean":  [df_n["M_Y1"], df_n["M_Y0"]],
        "sd":    [df_n["S_Y1"], df_n["S_Y0"]],
        "probs": [df_n["p1"],   1 - df_n["p1"]],
    })

    if lambda_ == 0:
        def pdf_Yt(x):
            return dmixnorm(x, comb["mean"], comb["sd"], comb["probs"])
        def pdf_Yt_up(x):
            return dmixnorm(x, [comb["mean"][0]], [comb["sd"][0]], [1.0])
        def pdf_Yt_dw(x):
            return dmixnorm(x, [comb["mean"][1]], [comb["sd"][1]], [1.0])
        def cdf_Yt(x):
            return pmixnorm(x, comb["mean"], comb["sd"], comb["probs"])
        def cdf_Yt_up(x):
            return pmixnorm(x, [comb["mean"][0]], [comb["sd"][0]], [1.0])
        def cdf_Yt_dw(x):
            return pmixnorm(x, [comb["mean"][1]], [comb["sd"][1]], [1.0])
    else:
        _pdf_Yt    = lambda x: dmixnorm(x, comb["mean"], comb["sd"], comb["probs"])
        _pdf_Yt_up = lambda x: dmixnorm(x, [comb["mean"][0]], [comb["sd"][0]], [1.0])
        _pdf_Yt_dw = lambda x: dmixnorm(x, [comb["mean"][1]], [comb["sd"][1]], [1.0])
        _cdf_Yt    = lambda x: pmixnorm(x, comb["mean"], comb["sd"], comb["probs"])
        _cdf_Yt_up = lambda x: pmixnorm(x, [comb["mean"][0]], [comb["sd"][0]], [1.0])
        _cdf_Yt_dw = lambda x: pmixnorm(x, [comb["mean"][1]], [comb["sd"][1]], [1.0])
        # Sugeno distortions
        pdf_Yt    = dsugeno(_pdf_Yt,    _cdf_Yt,    lambda_)
        pdf_Yt_up = dsugeno(_pdf_Yt_up, _cdf_Yt_up, lambda_)
        pdf_Yt_dw = dsugeno(_pdf_Yt_dw, _cdf_Yt_dw, lambda_)
        cdf_Yt    = psugeno(_cdf_Yt,    lambda_)
        cdf_Yt_up = psugeno(_cdf_Yt_up, lambda_)
        cdf_Yt_dw = psugeno(_cdf_Yt_dw, lambda_)

    def e_Rt_q(q: int, pdf) -> float:
        """E[R_t^q] via numerical integration over Y-space."""
        val, _ = sci_integrate.quad(
            lambda x: model.transform.iRY(x, df_n["Ct"]) ** q * pdf(x),
            -np.inf, np.inf,
        )
        return val

    # ---- Expected values ----
    df_n["e_Rt"]    = e_Rt_q(1, pdf_Yt)
    df_n["e_Rt_up"] = e_Rt_q(1, pdf_Yt_up)
    df_n["e_Rt_dw"] = e_Rt_q(1, pdf_Yt_dw)

    # ---- Variances ----
    df_n["v_Rt"]    = e_Rt_q(2, pdf_Yt)    - df_n["e_Rt"]    ** 2
    df_n["v_Rt_up"] = e_Rt_q(2, pdf_Yt_up) - df_n["e_Rt_up"] ** 2
    df_n["v_Rt_dw"] = e_Rt_q(2, pdf_Yt_dw) - df_n["e_Rt_dw"] ** 2

    # ---- Confidence intervals ----
    link = model.spec["transform"]["link"]

    df_n["ci_mix_lo"] = q_solar_ghi(ci,     df_n["Ct"], model.transform.alpha, model.transform.beta, cdf_Yt,    link)
    df_n["ci_mix_hi"] = q_solar_ghi(1 - ci, df_n["Ct"], model.transform.alpha, model.transform.beta, cdf_Yt,    link)
    df_n["ci_up_lo"]  = q_solar_ghi(ci,     df_n["Ct"], model.transform.alpha, model.transform.beta, cdf_Yt_up, link)
    df_n["ci_up_hi"]  = q_solar_ghi(1 - ci, df_n["Ct"], model.transform.alpha, model.transform.beta, cdf_Yt_up, link)
    df_n["ci_dw_lo"]  = q_solar_ghi(ci,     df_n["Ct"], model.transform.alpha, model.transform.beta, cdf_Yt_dw, link)
    df_n["ci_dw_hi"]  = q_solar_ghi(1 - ci, df_n["Ct"], model.transform.alpha, model.transform.beta, cdf_Yt_dw, link)

    # ---- PDF grid ----
    n_points = 100
    lower_Rt = df_n["Ct"] * model.transform.bounds("Kt")[0]
    upper_Rt = df_n["Ct"] * model.transform.bounds("Kt")[1]
    grid_x   = np.linspace(lower_Rt, upper_Rt, n_points + 2)[1:-1]   # exclude endpoints
    grid     = pd.DataFrame({"x": grid_x})

    def pdf_Rt(x, pdf_y):
        return d_solar_ghi(x, df_n["Ct"], model.transform.alpha, model.transform.beta, pdf_y, link)

    grid["pdf_Rt_mix"]    = pdf_Rt(grid["x"], pdf_Yt)
    grid["pdf_Rt_mix_up"] = pdf_Rt(grid["x"], pdf_Yt_up) * df_n["p1"]
    grid["pdf_Rt_mix_dw"] = pdf_Rt(grid["x"], pdf_Yt_dw) * (1 - df_n["p1"])

    # ---- PDF values at key points ----
    df_n["pdf_e_Rt"]      = pdf_Rt(df_n["e_Rt"],    pdf_Yt)
    df_n["pdf_ci_mix_lo"] = pdf_Rt(df_n["ci_mix_lo"], pdf_Yt)
    df_n["pdf_ci_mix_hi"] = pdf_Rt(df_n["ci_mix_hi"], pdf_Yt)
    df_n["pdf_ci_up_lo"]  = pdf_Rt(df_n["ci_up_lo"],  pdf_Yt_up)
    df_n["pdf_ci_up_hi"]  = pdf_Rt(df_n["ci_up_hi"],  pdf_Yt_up)
    df_n["pdf_ci_dw_lo"]  = pdf_Rt(df_n["ci_dw_lo"],  pdf_Yt_dw)
    df_n["pdf_ci_dw_hi"]  = pdf_Rt(df_n["ci_dw_hi"],  pdf_Yt_dw)

    # ---- Realised GHI ----
    match = model.data[model.data["date"] == df_n["date"]]["GHI"]
    df_n["Rt"] = float(match.iloc[0]) if len(match) > 0 else np.nan

    result = {"grid": grid, "df_n": df_n, "ci": ci}
    result["__class__"] = "solarModelForecast"
    return result


# ---------------------------------------------------------------------------
# solar_model_forecast
# ---------------------------------------------------------------------------

def solar_model_forecast(
    model,
    moments: pd.DataFrame,
    ci: float = 0.1,
    lambda_: float = 0.0,
) -> pd.DataFrame:
    """
    Iterate ``solar_model_predict`` over all rows of ``moments``.

    Parameters
    ----------
    model   : SolarModel
    moments : pd.DataFrame, one row per forecast date
    ci      : float, confidence level
    lambda_ : float, Sugeno parameter

    Returns
    -------
    pd.DataFrame — stacked ``df_n`` rows for all dates.
    """
    rows = []
    for i in range(len(moments)):
        try:
            result = solar_model_predict(
                model,
                moments=moments.iloc[i],
                ci=ci,
                lambda_=lambda_,
            )
            rows.append(result["df_n"])
        except Exception as e:
            warnings.warn(f"solar_model_forecast: row {i} failed — {e}")
            continue

    return pd.DataFrame(rows).reset_index(drop=True)


# ---------------------------------------------------------------------------
# solar_model_match_params
# ---------------------------------------------------------------------------

def solar_model_match_params(
    vec_params: dict | pd.Series,
    params: dict,
) -> dict:
    """
    Match a flat named parameter vector into a nested parameter dict,
    updating values in-place by name.

    Parameters
    ----------
    vec_params : dict or pd.Series
        Flat mapping ``{param_name: value}`` — e.g. ``{"theta": 1, "alpha1": 10}``.
    params : dict
        Nested parameter structure, e.g. ``model.coefficients``.
        Values must be dict-like or pd.Series with named entries.

    Returns
    -------
    dict — updated copy of ``params``.

    Examples
    --------
    >>> vec_params = {"theta": 1, "alpha1": 10}
    >>> updated = solar_model_match_params(vec_params, model.coefficients)
    """
    import copy
    params = copy.deepcopy(params)

    # Build a lookup: {param_name: block_key} for every name in every block
    names_map: dict[str, str] = {}
    for block_key, block_val in params.items():
        if hasattr(block_val, "keys"):
            for name in block_val.keys():
                names_map[name] = block_key
        elif hasattr(block_val, "index"):       # pd.Series
            for name in block_val.index:
                names_map[name] = block_key

    if isinstance(vec_params, pd.Series):
        vec_params = vec_params.to_dict()

    for name, value in vec_params.items():
        if name not in names_map:
            warnings.warn(
                f'solar_model_match_params: parameter "{name}" not found '
                f"in the model specification. Ignored!"
            )
            continue
        block_key = names_map[name]
        block = params[block_key]
        if isinstance(block, pd.Series):
            block[name] = value
        else:
            block[name] = value

    return params