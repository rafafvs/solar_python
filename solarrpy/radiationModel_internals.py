import warnings
import numpy as np
import pandas as pd
from scipy.integrate import quad

from .zzz import number_of_day
from .seasonalModel import SeasonalModel

"""
Solar model helper functions — Python translation.

Functions:
    create_monthly_sequence       — monthly date grid between t_now and t_hor
    martingale_method_seasonal    — mean-reversion θ via martingale estimator
    reparam_seasonal_function     — reparametrise seasonal variance to continuous time
    integral_sigma_numeric        — numerical integral of √(σ̄(s) · e^{-θ(T-s)})
    integral_sigma2_formula       — closed-form integral of σ̄²(s) · e^{-2θ(T-s)}

Version: 1.0.0
"""

def shift(arr, n):
    arr_shifted = np.empty_like(arr, dtype=float)
    if n == 0:
        arr_shifted[:] = arr
    if n > 0:
        arr_shifted[:n] = np.nan
        arr_shifted[n:] = arr[:-n]
    else:
        arr_shifted[n:] = np.nan
        arr_shifted[:n] = arr[-n:]
    return arr_shifted

def create_monthly_sequence(t_now, t_hor, last_day=False):
    t_now = pd.Timestamp(t_now)
    t_hor = pd.Timestamp(t_hor)

    # First day of the month for t_now and t_hor
    t_start = t_now.replace(day=1)
    t_end   = t_hor.replace(day=1)

    # Monthly sequence of first-of-month dates — mirrors seq.Date(..., by="1 month")
    months = pd.date_range(t_start, t_end, freq="MS")
    
    # dates_now: t_now, then last day of previous month for each subsequent first-of-month
    dates_now = [t_now] + [m - pd.Timedelta(days=m.day) for m in months[1:]]
    
    # dates_hor: last day of each month for all-but-last, then t_hor
    dates_hor = [m + pd.offsets.MonthEnd(0) for m in months[:-1]] + [t_hor]

    # Number of days in each sub-interval — mirrors difftime(dates_hor, dates_now)
    n_of_day = np.array([(h - s).days for s, h in zip(dates_now, dates_hor)], dtype=float)

    # Build base dataframe — Year from dates_now, Month from dates_hor
    df = pd.DataFrame({
        "Year":  [d.year  for d in dates_now],
        "Month": [d.month for d in dates_hor],
        "n":     n_of_day,
    })
    
    df["N"] = df["n"] + number_of_day(t_now) + np.concatenate([[0], np.cumsum(n_of_day[1:])])

    # n = N - n  (start index of each sub-interval)
    df["n"] = df["N"] - df["n"]

    # tau = total days t_now→t_hor + n[0]  (constant column)
    df["tau"] = (t_hor - t_now).days + df["n"].iloc[0]

    # Optional last-day split
    if last_day:
        last_row = df.iloc[[-1]].copy()
        last_row["n"] = last_row["tau"] - 1          # mirrors mutate(n = tau - 1)
        df.iloc[-1, df.columns.get_loc("N")] -= 1    # trim last row's N by 1
        df = pd.concat([df, last_row], ignore_index=True)

    df.attrs["last_day"] = last_day
    return df[['Year', 'Month', 'n', 'N', 'tau']]

# ---------------------------------------------------------------------------
# martingale_method_seasonal
# ---------------------------------------------------------------------------
def martingale_method_seasonal(Yt, Yt_bar, e_mu = 0):
    """
    Martingale estimation for mean-reversion parameter theta.
    Uses SeasonalModel to fit a time-varying variance to the squared differences.
    """
    Yt = np.asarray(Yt)
    Yt_bar = np.asarray(Yt_bar)
    
    if np.isscalar(e_mu):
        e_mu = np.full_like(Yt, e_mu)
    else:
        e_mu = np.asarray(e_mu)
        
    n = len(Yt)
    
    # Quadratic variation (dYt2)
    # Equivalent to (lag(Yt, 1) - lag(Yt, 2))^2
    dYt2 = (shift(Yt, 1) - shift(Yt, 2)) ** 2
    
    # Construct a DataFrame to feed into SeasonalModel
    data_df = pd.DataFrame({
        'dYt2': dYt2,
        'n': np.arange(1, n + 1)
    })
    
    # Drop NaNs for the fit
    fit_df = data_df.dropna(subset=['dYt2']).copy()
    
    # Initialize and fit the seasonal variance model
    seasonal_variance = SeasonalModel()
    seasonal_variance.fit(
        data=fit_df, 
        target_col='dYt2', 
        time_col='n', 
        include_intercept=True
    )
    
    pred_df = pd.DataFrame({'n': np.arange(n)})
    sigma2_bar = seasonal_variance.predict(data=pred_df, time_col='n')
    
    # Martingale estimation components
    # lag(Yt_bar, 1) - lag(Yt, 1)
    Y_est_L1 = (shift(Yt_bar, 1) - shift(Yt, 1)) / sigma2_bar.values
    dYt = Yt - Yt_bar
    dYt_L1 = shift(Yt, 1) - shift(Yt_bar, 1) - shift(e_mu, 1)
    
    df = pd.DataFrame({
        "Y_est_L1": Y_est_L1,
        "dYt":      dYt,
        "dYt_L1":   dYt_L1,
    }).dropna()
    
    # Create mask to drop NaNs
    a_n = np.sum(df["Y_est_L1"] * df["dYt"]) / np.sum(df["Y_est_L1"] * df["dYt_L1"])
    
    return -np.log(a_n)

# ---------------------------------------------------------------------------
# reparam_seasonal_function
# ---------------------------------------------------------------------------
def reparam_seasonal_function(par, theta, omega=2 * np.pi / 365):
    """
    Reparametrization for seasonal variance. Maps OLS parameters to continuous-time parameters.
    """
    a0, a1, a2 = par[0], par[1], par[2]
    
    # Correction for long term variance
    c0_long = a0 * 2 * theta
    c1_long = a1 * 2 * theta - omega * a2
    c2_long = a2 * 2 * theta + omega * a1
    
    gamma0_long = c0_long / (2 * theta)
    gamma1_long = (c1_long * 2 * theta + c2_long * omega) / (4 * theta**2 + omega**2)
    gamma2_long = (c2_long * 2 * theta - c1_long * omega) / (4 * theta**2 + omega**2)

    # Correction for short term variance
    alpha = 1 - np.exp(-2 * theta) * np.cos(omega)
    beta = np.exp(-2 * theta) * np.sin(omega)
    detM = alpha**2 + beta**2
    
    c0 = (2 * theta * a0) / (1 - np.exp(-2 * theta))
    c1 = ((2 * theta * alpha + omega * beta) * a1 + (2 * theta * beta - omega * alpha) * a2) / detM
    c2 = ((omega * alpha - 2 * theta * beta) * a1 + (omega * beta + 2 * theta * alpha) * a2) / detM
    
    gamma0 = c0 / (2 * theta)
    gamma1 = (c1 * 2 * theta + c2 * omega) / (4 * theta**2 + omega**2)
    gamma2 = (c2 * 2 * theta - c1 * omega) / (4 * theta**2 + omega**2)
    
    return {
        "alpha": alpha,
        "beta": beta,
        "detM": detM,
        "a_": np.array([a0, a1, a2]),
        "c_": np.array([c0, c1, c2]),
        "c_long": np.array([c0_long, c1_long, c2_long]),
        "gamma": np.array([gamma0, gamma1, gamma2]),
        "gamma_long": np.array([gamma0_long, gamma1_long, gamma2_long])
    }

# ---------------------------------------------------------------------------
# integral_sigma_numeric
# ---------------------------------------------------------------------------
def integral_sigma_numeric(theta, par, omega=2 * np.pi / 365):
    """
    Returns a function that computes the numerical integral for the square root of sigma_s * exp(-theta).
    """
    c0, c1, c2 = par[0], par[1], par[2]
    
    def seasonal_function(tau):
        return c0 + c1 * np.sin(omega * tau) + c2 * np.cos(omega * tau)
        
    def integrand(tau, T_val):
        val = seasonal_function(tau) * np.exp(-2 * theta * (T_val - tau))
        # Protect against negative values inside the sqrt due to numerical drift
        return np.sqrt(np.maximum(val, 0))
        
    def integrator(t, s, T_):
        t = np.atleast_1d(t)
        s = np.atleast_1d(s)
        T_ = np.atleast_1d(T_)
        
        result = np.zeros(len(T_))
        for i in range(len(T_)):
            res, _ = quad(integrand, t[i], s[i], args=(T_[i],))
            result[i] = res
        return result
        
    return integrator

# ---------------------------------------------------------------------------
# integral_sigma2_formula
# ---------------------------------------------------------------------------
def integral_sigma2_formula(theta, par, omega=2 * np.pi / 365):
    """
    Returns a vectorized function computing the exact formula integral for sigma_s^2 * exp(-2theta).
    """
    g0, g1, g2 = par[0], par[1], par[2]
    
    def f0(t, s, T_):
        return np.exp(-2 * theta * (T_ - s)) - np.exp(-2 * theta * (T_ - t))
        
    def f1(t, s, T_):
        return np.exp(-2 * theta * (T_ - s)) * np.sin(omega * s) - np.exp(-2 * theta * (T_ - t)) * np.sin(omega * t)
        
    def f2(t, s, T_):
        return np.exp(-2 * theta * (T_ - s)) * np.cos(omega * s) - np.exp(-2 * theta * (T_ - t)) * np.cos(omega * t)

    def integrator(t, s, T_):
        return g0 * f0(t, s, T_) + g1 * f1(t, s, T_) + g2 * f2(t, s, T_)
        
    return integrator