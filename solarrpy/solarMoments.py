import numpy as np
import pandas as pd
from typing import Optional, List, Dict, Callable, Union
from scipy import optimize
import warnings


def solarMoments_conditional(data: pd.DataFrame, theta: float = 0,
                             control_model) -> pd.DataFrame:
    """
    Compute the conditional moments.
    
    Parameters:
    -----------
    data : pd.DataFrame
        Data slot of a SolarModel object
    theta : float
        Shift parameter for Gaussian mixture residuals
    control_model : SolarModelSpec
        Solar model specification object
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with conditional moments
    
    Version: 1.0.0
    """
    data = data.copy()
    
    if not control_model.garch_variance:
        data['sigma'] = 1.0
    
    # Compute conditional moments
    data['theta'] = theta
    
    # Conditional expectation Yt
    data['e_Yt'] = data['Yt_bar'] + data['Yt_tilde_hat'] + data['Yt_tilde_uncond']
    
    # Conditional std. deviation Yt
    data['sd_Yt'] = data['sigma'] * data['sigma_bar'] * data['sigma_uncond']
    
    # Conditional moments Yt (state up)
    data['M_Y1'] = data['e_Yt'] + data['sd_Yt'] * (data['mu1'] + data['sd1'] * theta)
    data['S_Y1'] = data['sd_Yt'] * data['sd1']
    
    # Conditional moments Yt (state down)
    data['M_Y0'] = data['e_Yt'] + data['sd_Yt'] * (data['mu2'] + data['sd2'] * theta)
    data['S_Y0'] = data['sd_Yt'] * data['sd2']
    
    # Store only relevant variables
    return data[['date', 'Year', 'Month', 'Day', 'e_Yt', 'sd_Yt',
                 'M_Y1', 'S_Y1', 'M_Y0', 'S_Y0', 'theta']]


def solarMoments_unconditional(data: pd.DataFrame, ARMA, GARCH,
                               theta: float = 0) -> pd.DataFrame:
    """
    Compute the unconditional moments.
    
    Parameters:
    -----------
    data : pd.DataFrame
        Data slot of a SolarModel object
    ARMA : ARMAModel
        ARMA model object
    GARCH : sGARCH
        GARCH model object
    theta : float
        Shift parameter for Gaussian mixture residuals
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with unconditional moments
    
    Version: 1.0.0
    """
    data = data.copy()
    
    # ARMA variance
    arma_variance = np.sqrt(ARMA.variance(1000)[999])
    
    # GARCH long-term std. deviation
    GARCH_vol = GARCH.omega[0] / (1 - np.sum(GARCH.coefficients[1:]))
    
    # Compute unconditional moments
    data['theta'] = theta
    
    # Unconditional expectation Yt
    intercept = ARMA.intercept if hasattr(ARMA, 'intercept') else 0
    data['e_Yt'] = intercept + data['Yt_bar'] + data['Yt_tilde_uncond']
    
    # Unconditional std. deviation Yt
    data['sd_Yt'] = (data['sigma_bar'] * data['sigma_uncond'] * 
                     GARCH_vol * arma_variance)
    
    # Unconditional moments Yt (state up)
    data['M_Y1'] = data['e_Yt'] + data['sd_Yt'] * (data['mu1'] + data['sd1'] * theta)
    data['S_Y1'] = data['sd_Yt'] * data['sd1']
    
    # Unconditional moments Yt (state down)
    data['M_Y0'] = data['e_Yt'] + data['sd_Yt'] * (data['mu2'] + data['sd2'] * theta)
    data['S_Y0'] = data['sd_Yt'] * data['sd2']
    
    # Store only relevant variables
    return data[['date', 'Year', 'Month', 'Day', 'e_Yt', 'sd_Yt',
                 'M_Y1', 'S_Y1', 'M_Y0', 'S_Y0', 'theta']]


def solarMoments(t_now: str, t_hor: str, data: pd.DataFrame,
                ARMA, GARCH, NM_model, transform,
                theta: float = 0, quiet: bool = False) -> pd.DataFrame:
    """
    Compute the generic conditional moments of a SolarModel object.
    
    Parameters:
    -----------
    t_now : str
        Date for today
    t_hor : str
        Horizon date
    data : pd.DataFrame
        Data slot of a SolarModel object
    ARMA : ARMAModel
        ARMA model object
    GARCH : sGARCH
        GARCH model object
    NM_model : SolarMixture
        Normal mixture model object
    transform : SolarTransform
        Transform object
    theta : float
        Shift parameter for Gaussian mixture residuals
    quiet : bool
        Suppress info messages
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with computed moments
    
    Example:
    --------
    moments = solarMoments("2019-07-11", "2019-10-19", 
                          model.data, model.ARMA, model.GARCH,
                          model.NM_model, model.transform, theta=0)
    
    Version: 1.0.0
    """
    if not quiet:
        print(f"Forecast: {t_hor} given {t_now}")
    
    # Parse dates
    t_now = pd.to_datetime(t_now)
    t_hor = pd.to_datetime(t_hor)
    
    # Maximum order of AR / GARCH
    lag_max = max(ARMA.order[0], ARMA.order[1], 
                  GARCH.order[0], GARCH.order[1])
    
    # Filter data between (t_now - lag_max + 1) and t_hor
    data = data[
        (data['date'] >= (t_now - pd.Timedelta(days=lag_max-1))) &
        (data['date'] <= t_hor)
    ].copy()
    
    # Unknown data till maturity
    df_tT = data.iloc[lag_max:].copy().reset_index(drop=True)
    
    # Known data at time t_now
    df_t = data.iloc[:lag_max].copy()
    df_t = df_t.sort_values('date', ascending=False)
    
    # Forecasting horizon
    h = len(df_tT)
    
    # Store alpha and beta
    df_tT['alpha'] = transform.alpha
    df_tT['beta'] = transform.beta
    
    # ========================================================================
    # 0) Forecast mean and variance of Yt_tilde
    # ========================================================================
    # Companion matrix
    A = ARMA.A
    
    # Residuals vector for mean
    b = ARMA.b.reshape(-1, 1)
    
    # Residuals matrix for variance
    B = b @ b.T
    
    # Extract first component
    e1 = np.zeros((len(b), 1))
    e1[0] = 1
    
    # Conditioning values
    Y0 = df_t['Yt_tilde'].iloc[:ARMA.order[0]].values
    eps0 = np.array([])
    if ARMA.order[1] > 0:
        eps0 = df_t['eps'].iloc[:ARMA.order[1]].values
    
    # State vector
    Xt = np.concatenate([Y0, eps0]).reshape(-1, 1)
    
    # Initialize variables
    df_tT['psi_j'] = np.nan
    df_tT['psi2_j'] = np.nan
    df_tT['intercept'] = 0.0
    
    intercept = ARMA.intercept if hasattr(ARMA, 'intercept') else 0
    
    for j in range(h):
        # Pre-compute the powers
        A_pow_j = np.linalg.matrix_power(A, j + 1)
        A_pow_hj = np.linalg.matrix_power(A, h - j - 1)
        
        # Compute weights for expectations
        df_tT.loc[j, 'psi_j'] = float(e1.T @ A_pow_hj @ b)
        
        # Intercept contribution
        df_tT.loc[j, 'intercept'] = float(e1.T @ A_pow_j @ (e1 * intercept))
        
        # Variance
        df_tT.loc[j, 'psi2_j'] = float(e1.T @ A_pow_hj @ B @ A_pow_hj.T @ e1)
        
        # Forecasted value
        intercept_sum = df_tT['intercept'].iloc[:j+1].sum()
        df_tT.loc[j, 'Yt_tilde_hat'] = intercept + intercept_sum + float(e1.T @ A_pow_j @ Xt)
    
    # ========================================================================
    # 2) Forecast seasonal variance
    # ========================================================================
    df_tT['sigma_bar'] = df_tT['sigma_bar'] * df_tT['sigma_uncond']
    
    return solarMoments_path(df_tT, GARCH, NM_model, theta=theta)


def solarMoments_path(moments: pd.DataFrame, GARCH, NM_model,
                     theta: float = 0, B: Optional[np.ndarray] = None,
                     t_cond: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Condition the moments for a specific Bernoulli realization at time t_cond.
    
    Parameters:
    -----------
    moments : pd.DataFrame
        Moments dataframe or moments with psi_j
    GARCH : sGARCH
        GARCH model object
    NM_model : SolarMixture
        Normal mixture model
    theta : float
        Shift parameter
    B : np.ndarray, optional
        Conditioning value for Bernoulli state
    t_cond : List[str], optional
        Conditioning dates
    
    Returns:
    --------
    pd.DataFrame
        DataFrame with conditioned moments
    
    Version: 1.0.0
    """
    # Handle different input formats
    if isinstance(moments.iloc[0].get('psi_j'), list):
        df_tT = moments['psi_j'].iloc[0].copy()
    else:
        df_tT = moments.copy()
    
    # Forecasting horizon
    h = len(df_tT)
    
    # Add mixture parameters
    df_tT = df_tT.drop(['mu1', 'mu2', 'sd1', 'sd2', 'p1', 'p2'], 
                       axis=1, errors='ignore')
    df_tT = df_tT.merge(NM_model.coefficients, on='Month', how='left')
    
    # Apply theta shift
    if theta != 0:
        df_tT['mu1'] = df_tT['mu1'] + df_tT['sd1'] * theta
        df_tT['mu2'] = df_tT['mu2'] + df_tT['sd2'] * theta
    
    # Compute mixture moments
    df_tT['mean'] = (df_tT['mu1'] * df_tT['p1']) + (df_tT['mu2'] * df_tT['p2'])
    
    # Second moment
    df_tT['m2'] = ((df_tT['mu1']**2 + df_tT['sd1']**2) * df_tT['p1'] +
                   (df_tT['mu2']**2 + df_tT['sd2']**2) * df_tT['p2'])
    
    # Fourth moment
    df_tT['m4'] = ((3 * df_tT['sd1']**4 + 6 * df_tT['mu1']**2 * df_tT['sd1']**2 + 
                    df_tT['mu1']**4) * df_tT['p1'] +
                   (3 * df_tT['sd2']**4 + 6 * df_tT['mu2']**2 * df_tT['sd2']**2 + 
                    df_tT['mu2']**4) * df_tT['p2'])
    
    # Variance
    df_tT['variance'] = df_tT['m2'] - df_tT['mean']**2
    
    # ========================================================================
    # 1) Conditioning variable
    # ========================================================================
    if t_cond is not None and B is not None:
        t_cond = pd.to_datetime(t_cond)
        for t_idx, t_val in enumerate(t_cond):
            mask = df_tT['date'] == t_val
            
            # Second moment
            df_tT.loc[mask, 'm2'] = np.where(
                B[t_idx] == 1,
                df_tT.loc[mask, 'mu1']**2 + df_tT.loc[mask, 'sd1']**2,
                df_tT.loc[mask, 'mu2']**2 + df_tT.loc[mask, 'sd2']**2
            )
            
            # Fourth moment
            df_tT.loc[mask, 'm4'] = np.where(
                B[t_idx] == 1,
                (3 * df_tT.loc[mask, 'sd1']**4 + 
                 6 * df_tT.loc[mask, 'mu1']**2 * df_tT.loc[mask, 'sd1']**2 + 
                 df_tT.loc[mask, 'mu1']**4),
                (3 * df_tT.loc[mask, 'sd2']**4 + 
                 6 * df_tT.loc[mask, 'mu2']**2 * df_tT.loc[mask, 'sd2']**2 + 
                 df_tT.loc[mask, 'mu2']**4)
            )
            
            # Expectation
            df_tT.loc[mask, 'mean'] = np.where(
                B[t_idx] == 1,
                df_tT.loc[mask, 'mu1'],
                df_tT.loc[mask, 'mu2']
            )
            
            # Variance
            df_tT.loc[mask, 'variance'] = np.where(
                B[t_idx] == 1,
                df_tT.loc[mask, 'sd1']**2,
                df_tT.loc[mask, 'sd2']**2
            )
    
    # ========================================================================
    # 1) Forecast GARCH moments
    # ========================================================================
    if h == 1:
        df_tT['e_sigma2'] = df_tT['sigma']**2
        df_tT['e_sigma4'] = df_tT['sigma']**4
    else:
        # Second moment GARCH variance (exact)
        df_tT['e_sigma2'] = e_sigma2_h_mix(
            h - 1, GARCH.omega[0], GARCH.alpha, GARCH.beta,
            e_x2=df_tT['m2'].values, sigma2_0=df_tT['sigma'].iloc[0]**2
        )
        
        # Fourth moment GARCH variance (exact)
        df_tT['e_sigma4'] = e_sigma4_h_mix(
            h - 1, GARCH.omega[0], GARCH.alpha, GARCH.beta,
            e_x2=df_tT['m2'].values, e_x4=df_tT['m4'].values,
            sigma4_0=df_tT['sigma'].iloc[0]**4
        )
    
    # Variance of GARCH (exact)
    df_tT['v_sigma2'] = df_tT['e_sigma4'] - df_tT['e_sigma2']**2
    
    # First moment of GARCH std.dev (approximated)
    df_tT['e_sigma1'] = (df_tT['e_sigma2']**(1/2) - 
                         (1/8) * df_tT['v_sigma2'] / df_tT['e_sigma2']**(3/2))
    
    # Variance of GARCH std.dev
    df_tT['v_sigma'] = df_tT['v_sigma2'] / (4 * df_tT['e_sigma2'])
    
    # Compute covariances
    df_tT['cv_psi_ij'] = 0.0
    if h > 1:
        for i in range(1, len(df_tT)):
            cv_psi_ij = cov_sigma2_h_mix(
                i + 1, 1, GARCH.omega[0], GARCH.alpha, GARCH.beta,
                e_x2=df_tT['m2'].values, e_x4=df_tT['m4'].values,
                sigma2_0=df_tT['sigma'].iloc[0]**2
            )
            
            cv_psi_ij = cv_psi_ij / (4 * np.sqrt(
                df_tT['e_sigma2'].iloc[i] * df_tT['e_sigma2'].iloc[1:i+1].values
            ))
            
            cv_psi_ij = cv_psi_ij * (
                df_tT['psi_j'].iloc[i] * df_tT['psi_j'].iloc[1:i+1].values
            )
            
            df_tT.loc[i, 'cv_psi_ij'] = cv_psi_ij.sum()
    
    # ========================================================================
    # 3) Compute the series of psi
    # ========================================================================
    # Compute weights for expectations
    df_tT['psi_y'] = (df_tT['psi_j'] * df_tT['sigma_bar'] * 
                      df_tT['e_sigma1'] * df_tT['mean'])
    
    # Variance
    df_tT['psi2_y'] = (df_tT['psi2_j'] * df_tT['sigma_bar']**2 * 
                       (df_tT['e_sigma2'] * df_tT['m2'] - 
                        df_tT['e_sigma1']**2 * df_tT['mean']**2))
    
    # Forecasted mean value
    df_tT['e_Yt_tilde'] = df_tT['Yt_tilde_hat'] + df_tT['psi_y'].cumsum()
    
    # Mean value
    df_tT['e_Yt'] = (df_tT['e_Yt_tilde'] + df_tT['Yt_tilde_uncond'] + 
                     df_tT['Yt_bar'])
    
    # Forecasted total variance
    df_tT['Sigma2'] = df_tT['psi2_y'].cumsum()
    
    # Store data to compute covariances
    data_psi = df_tT[[
        'date', 'Month', 'e_sigma1', 'e_sigma2', 'e_sigma4', 'v_sigma',
        'v_sigma2', 'Yt_tilde_hat', 'Yt_tilde_uncond', 'Yt_bar',
        'psi_j', 'psi_y', 'psi2_j', 'psi2_y', 'mu1', 'mu2', 'sd1', 'sd2',
        'p1', 'p2', 'mean', 'm2', 'm4', 'variance', 'sigma_bar', 'sigma'
    ]].copy()
    
    # Last value
    df_T = df_tT.iloc[-1]
    
    # ========================================================================
    # 4) Approximate the multinomial mixture with a 2-component GM
    # ========================================================================
    # Approximated mixture means
    M_Y = np.array([
        df_T['Yt_tilde_hat'] + df_T['Yt_tilde_uncond'] + df_T['Yt_bar'] + 
        df_T['sigma_bar'] * df_T['e_sigma1'] * df_T['mu1'],
        df_T['Yt_tilde_hat'] + df_T['Yt_tilde_uncond'] + df_T['Yt_bar'] + 
        df_T['sigma_bar'] * df_T['e_sigma1'] * df_T['mu2']
    ])
    
    # Approximated mixture variances
    S2_Y = np.array([
        df_T['sd1']**2 * df_T['sigma_bar']**2 * df_T['e_sigma2'] + 
        df_T['v_sigma'] * df_T['sigma_bar']**2 * df_T['mu1']**2,
        df_T['sd2']**2 * df_T['sigma_bar']**2 * df_T['e_sigma2'] + 
        df_T['v_sigma'] * df_T['sigma_bar']**2 * df_T['mu2']**2
    ])
    
    if h > 1:
        # Add conditional covariances
        rho2_U = (
            np.sum((df_tT['cv_psi_ij'] * df_tT['mean']).iloc[:-1]) +
            df_T['cv_psi_ij'] * np.array([df_T['mu1'], df_T['mu2']])
        )
        
        # Update means
        M_Y += np.sum(df_tT['psi_y'].iloc[:-1])
        
        # Update variances
        S2_Y += np.sum(df_tT['psi2_y'].iloc[:-1]) + 2 * rho2_U
        
        # Update total variance
        df_tT['Sigma2'] = df_tT['Sigma2'] + (df_tT['cv_psi_ij'] * df_tT['mean']).cumsum()
    
    S_Y = np.sqrt(S2_Y)
    
    t_hor = df_T['date']
    
    # Get original moments data
    if 'Ct' in moments.columns:
        moments_T = moments[moments['date'] == t_hor].iloc[0]
    else:
        moments_T = df_T
    
    return pd.DataFrame({
        'date': [t_hor],
        'h': [h],
        'Year': [t_hor.year],
        'Month': [t_hor.month],
        'Day': [t_hor.day],
        'e_Yt': [df_T['e_Yt']],
        'sd_Yt': [df_T['sigma_bar'] * df_T['sigma'] if h == 1 else np.sqrt(df_T['Sigma2'])],
        'M_Y1': [M_Y[0]],
        'S_Y1': [S_Y[0]],
        'M_Y0': [M_Y[1]],
        'S_Y0': [S_Y[1]],
        'psi_j': [[data_psi]],
        'p1': [moments_T.get('p1', df_T['p1'])],
        'Ct': [moments_T.get('Ct', np.nan)],
        'theta': [theta if isinstance(theta, (int, float)) else theta[-1]],
        'GHI_bar': [moments_T.get('GHI_bar', np.nan)],
        'alpha': [moments_T.get('alpha', np.nan)],
        'beta': [moments_T.get('beta', np.nan)]
    })


# Helper functions (would need to be implemented)
def e_sigma2_h_mix(h, omega, alpha, beta, e_x2, sigma2_0):
    """Compute E[sigma^2_{t+h}] for mixture."""
    pass

def e_sigma4_h_mix(h, omega, alpha, beta, e_x2, e_x4, sigma4_0):
    """Compute E[sigma^4_{t+h}] for mixture."""
    pass

def cov_sigma2_h_mix(h1, h2, omega, alpha, beta, e_x2, e_x4, sigma2_0):
    """Compute Cov(sigma^2_{t+h1}, sigma^2_{t+h2}) for mixture."""
    pass
