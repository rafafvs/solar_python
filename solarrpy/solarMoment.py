import numpy as np
import pandas as pd
from typing import Optional, Callable, Tuple, Dict, List
from datetime import datetime, timedelta
from scipy import stats, integrate


class SolarMoment:
    """
    Compute the generic conditional moments of a SolarModel object.
    
    This class computes multi-step ahead forecasts and conditional moments
    for solar radiation models, including ARMA, GARCH, and mixture components.
    
    Example:
    --------
    # Create solar moment object
    sm = SolarMoment(model, t_now="2022-01-01", t_hor="2022-03-31")
    
    # Filter to compute moments
    sm.filter(theta=0)
    
    # Get moments
    moments = sm.data
    
    Version: 1.0.0
    """
    
    VERSION = "1.0.0"
    
    def __init__(self, model, t_now: str, t_hor: str):
        """
        Initialize a SolarMoment object.
        
        Parameters:
        -----------
        model : SolarModel
            Fitted solar model object
        t_now : str
            Current date (pricing date)
        t_hor : str
            Horizon date (forecast target)
        """
        # Parse dates
        t_now = pd.to_datetime(t_now)
        t_hor = pd.to_datetime(t_hor)
        self.date = t_hor
        
        # Maximum order of AR / GARCH
        lag_max = max(model.ARMA.order[0], model.ARMA.order[1],
                      model.GARCH.order[0], model.GARCH.order[1])
        
        # Filter data between (t_now - lag_max + 1) and t_hor
        data = model.data[
            (model.data['date'] >= (t_now - pd.Timedelta(days=lag_max-1))) &
            (model.data['date'] <= t_hor)
        ].copy()
        
        # Check if data above t_hor are known
        if data['date'].max() < t_hor:
            # Add unknown variables
            new_dates_range = pd.date_range(
                start=data['date'].max() + pd.Timedelta(days=1),
                end=t_hor,
                freq='D'
            )
            
            new_dates = pd.DataFrame({
                'date': new_dates_range,
                'Year': new_dates_range.year,
                'Month': new_dates_range.month,
                'Day': new_dates_range.day,
                'n': new_dates_range.dayofyear,
                'isTrain': False,
                'weights': 0.0
            })
            
            # Add predictions
            new_dates['Ct'] = model.seasonal_model_Ct.predict(new_dates['n'].values)
            new_dates['Yt_bar'] = model.seasonal_model_Yt.predict(new_dates['n'].values)
            new_dates['sigma_bar'] = np.sqrt(
                model.seasonal_variance.predict(new_dates['n'].values)
            )
            new_dates['GHI_bar'] = model.transform.iRY(
                new_dates['Yt_bar'].values,
                new_dates['Ct'].values
            )
            
            # Merge with monthly data
            new_dates = new_dates.merge(
                model.monthly_data[['Month', 'Yt_tilde_uncond', 'sigma_uncond',
                                   'mu1', 'mu2', 'sd1', 'sd2', 'p1', 'p2']],
                on='Month',
                how='left'
            )
            
            data = pd.concat([data, new_dates], ignore_index=True)
        
        # Forecast seasonal variance
        data['sigma_bar'] = data['sigma_bar'] * data['sigma_uncond']
        
        # Store bounds parameters
        Ct = data['Ct'].iloc[-1]
        alpha = model.transform.alpha
        beta = model.transform.beta
        
        self._bounds = {'Ct': Ct, 'alpha': alpha, 'beta': beta}
        self._R_min_max = {
            'R_min': Ct * (1 - alpha - beta),
            'R_max': Ct * (1 - alpha)
        }
        self._link = model.transform.link
        
        # ========================================================================
        # Initialize derived columns
        # ========================================================================
        weights = data.iloc[lag_max:][
            ['date', 'Month', 'Yt_tilde_hat', 'sigma_bar', 'Yt_bar',
             'Yt_tilde_uncond', 'mu1', 'mu2', 'sd1', 'sd2', 'p1', 'p2']
        ].copy().reset_index(drop=True)
        
        # Step ahead
        weights['step'] = np.arange(1, len(weights) + 1)
        
        # ARMA intercept
        weights['ARMA_intercept'] = 0.0
        
        # ARMA weights
        weights['psi_j'] = 0.0
        weights['psi2_j'] = 0.0
        
        # GARCH moments
        weights['e_sigma1'] = np.nan
        weights['e_sigma2'] = np.nan
        weights['e_sigma4'] = np.nan
        
        # GARCH variances
        weights['v_sigma'] = 0.0
        weights['v_sigma2'] = 0.0
        
        # GARCH covariances
        weights['psi_hs'] = 0.0
        
        # Mixture moments
        weights['e_u'] = np.nan
        weights['e_u2'] = np.nan
        weights['e_u4'] = np.nan
        
        # Mixture variance
        weights['v_u'] = np.nan
        
        # Cumulated moments of U_{t+h}
        weights['e_Uh'] = 0.0
        weights['v_Uh'] = 0.0
        
        # Expectation of Yt (unconditional)
        weights['e_Yt'] = 0.0
        
        # Variance of Yt (unconditional)
        weights['Sigma2'] = 1.0
        
        self.weights = weights
        
        # ========================================================================
        # Known data at time t_now (used as vector in state-space forecast)
        # ========================================================================
        df_t = data.iloc[:lag_max].copy()
        df_t = df_t.sort_values('date', ascending=False)
        
        # State vector
        AR_order = model.ARMA.order[0]
        Y0 = df_t['Yt_tilde'].iloc[:AR_order].values
        
        eps0 = np.array([])
        MA_order = model.ARMA.order[1]
        if MA_order > 0:
            eps0 = df_t['eps'].iloc[:MA_order].values
        
        # Store state vector
        self.X0 = np.concatenate([Y0, eps0])
        
        # Store GARCH initial conditions
        self.H0 = np.array([
            data['sigma'].iloc[lag_max-1]**2,
            data['eps_tilde'].iloc[lag_max-1]
        ])
        
        # ========================================================================
        # Initialize moments data
        # ========================================================================
        df_T = data[data['date'] == t_hor].iloc[0]
        
        # Number of steps ahead
        h = len(weights)
        
        # Final object
        self._data = pd.DataFrame({
            'date': [t_hor],
            'lag_max': [lag_max],
            'h': [h],
            'Year': [t_hor.year],
            'Month': [t_hor.month],
            'Day': [t_hor.day],
            'e_Yt': [0.0],
            'sd_Yt': [1.0],
            'M_Y1': [0.0],
            'S_Y1': [1.0],
            'M_Y0': [1.0],
            'S_Y0': [0.0],
            'p1': [df_T['p1']],
            'GHI_bar': [df_T['GHI_bar']]
        })
        
        # ========================================================================
        # Store ARMA matrices and vectors
        # ========================================================================
        self._ARMA = {
            'intercept': model.ARMA.intercept if hasattr(model.ARMA, 'intercept') else 0,
            'A': model.ARMA.A,
            'b': model.ARMA.b
        }
        
        # Store GARCH parameters
        e1 = np.zeros(len(model.GARCH.d))
        e1[0] = 1.0
        
        self._GARCH = {
            'A': model.GARCH.A,
            'b': model.GARCH.b,
            'd': model.GARCH.d,
            'beta': model.GARCH.beta,
            'e1': e1
        }
    
    def filter(self, theta: float = 0, B: Optional[np.ndarray] = None,
               t_cond: Optional[np.ndarray] = None):
        """
        Filter to compute all moments.
        
        Parameters:
        -----------
        theta : float
            Shift parameter for mixture
        B : np.ndarray, optional
            Conditioning Bernoulli variable
        t_cond : np.ndarray, optional
            Conditioning dates
        """
        self.filter_ARMA()
        self.filter_NM(theta, B, t_cond)
        self.filter_GARCH()
        self.filter_weights()
    
    def filter_ARMA(self):
        """Filter ARMA component to compute psi weights."""
        h = len(self.weights)
        
        # ARMA forecast
        df_tT = ARMA_forecast(
            h, self.X0,
            self._ARMA['A'],
            self._ARMA['b'],
            intercept=self._ARMA['intercept']
        )['weights'][0]
        
        # Store ARMA summations
        self.weights['psi_j'] = df_tT['psi_j']
        self.weights['psi2_j'] = df_tT['psi2_j']
        self.weights['Yt_tilde_hat'] = df_tT['Yt_tilde_hat']
    
    def filter_NM(self, theta: float = 0, B: Optional[np.ndarray] = None,
                  t_cond: Optional[np.ndarray] = None):
        """
        Filter Normal Mixture component.
        
        Parameters:
        -----------
        theta : float
            Shift parameter for mixture means
        B : np.ndarray, optional
            Conditioning Bernoulli variable
        t_cond : np.ndarray, optional
            Conditioning dates
        """
        df_tT = self.weights.copy()
        
        # ====================================================================
        # 1) Conditioning variable
        # ====================================================================
        B1 = np.ones(len(df_tT))
        B0 = np.ones(len(df_tT))
        
        if t_cond is not None and B is not None:
            if len(t_cond) == len(B):
                idx_cond = df_tT['date'].isin(t_cond)
                B1[idx_cond] = B / df_tT.loc[idx_cond, 'p1']
                B0[idx_cond] = (1 - B) / df_tT.loc[idx_cond, 'p2']
        
        # ====================================================================
        # 2) Mixture moments
        # ====================================================================
        if theta != 0:
            df_tT['mu1'] = df_tT['mu1'] + df_tT['sd1'] * theta
            df_tT['mu2'] = df_tT['mu2'] + df_tT['sd2'] * theta
        
        # Compute mixture moments
        # First moment
        df_tT['e_u'] = (df_tT['mu1'] * df_tT['p1']) * B1 + \
                       (df_tT['mu2'] * df_tT['p2']) * B0
        
        # Second moment
        df_tT['e_u2'] = ((df_tT['mu1']**2 + df_tT['sd1']**2) * df_tT['p1']) * B1 + \
                        ((df_tT['mu2']**2 + df_tT['sd2']**2) * df_tT['p2']) * B0
        
        # Fourth moment
        df_tT['e_u4'] = ((3 * df_tT['sd1']**4 + 6 * df_tT['mu1']**2 * df_tT['sd1']**2 +
                         df_tT['mu1']**4) * df_tT['p1']) * B1 + \
                        ((3 * df_tT['sd2']**4 + 6 * df_tT['mu2']**2 * df_tT['sd2']**2 +
                         df_tT['mu2']**4) * df_tT['p2']) * B0
        
        # Variance
        df_tT['v_u'] = df_tT['e_u2'] - df_tT['e_u']**2
        
        self.weights = df_tT
    
    def filter_GARCH(self):
        """Filter GARCH component to compute variance moments."""
        h = len(self.weights)
        
        # Forecast GARCH moments
        df_tT = sGARCH_forecast(
            h,
            self._GARCH['A'],
            self._GARCH['b'],
            self._GARCH['d'],
            self._GARCH['e1'],
            self.H0,
            self.weights['e_u2'].values,
            self.weights['e_u4'].values
        )
        
        # Store GARCH moments (exact)
        self.weights['e_sigma1'] = df_tT['e_sigma']
        self.weights['e_sigma2'] = df_tT['e_sigma2']
        self.weights['e_sigma4'] = df_tT['e_sigma4']
        self.weights['v_sigma2'] = df_tT['v_sigma2']
        self.weights['v_sigma'] = df_tT['v_sigma']
        
        psi_j = self.weights['psi_j'].values
        
        # Compute covariances
        for j in range(h):
            if j < h - 1:
                cv_sigma_j = df_tT['cv_sigma'][j]
                psi_hs_j = np.sum(cv_sigma_j * psi_j[j] * np.delete(psi_j, j))
                self.weights.loc[j, 'psi_hs'] = psi_hs_j
    
    def filter_weights(self):
        """Filter and compute final conditional moments."""
        df_tT = self.weights.copy()
        h = len(df_tT)
        
        # ====================================================================
        # 3) Compute the series of psi
        # ====================================================================
        # Compute weights for expectations
        df_tT['e_Uh'] = (df_tT['psi_j'] * df_tT['sigma_bar'] *
                         df_tT['e_sigma1'] * df_tT['e_u'])
        
        # Variance
        df_tT['v_Uh'] = (df_tT['psi2_j'] * df_tT['sigma_bar']**2 *
                        (df_tT['e_sigma2'] * df_tT['e_u2'] -
                         df_tT['e_sigma1']**2 * df_tT['e_u']**2))
        
        # Forecasted total expectation Yt_tilde
        df_tT['e_Yt_tilde'] = df_tT['Yt_tilde_hat'] + df_tT['e_Uh'].cumsum()
        
        # Forecasted total expectation Yt
        df_tT['e_Yt'] = (df_tT['e_Yt_tilde'] + df_tT['Yt_tilde_uncond'] +
                         df_tT['Yt_bar'])
        
        # Forecasted total variance of Yt
        df_tT['Sigma2'] = df_tT['v_Uh'].cumsum()
        
        # ====================================================================
        # 4) Approximate the multinomial mixture with a 2-component GM
        # ====================================================================
        # Last value
        df_T = df_tT.iloc[-1]
        
        # Next step parameters
        # Approximated mixture means
        M_Y = (df_T['Yt_tilde_hat'] + df_T['Yt_tilde_uncond'] + df_T['Yt_bar'] +
               df_T['sigma_bar'] * df_T['e_sigma1'] * np.array([df_T['mu1'], df_T['mu2']]))
        
        # Approximated mixture variances
        S2_Y = (np.array([df_T['sd1']**2, df_T['sd2']**2]) *
                df_T['sigma_bar']**2 * df_T['e_sigma2'] +
                df_T['v_sigma'] * df_T['sigma_bar']**2 *
                np.array([df_T['mu1']**2, df_T['mu2']**2]))
        
        if h > 1:
            # Add conditional covariances
            rho2_U = (np.sum(df_tT['psi_hs'].iloc[:-1] * df_tT['e_u'].iloc[:-1]) +
                     df_T['psi_hs'] * np.array([df_T['mu1'], df_T['mu2']]))
            
            # Approximated mixture means
            M_Y = M_Y + np.sum(df_tT['e_Uh'].iloc[:-1])
            
            # Approximated mixture variances
            S2_Y = S2_Y + np.sum(df_tT['v_Uh'].iloc[:-1]) + 2 * rho2_U
            
            # Update total variance
            df_tT['Sigma2'] = df_tT['Sigma2'] + (df_tT['psi_hs'] * df_tT['e_u']).cumsum()
        
        S_Y = np.sqrt(S2_Y)
        
        # Update total moments
        self._data['e_Yt'] = df_tT['e_Yt'].iloc[-1]
        
        if h == 1:
            self._data['sd_Yt'] = df_T['sigma_bar'] * df_T['e_sigma1']
        else:
            self._data['sd_Yt'] = np.sqrt(df_tT['Sigma2'].iloc[-1])
        
        # Update mixture moments
        self._data['M_Y1'] = M_Y[0]
        self._data['M_Y0'] = M_Y[1]
        self._data['S_Y1'] = S_Y[0]
        self._data['S_Y0'] = S_Y[1]
        
        # Update weights
        self.weights = df_tT
    
    def pdf_Y(self, type: str = "mix") -> Callable:
        """
        Get PDF function for transformed variable Y.
        
        Parameters:
        -----------
        type : str
            Type: 'mix', 'up', or 'dw'
        
        Returns:
        --------
        Callable
            PDF function
        """
        mom = self._data.iloc[0]
        
        if type == "mix":
            means = np.array([mom['M_Y1'], mom['M_Y0']])
            sd = np.array([mom['S_Y1'], mom['S_Y0']])
            probs = np.array([mom['p1'], 1 - mom['p1']])
            return lambda x: dmixnorm(x, means, sd, probs)
        
        elif type == "up":
            return lambda x: stats.norm.pdf(x, mom['M_Y1'], mom['S_Y1'])
        
        elif type == "dw":
            return lambda x: stats.norm.pdf(x, mom['M_Y0'], mom['S_Y0'])
        
        else:
            raise ValueError("type must be 'mix', 'up', or 'dw'")
    
    def cdf_Y(self, type: str = "mix") -> Callable:
        """
        Get CDF function for transformed variable Y.
        
        Parameters:
        -----------
        type : str
            Type: 'mix', 'up', or 'dw'
        
        Returns:
        --------
        Callable
            CDF function
        """
        mom = self._data.iloc[0]
        
        if type == "mix":
            means = np.array([mom['M_Y1'], mom['M_Y0']])
            sd = np.array([mom['S_Y1'], mom['S_Y0']])
            probs = np.array([mom['p1'], 1 - mom['p1']])
            return lambda x: pmixnorm(x, means, sd, probs)
        
        elif type == "up":
            return lambda x: stats.norm.cdf(x, mom['M_Y1'], mom['S_Y1'])
        
        elif type == "dw":
            return lambda x: stats.norm.cdf(x, mom['M_Y0'], mom['S_Y0'])
        
        else:
            raise ValueError("type must be 'mix', 'up', or 'dw'")
    
    def pdf_R(self, type: str = "mix") -> Callable:
        """
        Get PDF function for solar radiation R.
        
        Parameters:
        -----------
        type : str
            Type: 'mix', 'up', or 'dw'
        
        Returns:
        --------
        Callable
            PDF function for R
        """
        pdf_Y = self.pdf_Y(type)
        Ct = self._bounds['Ct']
        alpha = self._bounds['alpha']
        beta = self._bounds['beta']
        link = self._link
        
        return lambda x: dsolarGHI(x, Ct, alpha, beta, pdf_Y, link=link)
    
    def cdf_R(self, type: str = "mix") -> Callable:
        """
        Get CDF function for solar radiation R.
        
        Parameters:
        -----------
        type : str
            Type: 'mix', 'up', or 'dw'
        
        Returns:
        --------
        Callable
            CDF function for R
        """
        cdf_Y = self.cdf_Y(type)
        Ct = self._bounds['Ct']
        alpha = self._bounds['alpha']
        beta = self._bounds['beta']
        link = self._link
        
        return lambda x: psolarGHI(x, Ct, alpha, beta, cdf_Y, link=link)
    
    def Q_R(self, type: str = "mix") -> Callable:
        """
        Get quantile function for solar radiation R.
        
        Parameters:
        -----------
        type : str
            Type: 'mix', 'up', or 'dw'
        
        Returns:
        --------
        Callable
            Quantile function for R
        """
        cdf_Y = self.cdf_Y(type)
        Ct = self._bounds['Ct']
        alpha = self._bounds['alpha']
        beta = self._bounds['beta']
        link = self._link
        
        return lambda p: qsolarGHI(p, Ct, alpha, beta, cdf_Y, link=link)
    
    def __repr__(self) -> str:
        """Print method for SolarMoment class."""
        h = len(self.weights)
        t_hor = self._data['date'].iloc[0]
        t_now = t_hor - pd.Timedelta(days=h)
        
        lines = [
            "------------- Solar Moments -------------",
            f"Time to maturity: {h}",
            f"t_now: {t_now.date()}",
            f"t_hor: {t_hor.date()}"
        ]
        
        return "\n".join(lines)
    
    # ==================== Properties ====================
    
    @property
    def data(self) -> pd.DataFrame:
        """Get moments data."""
        return self._data
    
    @property
    def ARMA(self) -> Dict:
        """Get ARMA parameters."""
        return self._ARMA
    
    @property
    def GARCH(self) -> Dict:
        """Get GARCH parameters."""
        return self._GARCH
    
    @property
    def R_min_max(self) -> Dict:
        """Get radiation bounds."""
        return self._R_min_max


# Helper functions (would need to be implemented)
def ARMA_forecast(h, X0, A, b, intercept):
    """Forecast ARMA process."""
    pass

def sGARCH_forecast(h, A, b, d, e1, H0, e_u2, e_u4):
    """Forecast GARCH variance."""
    pass

def dmixnorm(x, means, sd, probs):
    """PDF of Gaussian mixture."""
    pass

def pmixnorm(x, means, sd, probs):
    """CDF of Gaussian mixture."""
    pass

def dsolarGHI(x, Ct, alpha, beta, pdf_Y, link):
    """PDF for solar GHI."""
    pass

def psolarGHI(x, Ct, alpha, beta, cdf_Y, link):
    """CDF for solar GHI."""
    pass

def qsolarGHI(p, Ct, alpha, beta, cdf_Y, link):
    """Quantile function for solar GHI."""
    pass
