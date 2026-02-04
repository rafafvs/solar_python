import numpy as np
import pandas as pd
from typing import Optional, List, Dict, Tuple, Callable
from datetime import datetime
from scipy import stats, integrate, optimize
import warnings
import copy


class SolarMixture:
    """
    Monthly Gaussian mixture with two components.
    
    This class manages 12 separate Gaussian mixture models (one per month)
    for modeling solar radiation residuals.
    
    Example:
    --------
    # Create solar mixture object
    sm = SolarMixture(components=2, maxit=5000, abstol=1e-8)
    
    # Fit the model
    sm.fit(x=residuals, date=dates, weights=weights)
    
    # Get parameters
    params = sm.coefficients
    std_errs = sm.std_errors
    
    Version: 1.0.0
    """
    
    VERSION = "1.0.0"
    
    def __init__(self, components: int = 2, maxit: int = 5000, 
                 maxrestarts: int = 500, abstol: float = 1e-8):
        """
        Initialize a SolarMixture object.
        
        Parameters:
        -----------
        components : int
            Number of mixture components
        maxit : int
            Maximum number of iterations
        maxrestarts : int
            Maximum number of restarts
        abstol : float
            Absolute tolerance for convergence
        """
        self.components = components
        self.maxit = maxit
        self.maxrestarts = maxrestarts
        self.abstol = abstol
        
        # Monthly parameter functions
        self.mu1 = None
        self.mu2 = None
        self.sd1 = None
        self.sd2 = None
        self.prob = None
        
        # Private attributes
        self._x = None
        self._w = None
        self._date = None
        self._date_month = {}
        self._model = {}
        self._use_empiric = False
    
    def fit(self, x: np.ndarray, date: np.ndarray, 
            weights: Optional[np.ndarray] = None,
            method: str = "mixtools",
            mu_target: np.ndarray = None,
            var_target: np.ndarray = None):
        """
        Fit the Gaussian mixture parameters for each month.
        
        Parameters:
        -----------
        x : np.ndarray
            Time series data
        date : np.ndarray
            Date vector
        weights : np.ndarray, optional
            Observation weights (0 excludes, others included with unit weight)
        method : str
            Package for fitting: 'mclust' or 'mixtools'
        mu_target : np.ndarray, optional
            Target means to match (length 12)
        var_target : np.ndarray, optional
            Target variances to match (length 12)
        """
        # Store data
        self._x = np.array(x)
        self._date = pd.to_datetime(date)
        
        # Set weights
        if weights is None:
            self._w = np.ones(len(x))
        else:
            self._w = np.where(weights == 0, 0, 1)
        
        # Default targets
        if mu_target is None:
            mu_target = np.full(12, np.nan)
        if var_target is None:
            var_target = np.full(12, np.nan)
        
        # Get data DataFrame
        data = self.data
        
        # Fit Gaussian Mixture for each month
        for m in range(1, 13):
            data_month = data[data['Month'] == m].copy()
            w = data_month['w'].values
            eps = data_month['x'].values
            
            # Store monthly dates
            self._date_month[m] = data_month['date'].values
            
            # Fit Gaussian Mixture model
            GM_model = GaussianMixture(
                maxit=self.maxit,
                maxrestarts=self.maxrestarts,
                abstol=self.abstol,
                components=self.components
            )
            
            GM_model.fit(
                eps, w,
                method=method,
                mu_target=mu_target[m-1],
                var_target=var_target[m-1]
            )
            
            # Store model
            self._model[m] = copy.deepcopy(GM_model)
        
        # Initialize monthly parameter functions
        coefs = self.coefficients
        self.mu1 = MonthlyParams(coefs['mu1'].values)
        self.mu2 = MonthlyParams(coefs['mu2'].values)
        self.sd1 = MonthlyParams(coefs['sd1'].values)
        self.sd2 = MonthlyParams(coefs['sd2'].values)
        self.prob = MonthlyParams(coefs['p1'].values)
    
    def update(self, means: Optional[np.ndarray] = None,
               sd: Optional[np.ndarray] = None,
               p: Optional[np.ndarray] = None):
        """
        Update means, standard deviations, and probabilities.
        
        Parameters:
        -----------
        means : np.ndarray, optional
            Matrix of mean parameters (12 x components)
        sd : np.ndarray, optional
            Matrix of std deviation parameters (12 x components)
        p : np.ndarray, optional
            Matrix of probability parameters (12 x components)
        """
        # Use current values if not provided
        if means is None:
            means = self.means
        if sd is None:
            sd = self.sd
        if p is None:
            p = self.p
        
        # Update parameters for each month
        for m in range(1, 13):
            self._model[m].update(
                means=means[m-1, :],
                sd=sd[m-1, :],
                p=p[m-1, :]
            )
        
        # Update monthly parameter functions
        coefs = self.coefficients
        self.mu1 = MonthlyParams(coefs['mu1'].values)
        self.mu2 = MonthlyParams(coefs['mu2'].values)
        self.sd1 = MonthlyParams(coefs['sd1'].values)
        self.sd2 = MonthlyParams(coefs['sd2'].values)
        self.prob = MonthlyParams(coefs['p1'].values)
    
    def update_logLik(self):
        """Apply update_logLik() method to all GaussianMixture models."""
        for m in range(1, 13):
            self._model[m].update_logLik()
    
    def update_empiric_parameters(self):
        """Apply update_empiric_parameters() to all GaussianMixture models."""
        for m in range(1, 13):
            self._model[m].update_empiric_parameters()
    
    def filter(self):
        """Apply filter() method to all GaussianMixture models."""
        for m in range(1, 13):
            self._model[m].filter()
    
    def Hessian(self):
        """Apply Hessian() method to all GaussianMixture models."""
        for m in range(1, 13):
            self._model[m].Hessian()
    
    def use_empiric_parameters(self):
        """Toggle between empiric and EM parameters."""
        for m in range(1, 13):
            self._model[m].use_empiric_parameters()
        self._use_empiric = self._model[1].use_empiric
    
    def logLik(self, x: np.ndarray, date: np.ndarray) -> float:
        """
        Compute log-likelihood.
        
        Parameters:
        -----------
        x : np.ndarray
            Data vector
        date : np.ndarray
            Date vector
        
        Returns:
        --------
        float
            Total log-likelihood
        """
        df = pd.DataFrame({
            'date': pd.to_datetime(date),
            'Month': pd.to_datetime(date).month,
            'x': x
        })
        
        log_lik = 0
        for m in range(1, 13):
            data_month = df[df['Month'] == m]
            if len(data_month) > 0:
                log_lik += self._model[m].logLik(data_month['x'].values)
        
        return log_lik
    
    def grades(self, x: np.ndarray, date: np.ndarray) -> np.ndarray:
        """
        Compute probability integral transform (grades).
        
        Parameters:
        -----------
        x : np.ndarray
            Data vector
        date : np.ndarray
            Date vector
        
        Returns:
        --------
        np.ndarray
            Grades (uniform if model is correct)
        """
        df = pd.DataFrame({
            'date': pd.to_datetime(date),
            'Month': pd.to_datetime(date).month,
            'x': x,
            'u': np.nan
        })
        
        for m in range(1, 13):
            idx_month = df['Month'] == m
            if not idx_month.any():
                continue
            
            # Extract parameters
            means = self._model[m].means
            sd = self._model[m].sd
            p = self._model[m].p
            
            # Monthly CDF
            def cdf_M(x_val):
                return pmixnorm(x_val, means, sd, p)
            
            # Compute grades
            df.loc[idx_month, 'u'] = [cdf_M(val) for val in x[idx_month]]
        
        return df['u'].values
    
    def VaR(self, date: np.ndarray, alpha: np.ndarray = np.array([0.05])) -> np.ndarray:
        """
        Compute Value at Risk at specified confidence levels.
        
        Parameters:
        -----------
        date : np.ndarray
            Date vector
        alpha : np.ndarray
            Confidence levels for VaR
        
        Returns:
        --------
        np.ndarray
            VaR matrix (n_dates x n_alphas)
        """
        df = pd.DataFrame({
            'date': pd.to_datetime(date),
            'Month': pd.to_datetime(date).month
        })
        
        # Initialize VaR matrix
        VaR_alpha = np.zeros((len(df), len(alpha)))
        alpha = np.sort(alpha)
        
        for m in range(1, 13):
            idx_month = df['Month'] == m
            if not idx_month.any():
                continue
            
            # Extract parameters
            means = self._model[m].means
            sd = self._model[m].sd
            p = self._model[m].p
            
            # Monthly quantile function
            def Q_M(q):
                return qmixnorm(q, means, sd, p)
            
            # Compute VaR for this month
            VaR_alpha_M = np.array([Q_M(a) for a in alpha])
            VaR_alpha[idx_month, :] = VaR_alpha_M
        
        return VaR_alpha
    
    def ES(self, date: np.ndarray, alpha: np.ndarray = np.array([0.05])) -> np.ndarray:
        """
        Compute Expected Shortfall at specified confidence levels.
        
        Parameters:
        -----------
        date : np.ndarray
            Date vector
        alpha : np.ndarray
            Confidence levels for ES
        
        Returns:
        --------
        np.ndarray
            ES matrix (n_dates x n_alphas)
        """
        df = pd.DataFrame({
            'date': pd.to_datetime(date),
            'Month': pd.to_datetime(date).month
        })
        
        # Initialize ES matrix
        ES_alpha = np.zeros((len(df), len(alpha)))
        
        for m in range(1, 13):
            idx_month = df['Month'] == m
            if not idx_month.any():
                continue
            
            # Extract parameters
            means = self._model[m].means
            sd = self._model[m].sd
            p = self._model[m].p
            
            # Density function
            def f_R(x_val):
                return dmixnorm(x_val, means, sd, p)
            
            # Monthly quantile function
            def Q_M(q):
                return qmixnorm(q, means, sd, p)
            
            VaR_vals = np.array([Q_M(a) for a in alpha])
            
            # Expected Shortfall
            def compute_ES(VaR_val, alpha_val):
                result = integrate.quad(
                    lambda x: x * f_R(x),
                    -np.inf, VaR_val,
                    limit=50
                )
                return result[0] / alpha_val
            
            ES_alpha_M = np.array([compute_ES(v, a) for v, a in zip(VaR_vals, alpha)])
            ES_alpha[idx_month, :] = ES_alpha_M
        
        return ES_alpha
    
    def __repr__(self) -> str:
        """Print method for SolarMixture class."""
        lines = []
        for m in range(1, 13):
            month_name = pd.Timestamp(2020, m, 1).strftime('%B')
            lines.append(self._model[m].__repr__(month_name))
        return "\n".join(lines)
    
    # ==================== Properties ====================
    
    @property
    def data(self) -> pd.DataFrame:
        """
        Get data tibble with date, Month, x, and weights.
        
        Returns:
        --------
        pd.DataFrame
            DataFrame with columns: date, Month, x, w
        """
        return pd.DataFrame({
            'date': self._date,
            'Month': self._date.month,
            'x': self._x,
            'w': self._w
        })
    
    @property
    def means(self) -> np.ndarray:
        """
        Matrix of means (months x components).
        
        Returns:
        --------
        np.ndarray
            Mean parameters for all months
        """
        M = len(self._model)
        means = np.zeros((M, self.components))
        
        for m in range(1, M + 1):
            means[m-1, :] = self._model[m].means
        
        return means
    
    @property
    def sd(self) -> np.ndarray:
        """
        Matrix of standard deviations (months x components).
        
        Returns:
        --------
        np.ndarray
            Std deviation parameters for all months
        """
        M = len(self._model)
        sd = np.zeros((M, self.components))
        
        for m in range(1, M + 1):
            sd[m-1, :] = self._model[m].sd
        
        return sd
    
    @property
    def p(self) -> np.ndarray:
        """
        Matrix of probabilities (months x components).
        
        Returns:
        --------
        np.ndarray
            Probability parameters for all months
        """
        M = len(self._model)
        p = np.zeros((M, self.components))
        
        for m in range(1, M + 1):
            p[m-1, :] = self._model[m].p
        
        return p
    
    @property
    def model(self) -> Dict[int, 'GaussianMixture']:
        """Get dictionary of GaussianMixture models (one per month)."""
        return self._model
    
    @property
    def use_empiric(self) -> bool:
        """Check if empiric parameters are currently used."""
        return self._use_empiric
    
    @property
    def loglik(self) -> float:
        """Total log-likelihood across all months."""
        return sum(self._model[m].loglik for m in range(1, 13))
    
    @property
    def fitted(self) -> pd.DataFrame:
        """
        Get fitted classifications.
        
        Returns:
        --------
        pd.DataFrame
            DataFrame with date, B (classification), uncertainty
        """
        df_list = []
        for m in range(1, 13):
            if m in self._date_month:
                fitted_m = self._model[m].fitted
                fitted_m['date'] = self._date_month[m]
                df_list.append(fitted_m[['date', 'B', 'uncertainty']])
        
        df_fitted = pd.concat(df_list, ignore_index=True)
        return df_fitted.sort_values('date').reset_index(drop=True)
    
    @property
    def moments(self) -> pd.DataFrame:
        """
        Get theoretical moments for all months.
        
        Returns:
        --------
        pd.DataFrame
            DataFrame with Month, mean, variance, skewness, kurtosis, nobs, loglik
        """
        moments_list = []
        for m in range(1, 13):
            mom = self._model[m].moments.copy()
            mom['Month'] = m
            mom['loglik'] = self._model[m].loglik
            moments_list.append(mom)
        
        return pd.concat(moments_list, ignore_index=True)
    
    @property
    def coefficients(self) -> pd.DataFrame:
        """
        Get fitted parameters for all months.
        
        Returns:
        --------
        pd.DataFrame
            DataFrame with Month and mixture parameters
        """
        coef_list = []
        for m in range(1, 13):
            coef = self._model[m].model.copy()
            coef['Month'] = m
            coef_list.append(coef)
        
        return pd.concat(coef_list, ignore_index=True)
    
    @property
    def std_errors(self) -> pd.DataFrame:
        """
        Get standard errors for all months.
        
        Returns:
        --------
        pd.DataFrame
            DataFrame with Month and std errors
        """
        std_err_list = []
        for m in range(1, 13):
            se = pd.DataFrame(self._model[m].std_errors, index=[0])
            se['Month'] = m
            std_err_list.append(se)
        
        return pd.concat(std_err_list, ignore_index=True)
    
    @property
    def summary(self) -> pd.DataFrame:
        """Get summary statistics for all months."""
        summary_list = []
        for m in range(1, 13):
            summ = self._model[m].summary.copy()
            summ['Month'] = m
            summary_list.append(summ)
        
        return pd.concat(summary_list, ignore_index=True)


class MonthlyParams:
    """
    Monthly parameter interpolation function.
    
    This class provides methods to interpolate parameters across months.
    """
    
    def __init__(self, values: np.ndarray):
        """
        Initialize with 12 monthly values.
        
        Parameters:
        -----------
        values : np.ndarray
            Array of 12 monthly parameter values
        """
        self.values = np.array(values)
        if len(self.values) != 12:
            raise ValueError("Must provide exactly 12 monthly values")
    
    def __call__(self, month: int) -> float:
        """
        Get parameter value for a given month.
        
        Parameters:
        -----------
        month : int
            Month number (1-12)
        
        Returns:
        --------
        float
            Parameter value for that month
        """
        if not 1 <= month <= 12:
            raise ValueError("Month must be between 1 and 12")
        return self.values[month - 1]
    
    def interpolate(self, dates: np.ndarray) -> np.ndarray:
        """
        Interpolate parameter values for given dates.
        
        Parameters:
        -----------
        dates : np.ndarray
            Array of dates
        
        Returns:
        --------
        np.ndarray
            Interpolated parameter values
        """
        dates = pd.to_datetime(dates)
        months = dates.month
        return self.values[months - 1]


def solarMixture_moments_match(coefficients: pd.DataFrame,
                               e_target: np.ndarray,
                               v_target: np.ndarray,
                               sk_target: np.ndarray,
                               kt_target: np.ndarray,
                               x: List[np.ndarray]) -> pd.DataFrame:
    """
    Correct moments to ensure moment matching.
    
    Parameters:
    -----------
    coefficients : pd.DataFrame
        Initial mixture parameters
    e_target : np.ndarray
        Target means (length 12)
    v_target : np.ndarray
        Target variances (length 12)
    sk_target : np.ndarray
        Target skewness (length 12)
    kt_target : np.ndarray
        Target kurtosis (length 12)
    x : List[np.ndarray]
        List of 12 monthly data arrays
    
    Returns:
    --------
    pd.DataFrame
        Optimized coefficients
    """
    def loss_function(params, prob, e_t, v_t, sk_t, kt_t, x_data):
        """Loss function for optimization."""
        # Extract parameters
        n_comp = len(params) // 2
        means = params[:n_comp]
        sigma = params[n_comp:]
        probs = np.array([prob, 1 - prob])
        
        # Compute moments
        mom = GM_moments(means, sigma, probs)
        
        loss = 0
        # Add losses for target moments
        if not np.isnan(e_t):
            loss += abs(mom['mean'] - e_t)
        if not np.isnan(v_t):
            loss += abs(mom['variance'] - v_t)
        if not np.isnan(sk_t):
            loss += abs(mom['skewness'] - sk_t)
        if not np.isnan(kt_t):
            loss += abs(mom['kurtosis'] - kt_t)
        
        # Penalize heavily
        loss = 10000 * loss ** 2
        
        return loss
    
    # Optimize for each month
    opt_coefficients = coefficients.copy()
    
    for nmonth in range(len(coefficients)):
        # Initial parameters
        row = coefficients.iloc[nmonth]
        params_init = np.array([row['mu1'], row['mu2'], row['sd1'], row['sd2']])
        prob = row['p1']
        
        # Optimize
        result = optimize.minimize(
            loss_function,
            params_init,
            args=(prob, e_target[nmonth], v_target[nmonth],
                  sk_target[nmonth], kt_target[nmonth], x[nmonth]),
            method='Nelder-Mead'
        )
        
        # Update coefficients
        opt_coefficients.loc[nmonth, ['mu1', 'mu2', 'sd1', 'sd2']] = result.x
    
    # Ensure probabilities sum to 1
    opt_coefficients['p2'] = 1 - opt_coefficients['p1']
    
    return opt_coefficients


def solarModel_mvmixture(model_Ct, model_GHI):
    """
    Fit multivariate Gaussian mixture for clear sky and GHI models.
    
    Parameters:
    -----------
    model_Ct : SolarModel
        Solar model for clear sky
    model_GHI : SolarModel
        Solar model for GHI
    
    Returns:
    --------
    dict
        Dictionary with updated models
    """
    # Extract bivariate dataset
    data_GHI = model_GHI.data[['date', 'Month', 'u_tilde']].copy()
    data_GHI.columns = ['date', 'Month', 'GHI']
    
    data_Ct = model_Ct.data[['date', 'Month', 'u_tilde']].copy()
    data_Ct.columns = ['date', 'Month', 'Ct']
    
    # Merge datasets
    data = data_GHI.merge(data_Ct, on=['date', 'Month'], how='inner')
    
    # Remove outliers
    outliers_date_Ct = model_Ct._outliers.get('date', [])
    outliers_date_GHI = model_GHI._outliers.get('date', [])
    outliers_date = np.concatenate([outliers_date_Ct, outliers_date_GHI])
    data = data[~data['date'].isin(outliers_date)]
    
    # Extract Gaussian mixture parameters
    NM_model_GHI = copy.deepcopy(model_GHI._NM_model)
    NM_model_GHI.rho_up = np.zeros(12)
    NM_model_GHI.rho_dw = np.zeros(12)
    
    NM_model_Ct = copy.deepcopy(model_Ct._NM_model)
    NM_model_Ct.rho_up = np.zeros(12)
    NM_model_Ct.rho_dw = np.zeros(12)
    
    # Fit for each month
    data_months = []
    for m in range(1, 13):
        data_m = data[data['Month'] == m].copy()
        
        if len(data_m) == 0:
            continue
        
        # Monthly data (bivariate)
        eps = data_m[['GHI', 'Ct']].values
        
        # Fit multivariate Gaussian mixture
        gm = mvgaussian_mixture(eps, components=2)
        
        # Update Gaussian mixture parameters (GHI)
        NM_model_GHI.mu_up[m-1] = gm['params']['means'][0, 1]
        NM_model_GHI.mu_dw[m-1] = gm['params']['means'][1, 1]
        NM_model_GHI.sd_up[m-1] = np.sqrt(gm['params']['sigma2'][0, 1])
        NM_model_GHI.sd_dw[m-1] = np.sqrt(gm['params']['sigma2'][1, 1])
        NM_model_GHI.p_up[m-1] = gm['params']['p'][0]
        NM_model_GHI.p_dw[m-1] = gm['params']['p'][1]
        NM_model_GHI.rho_up[m-1] = gm['params']['rho'][0]
        NM_model_GHI.rho_dw[m-1] = gm['params']['rho'][1]
        
        # Update Gaussian mixture parameters (Ct)
        NM_model_Ct.mu_up[m-1] = gm['params']['means'][0, 0]
        NM_model_Ct.mu_dw[m-1] = gm['params']['means'][1, 0]
        NM_model_Ct.sd_up[m-1] = np.sqrt(gm['params']['sigma2'][0, 0])
        NM_model_Ct.sd_dw[m-1] = np.sqrt(gm['params']['sigma2'][1, 0])
        NM_model_Ct.p_up[m-1] = gm['params']['p'][0]
        NM_model_Ct.p_dw[m-1] = gm['params']['p'][1]
        NM_model_Ct.rho_up[m-1] = gm['params']['rho'][0]
        NM_model_Ct.rho_dw[m-1] = gm['params']['rho'][1]
        
        # Store Bernoulli classification
        data_m['B'] = gm['B_hat']['B1']
        data_months.append(data_m[['date', 'B']])
    
    # Combine fitted Bernoulli series
    data_months = pd.concat(data_months, ignore_index=True)
    
    # Update models
    model_Ct._data = model_Ct._data.drop('B', axis=1, errors='ignore')
    model_Ct._data = model_Ct._data.merge(data_months, on='date', how='left')
    model_Ct._data['B'] = model_Ct._data['B'].fillna(0)
    model_Ct._NM_model = NM_model_Ct
    
    model_GHI._data = model_GHI._data.drop('B', axis=1, errors='ignore')
    model_GHI._data = model_GHI._data.merge(data_months, on='date', how='left')
    model_GHI._data['B'] = model_GHI._data['B'].fillna(0)
    model_GHI._NM_model = NM_model_GHI
    
    return {
        'model_Ct': model_Ct,
        'model_GHI': model_GHI,
        'class': ['solarModelMixture', 'dict']
    }


# Helper functions (would need to be implemented separately)
def pmixnorm(x, means, sd, p):
    """CDF of Gaussian mixture."""
    pass

def qmixnorm(q, means, sd, p):
    """Quantile function of Gaussian mixture."""
    pass

def dmixnorm(x, means, sd, p):
    """PDF of Gaussian mixture."""
    pass

def GM_moments(means, sigma, probs):
    """Compute moments of Gaussian mixture."""
    pass

def mvgaussian_mixture(x, components):
    """Fit multivariate Gaussian mixture."""
    pass
