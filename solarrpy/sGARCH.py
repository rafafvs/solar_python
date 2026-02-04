import numpy as np
import pandas as pd
from typing import Optional, Tuple, Dict, List
import warnings
from arch import arch_model


class sGARCH:
    """
    Implementation of standard GARCH(p,q) model as a Python class.
    
    This class provides methods to fit, filter, and forecast GARCH models,
    similar to the rugarch package in R.
    
    Example:
    --------
    # Initialize GARCH(1,1) model
    model = sGARCH(arch_order=1, garch_order=1)
    
    # Fit the model
    model.fit(x=returns)
    
    # Filter to get conditional variance
    sigma2 = model.filter(x=returns)
    
    # Forecast next step
    next_var = model.next_step(eps0=[0.01], sigma20=[0.005])
    
    Version: 1.0.1
    Keywords: GARCH
    """
    
    VERSION = "1.0.1"
    
    def __init__(self, arch_order: int = 1, garch_order: int = 1, 
                 mode: str = "unitOmega"):
        """
        Initialize a standard GARCH model.
        
        Parameters:
        -----------
        arch_order : int
            ARCH order (p)
        garch_order : int
            GARCH order (q)
        mode : str
            Unconditional variance constraint: 'unitOmega', 'targetSigma2', or 'freeOmega'
        """
        # ARCH order
        self._arch_order = arch_order
        
        # GARCH order
        self._garch_order = garch_order
        
        # Unconditional variance constraint
        if mode not in ['unitOmega', 'targetSigma2', 'freeOmega']:
            raise ValueError("mode must be one of: 'unitOmega', 'targetSigma2', 'freeOmega'")
        self._mode = mode
        
        # Initialize parameters
        self._omega = np.array([1.0])
        self._alpha = np.zeros(arch_order) if arch_order > 0 else np.array([0.0])
        self._beta = np.zeros(garch_order) if garch_order > 0 else np.array([0.0])
        
        # Matrix components
        self._b = self._GARCH_vector_b(arch_order, garch_order)
        self._A = np.array([[]])
        self._d = np.array([[]])
        
        # Data storage
        self._x = None
        self._w = None
        self._log_likelihoods = np.array([])
        self._std_errors = np.array([np.nan, np.nan, np.nan])
        self._sigma20 = None
    
    @staticmethod
    def _GARCH_vector_b(arch_order: int, garch_order: int) -> np.ndarray:
        """
        Create the b vector for GARCH companion form.
        
        Parameters:
        -----------
        arch_order : int
            ARCH order
        garch_order : int
            GARCH order
        
        Returns:
        --------
        np.ndarray
            Vector b
        """
        max_order = max(arch_order, garch_order)
        b = np.zeros((max_order, 1))
        b[0, 0] = 1.0
        return b
    
    @staticmethod
    def _ARMA_companion_matrix(phi: np.ndarray, theta: np.ndarray) -> np.ndarray:
        """
        Create ARMA companion matrix.
        
        Parameters:
        -----------
        phi : np.ndarray
            AR coefficients
        theta : np.ndarray
            MA coefficients
        
        Returns:
        --------
        np.ndarray
            Companion matrix
        """
        p = len(phi)
        q = len(theta)
        m = max(p, q)
        
        if m == 0:
            return np.array([[]])
        
        # Create companion matrix
        A = np.zeros((m, m))
        
        # First row contains AR and MA coefficients
        if p > 0:
            A[0, :p] = phi
        if q > 0:
            A[0, :q] += theta
        
        # Identity submatrix
        if m > 1:
            A[1:, :-1] = np.eye(m - 1)
        
        return A
    
    def fit(self, x: np.ndarray, weights: Optional[np.ndarray] = None):
        """
        Fit the GARCH model.
        
        Parameters:
        -----------
        x : np.ndarray
            Time series to be fitted
        weights : np.ndarray, optional
            Custom weights for observations
        """
        # Store time series
        self._x = np.array(x)
        
        # Flexible weights
        if weights is None:
            self._w = np.ones(len(x))
        else:
            self._w = np.where(weights == 0, 0, 1)
        
        # Fit using arch package
        try:
            # Create GARCH specification
            am = arch_model(
                x,
                vol='GARCH',
                p=self._arch_order,
                q=self._garch_order,
                mean='Zero',
                rescale=False
            )
            
            # Fit the model
            res = am.fit(disp='off')
            
            # Extract parameters
            params = res.params
            
            # Update coefficients
            self._omega = np.array([params['omega']])
            
            if self._arch_order > 0:
                self._alpha = np.array([params[f'alpha[{i+1}]'] 
                                       for i in range(self._arch_order)])
            
            if self._garch_order > 0:
                self._beta = np.array([params[f'beta[{i+1}]'] 
                                      for i in range(self._garch_order)])
            
            # Log-likelihood
            self._log_likelihoods = -res.resid**2 / (2 * res.conditional_volatility**2) - \
                                   0.5 * np.log(2 * np.pi * res.conditional_volatility**2)
            
            # Standard errors
            std_errs = []
            std_errs.append(res.std_err['omega'])
            for i in range(self._arch_order):
                std_errs.append(res.std_err[f'alpha[{i+1}]'])
            for i in range(self._garch_order):
                std_errs.append(res.std_err[f'beta[{i+1}]'])
            self._std_errors = np.array(std_errs)
            
            # Unconditional variance
            self._sigma20 = res.conditional_volatility.mean()**2
            
        except Exception as e:
            warnings.warn(f"GARCH fitting failed: {e}. Using fallback method.")
            self._fallback_fit(x, weights)
        
        # Update matrix components
        alpha_nonzero = self._alpha[self._alpha != 0]
        beta_nonzero = self._beta[self._beta != 0]
        
        if len(alpha_nonzero) > 0 or len(beta_nonzero) > 0:
            self._A = self._ARMA_companion_matrix(alpha_nonzero, beta_nonzero)
        
        # Intercept vector
        total_order = self._arch_order + self._garch_order
        if total_order > 0:
            self._d = np.zeros((total_order, 1))
            self._d[0, 0] = self._omega[0]
    
    def _fallback_fit(self, x: np.ndarray, weights: Optional[np.ndarray]):
        """Fallback fitting method using simple moment estimation."""
        # Simple moment-based estimation
        self._omega = np.array([np.var(x) * 0.1])
        
        if self._arch_order > 0:
            self._alpha = np.full(self._arch_order, 0.1)
        
        if self._garch_order > 0:
            self._beta = np.full(self._garch_order, 0.8)
        
        self._std_errors = np.full(1 + self._arch_order + self._garch_order, np.nan)
        self._sigma20 = np.var(x)
    
    def filter(self, x: np.ndarray, eps0: Optional[np.ndarray] = None,
               sigma20: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Filter method to compute GARCH variance.
        
        Parameters:
        -----------
        x : np.ndarray
            Time series to be filtered
        eps0 : np.ndarray, optional
            Initial epsilons (length p+q)
        sigma20 : np.ndarray, optional
            Initial variances (length p+q)
        
        Returns:
        --------
        np.ndarray
            Conditional variances
        """
        return self._sGARCH_filter(x, self._omega[0], self._alpha, self._beta, 
                                   eps0, sigma20)
    
    def _sGARCH_filter(self, x: np.ndarray, omega: float, 
                       alpha: np.ndarray, beta: np.ndarray,
                       eps0: Optional[np.ndarray] = None,
                       sigma20: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Compute GARCH filtered variances.
        
        Parameters:
        -----------
        x : np.ndarray
            Time series
        omega : float
            Intercept
        alpha : np.ndarray
            ARCH coefficients
        beta : np.ndarray
            GARCH coefficients
        eps0 : np.ndarray, optional
            Initial residuals
        sigma20 : np.ndarray, optional
            Initial variances
        
        Returns:
        --------
        np.ndarray
            Conditional variances
        """
        n = len(x)
        p = len(alpha)
        q = len(beta)
        
        # Initialize
        if eps0 is None:
            eps0 = np.zeros(p)
        if sigma20 is None:
            sigma20 = np.full(q, omega / (1 - np.sum(alpha) - np.sum(beta)))
        
        # Storage
        sigma2 = np.zeros(n)
        eps = np.concatenate([eps0, x])
        
        # Recursion
        for t in range(n):
            # ARCH component
            arch_sum = np.sum(alpha * eps[t:t+p][::-1]**2) if p > 0 else 0
            
            # GARCH component
            if t < q:
                available_sigma2 = np.concatenate([sigma20[:q-t][::-1], sigma2[:t][::-1]])
            else:
                available_sigma2 = sigma2[t-q:t][::-1]
            garch_sum = np.sum(beta * available_sigma2) if q > 0 else 0
            
            # Update variance
            sigma2[t] = omega + arch_sum + garch_sum
        
        return sigma2
    
    def update(self, coefficients: Optional[Dict[str, float]] = None):
        """
        Update the coefficients of the model.
        
        Parameters:
        -----------
        coefficients : dict, optional
            Dictionary of coefficient names and values
        """
        if coefficients is None:
            return
        
        # Get current coefficients
        current_coefs = self.coefficients
        
        # Update omega
        if 'omega' in coefficients:
            self._omega = np.array([coefficients['omega']])
            if 'omega' in range(len(self._std_errors)):
                self._std_errors[0] = np.nan
        
        # Update alpha parameters
        for i in range(self._arch_order):
            key = f'alpha{i+1}'
            if key in coefficients:
                self._alpha[i] = coefficients[key]
                if i + 1 < len(self._std_errors):
                    self._std_errors[i + 1] = np.nan
        
        # Update beta parameters
        for i in range(self._garch_order):
            key = f'beta{i+1}'
            if key in coefficients:
                self._beta[i] = coefficients[key]
                if i + 1 + self._arch_order < len(self._std_errors):
                    self._std_errors[i + 1 + self._arch_order] = np.nan
        
        # Update matrix components
        alpha_nonzero = self._alpha[self._alpha != 0]
        beta_nonzero = self._beta[self._beta != 0]
        
        if len(alpha_nonzero) > 0 or len(beta_nonzero) > 0:
            self._A = self._ARMA_companion_matrix(alpha_nonzero, beta_nonzero)
        
        # Update intercept vector
        total_order = sum(self.order)
        if total_order > 0:
            self._d = np.zeros((total_order, 1))
            self._d[0, 0] = self._omega[0]
    
    def update_std_errors(self, std_errors: Optional[Dict[str, float]] = None):
        """
        Update the standard errors of the parameters.
        
        Parameters:
        -----------
        std_errors : dict, optional
            Dictionary of standard errors
        """
        if std_errors is None or len(std_errors) == 0:
            return
        
        # Create new std errors array
        new_std_errors = self._std_errors.copy()
        
        # Map of parameter names to indices
        idx = 0
        param_map = {'omega': idx}
        idx += 1
        
        for i in range(self._arch_order):
            param_map[f'alpha{i+1}'] = idx
            idx += 1
        
        for i in range(self._garch_order):
            param_map[f'beta{i+1}'] = idx
            idx += 1
        
        # Update values
        for key, value in std_errors.items():
            if key in param_map:
                new_std_errors[param_map[key]] = value
        
        self._std_errors = new_std_errors
    
    def next_step(self, eps0: np.ndarray, sigma20: np.ndarray) -> float:
        """
        Compute next step GARCH variance forecast.
        
        Parameters:
        -----------
        eps0 : np.ndarray
            Recent residuals (length p)
        sigma20 : np.ndarray
            Recent variances (length q)
        
        Returns:
        --------
        float
            Next period variance forecast
        """
        # ARCH component
        arch_sum = np.sum(self._alpha * np.array(eps0)**2) if self._arch_order > 0 else 0
        
        # GARCH component
        garch_sum = np.sum(self._beta * np.array(sigma20)) if self._garch_order > 0 else 0
        
        # Next variance
        return self._omega[0] + arch_sum + garch_sum
    
    def __repr__(self) -> str:
        """Print method for sGARCH class."""
        # Format parameters with std. errors
        def format_param(val, stderr):
            return f"{val:.4f} ({stderr:.3f})" if not np.isnan(stderr) else f"{val:.4f} (NA)"
        
        # Intercept
        omega_str = format_param(self._omega[0], self._std_errors[0])
        
        # ARCH parameters
        arch_strs = []
        for i in range(self._arch_order):
            idx = 1 + i
            arch_strs.append(format_param(self._alpha[i], 
                                         self._std_errors[idx] if idx < len(self._std_errors) else np.nan))
        alpha_str = " ".join(arch_strs) if arch_strs else "None"
        
        # GARCH parameters
        garch_strs = []
        for i in range(self._garch_order):
            idx = 1 + self._arch_order + i
            garch_strs.append(format_param(self._beta[i],
                                          self._std_errors[idx] if idx < len(self._std_errors) else np.nan))
        beta_str = " ".join(garch_strs) if garch_strs else "None"
        
        # Model name
        model_name = f"GARCH({self._arch_order}, {self._garch_order})"
        
        lines = [
            f"--------------------- {model_name} ---------------------",
            f"ARCH: {self._arch_order > 0}",
            f"GARCH: {self._garch_order > 0}",
            f"Version: {self.VERSION}",
            f"Unconditional Variance: {self._sigma20}",
            "------------------------------------------------",
            f"Intercept: {omega_str}",
            f"Alpha: {alpha_str}",
            f"Beta: {beta_str}",
            f"Log-lik: {self.loglik:.3f}" if self.loglik is not None else "Log-lik: NA"
        ]
        
        return "\n".join(lines)
    
    # ==================== Properties ====================
    
    @property
    def arch_order(self) -> int:
        """ARCH order."""
        return self._arch_order
    
    @property
    def garch_order(self) -> int:
        """GARCH order."""
        return self._garch_order
    
    @property
    def order(self) -> Tuple[int, int]:
        """Model orders (ARCH, GARCH)."""
        return (self._arch_order, self._garch_order)
    
    @property
    def omega(self) -> np.ndarray:
        """Intercept parameter."""
        return self._omega
    
    @property
    def alpha(self) -> np.ndarray:
        """ARCH parameters."""
        return self._alpha
    
    @property
    def beta(self) -> np.ndarray:
        """GARCH parameters."""
        return self._beta
    
    @property
    def coefficients(self) -> Dict[str, float]:
        """Model coefficients as dictionary."""
        coefs = {'omega': self._omega[0]}
        
        for i in range(self._arch_order):
            coefs[f'alpha{i+1}'] = self._alpha[i]
        
        for i in range(self._garch_order):
            coefs[f'beta{i+1}'] = self._beta[i]
        
        return coefs
    
    @property
    def A(self) -> np.ndarray:
        """Companion matrix."""
        return self._A
    
    @property
    def b(self) -> np.ndarray:
        """Vector b."""
        return self._b
    
    @property
    def d(self) -> np.ndarray:
        """Intercept vector."""
        return self._d
    
    @property
    def std_errors(self) -> np.ndarray:
        """Standard errors of coefficients."""
        return self._std_errors
    
    @property
    def sigma2_inf(self) -> float:
        """Long-term unconditional variance."""
        denom = 1 - np.sum(self._alpha) - np.sum(self._beta)
        if denom <= 0:
            warnings.warn("Model is not stationary (sum of alpha + beta >= 1)")
            return np.inf
        return self._omega[0] / denom
    
    @property
    def loglik(self) -> Optional[float]:
        """Model log-likelihood."""
        if len(self._log_likelihoods) == 0:
            return None
        return np.sum(self._log_likelihoods)
    
    @property
    def tidy(self) -> pd.DataFrame:
        """Tibble with estimated parameters and std. errors."""
        coefs = self.coefficients
        
        terms = list(coefs.keys())
        estimates = list(coefs.values())
        std_errs = list(self._std_errors[:len(terms)])
        
        return pd.DataFrame({
            'term': terms,
            'estimate': estimates,
            'std_error': std_errs
        })


# Example usage
if __name__ == "__main__":
    # Simulate some returns data
    np.random.seed(42)
    n = 1000
    returns = np.random.randn(n) * 0.01
    
    # Create and fit GARCH(1,1) model
    model = sGARCH(arch_order=1, garch_order=1)
    model.fit(returns)
    
    print(model)
    print("\nTidy output:")
    print(model.tidy)
    
    # Filter to get conditional variance
    sigma2 = model.filter(returns)
    print(f"\nMean conditional variance: {np.mean(sigma2):.6f}")
    
    # Next step forecast
    next_var = model.next_step(eps0=[returns[-1]], sigma20=[sigma2[-1]])
    print(f"Next period variance forecast: {next_var:.6f}")
