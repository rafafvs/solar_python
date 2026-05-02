import numpy as np
from scipy.stats import norm
import warnings

class BoundTransform:
    """Bounded transformation functions mapping risk drivers to unbounded space."""
    
    def __init__(self, alpha=0.0, beta=1.0, link="invgumbel"):
        if alpha < 0:
            raise ValueError("Alpha cannot be lower than zero.")
        if alpha + beta > 1:
            raise ValueError("`alpha + beta` cannot be greater than one.")
            
        self._alpha = alpha
        self._beta = beta
        self.epsilon = 0.0
        self.set_transform(link)

    def set_transform(self, link):
        """Set the transform function g, g_inv, g_prime."""
        valid_links = ["invgumbel", "gumbel", "logis", "norm"]
        if link not in valid_links:
            raise ValueError(f"Link must be one of {valid_links}")
            
        self._link = link
        
        # We use numpy for vectorized operations over pandas Series or numpy arrays
        if link == "invgumbel":
            self._g = lambda x: np.log(-np.log(x))
            self._ig = lambda x: np.exp(-np.exp(x))
            self._g_prime = lambda x: 1 / (x * np.log(x))
            self._monotonicity = "decreasing"
        elif link == "gumbel":
            self._g = lambda x: -np.log(-np.log(x))
            self._ig = lambda x: np.exp(-np.exp(-x))
            self._g_prime = lambda x: -1 / (x * np.log(x))
            self._monotonicity = "increasing"
        elif link == "logis":
            self._g = lambda x: np.log(x / (1 - x))
            self._ig = lambda x: 1 / (1 + np.exp(-x))
            self._g_prime = lambda x: 1 / (x * (1 - x))
            self._monotonicity = "increasing"
        elif link == "norm":
            self._g = lambda x: norm.ppf(x)
            self._ig = lambda x: norm.cdf(x)
            self._g_prime = lambda x: 1 / norm.pdf(norm.ppf(x))
            self._monotonicity = "increasing"

    def X_prime(self, Xt):
        """Map the risk-driver Xt into the normalized risk-driver Xt_prime."""
        # Convert to numpy array for fast conditional bounding
        Xt = np.atleast_1d(Xt).astype(float)
        
        # Apply epsilon bounds to prevent log(0) or division issues
        upper_bound = self._alpha + self._beta
        Xt = np.where(Xt >= upper_bound, upper_bound - self.epsilon, Xt)
        Xt = np.where(Xt <= self._alpha, self._alpha + self.epsilon, Xt)
        
        res = (Xt - self._alpha) / self._beta
        return res[0] if res.size == 1 else res

    def iX_prime(self, Xt_prime):
        """Map the normalized variable Xt_prime to the risk driver Xt."""
        return self._alpha + self._beta * Xt_prime

    def Y(self, Xt_prime):
        """Map the normalized risk driver Xt_prime into the transformed variable Yt."""
        return self.g(Xt_prime)

    def iY(self, Yt):
        """Map the transformed variable Yt into the normalized risk driver Xt_prime."""
        return self.ig(Yt)

    def g(self, X_prime):
        return self._g(X_prime)

    def ig(self, Yt):
        return self._ig(Yt)

    def g_prime(self, X_prime):
        return self._g_prime(X_prime)

    def fit(self, Xt, epsilon=0.01, min_pos=1, max_pos=1):
        """Fit the best parameters alpha and beta from a given time series."""
        x = np.asarray(Xt)
        
        sorted_x = np.sort(x)
        
        xmin_val = sorted_x[min_pos - 1]
        xmax_val = sorted_x[-max_pos]
        
        self.epsilon = xmin_val * epsilon
        self._alpha = xmin_val - self.epsilon
        self._beta = (xmax_val - xmin_val) + 2 * self.epsilon
        
        return {
            "alpha": self._alpha,
            "beta": self._beta,
            "epsilon": self.epsilon,
            "Xt_min": xmin_val,
            "Xt_max": xmax_val
        }

    def bounds(self, target="Xt"):
        """Compute the bounds for the transformed variables."""
        valid_targets = ["Xt", "Yt", "Kt"]
        if target not in valid_targets:
            raise ValueError(f"target must be one of {valid_targets}")
            
        if target == "Yt":
            return {"Yt_min": -np.inf, "Yt_max": np.inf}
        elif target == "Kt":
            return {
                "Kt_min": 1 - self._alpha - self._beta, 
                "Kt_max": 1 - self._alpha
            }
        else: # target == "Xt"
            return {
                "Xt_min": self._alpha, 
                "Xt_max": self._alpha + self._beta
            }

    def update(self, alpha=None, beta=None):
        if alpha is None:
            alpha = self._alpha
        if beta is None:
            beta = self._beta
            
        if alpha < 0:
            warnings.warn("Error: alpha is lower than 0")
            return
        if alpha + beta > 1:
            warnings.warn("Error: alpha + beta is greater than 1")
            return
            
        self._alpha = alpha
        self._beta = beta

    @property
    def alpha(self):
        return self._alpha

    @property
    def beta(self):
        return self._beta

    @property
    def link(self):
        return self._link

    @property
    def monotonicity(self):
        return self._monotonicity