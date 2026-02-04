import numpy as np
from typing import Optional, Callable
from scipy import stats


class BoundTransform:
    """
    Base class for bounded transformations.
    
    This is the parent class that solarTransform inherits from.
    """
    
    def __init__(self, alpha: float = 0, beta: float = 1, link: str = "invgumbel"):
        """
        Initialize a bounded transform.
        
        Parameters:
        -----------
        alpha : float
            Lower bound parameter
        beta : float
            Width parameter
        link : str
            Link function type
        """
        self._alpha = alpha
        self._beta = beta
        self._link = link
        self._link_function = None
        self._inv_link_function = None
    
    def set_transform(self, link: str):
        """
        Set the transformation link function.
        
        Parameters:
        -----------
        link : str
            Link function: 'invgumbel', 'gumbel', 'logit', 'norm'
        """
        self._link = link
        
        if link == "invgumbel":
            # Inverse Gumbel (Type II extreme value)
            self._link_function = lambda x: -np.log(-np.log(1 - x))
            self._inv_link_function = lambda y: 1 - np.exp(-np.exp(-y))
        elif link == "gumbel":
            # Gumbel (Type I extreme value)
            self._link_function = lambda x: -np.log(-np.log(x))
            self._inv_link_function = lambda y: np.exp(-np.exp(-y))
        elif link == "logit":
            # Logistic
            self._link_function = lambda x: np.log(x / (1 - x))
            self._inv_link_function = lambda y: 1 / (1 + np.exp(-y))
        elif link == "norm":
            # Normal (probit)
            self._link_function = lambda x: stats.norm.ppf(x)
            self._inv_link_function = lambda y: stats.norm.cdf(y)
        else:
            raise ValueError(f"Unknown link function: {link}")
    
    def X_prime(self, Xt):
        """
        Normalize X_t to (0,1) range.
        
        Parameters:
        -----------
        Xt : array-like
            Risk driver in (alpha, alpha + beta)
        
        Returns:
        --------
        array-like
            Normalized X'_t in (0, 1)
        """
        return (Xt - self._alpha) / self._beta
    
    def iX_prime(self, Xt_prime):
        """
        Inverse normalization from (0,1) to (alpha, alpha + beta).
        
        Parameters:
        -----------
        Xt_prime : array-like
            Normalized risk driver in (0, 1)
        
        Returns:
        --------
        array-like
            X_t in (alpha, alpha + beta)
        """
        return self._alpha + self._beta * Xt_prime
    
    def Y(self, Xt_prime):
        """
        Apply link function: g(X'_t).
        
        Parameters:
        -----------
        Xt_prime : array-like
            Normalized risk driver in (0, 1)
        
        Returns:
        --------
        array-like
            Transformed variable Y_t in (-∞, ∞)
        """
        return self._link_function(Xt_prime)
    
    def iY(self, Yt):
        """
        Apply inverse link function: g^(-1)(Y_t).
        
        Parameters:
        -----------
        Yt : array-like
            Transformed variable in (-∞, ∞)
        
        Returns:
        --------
        array-like
            Normalized X'_t in (0, 1)
        """
        return self._inv_link_function(Yt)
    
    def update(self, alpha: float, beta: float):
        """
        Update transformation parameters.
        
        Parameters:
        -----------
        alpha : float
            New alpha parameter
        beta : float
            New beta parameter
        """
        if alpha < 0:
            raise ValueError("Alpha cannot be lower than zero.")
        if alpha + beta > 1:
            raise ValueError("`alpha + beta` cannot be greater than one.")
        
        self._alpha = alpha
        self._beta = beta
    
    def fit(self, Xt, threshold: float = 0.01, min_pos: int = 1, max_pos: int = 1):
        """
        Fit transformation parameters from data.
        
        Parameters:
        -----------
        Xt : array-like
            Risk driver values
        threshold : float
            Threshold for outlier detection
        min_pos : int
            Position of minimum (1 = smallest value)
        max_pos : int
            Position of maximum (1 = largest value)
        
        Returns:
        --------
        dict
            Dictionary with fitted parameters
        """
        Xt = np.array(Xt)
        
        # Sort to find min/max at specified positions
        sorted_Xt = np.sort(Xt)
        
        # Get min and max based on position
        Xt_min = sorted_Xt[min_pos - 1]
        Xt_max = sorted_Xt[-(max_pos)]
        
        # Estimate alpha and beta
        alpha = Xt_min - threshold
        beta = Xt_max - Xt_min + 2 * threshold
        
        # Ensure constraints
        if alpha < 0:
            alpha = 0
        if alpha + beta > 1:
            beta = 1 - alpha
        
        return {
            'alpha': alpha,
            'beta': beta,
            'Xt_min': Xt_min,
            'Xt_max': Xt_max
        }
    
    @property
    def alpha(self):
        """Get alpha parameter."""
        return self._alpha
    
    @property
    def beta(self):
        """Get beta parameter."""
        return self._beta
    
    @property
    def link(self):
        """Get link function name."""
        return self._link


class SolarTransform(BoundTransform):
    """
    Solar Model transformation functions.
    
    This class provides transformations between solar radiation (R_t),
    risk drivers (X_t), and transformed variables (Y_t).
    
    Example:
    --------
    st = SolarTransform(alpha=0.05, beta=0.90, link='invgumbel')
    Xt = st.X(Rt=5.0, Ct=6.0)
    Yt = st.RY(Rt=5.0, Ct=6.0)
    Rt_back = st.iRY(Yt, Ct=6.0)
    
    Version: 1.0.1
    """
    
    VERSION = "1.0.1"
    
    def __init__(self, alpha: float = 0, beta: float = 1, link: str = "invgumbel"):
        """
        Initialize a SolarTransform object.
        
        Parameters:
        -----------
        alpha : float
            α transformation parameter (lower bound shift)
        beta : float
            β transformation parameter (range width)
        link : str
            Link function: 'invgumbel', 'gumbel', 'logit', 'norm'
        
        Raises:
        -------
        ValueError
            If alpha < 0 or alpha + beta > 1
        """
        # Control parameters
        if alpha < 0:
            raise ValueError("Alpha cannot be lower than zero.")
        if alpha + beta > 1:
            raise ValueError("`alpha + beta` cannot be greater than one.")
        
        # Initialize parent class
        super().__init__(alpha, beta, link)
        
        # Additional parameter
        self.epsilon = 0
    
    def X(self, Rt, Ct):
        """
        Map solar radiation R_t to risk driver X_t.
        
        The function computes:
            X(R_t) = 1 - R_t / C_t
        
        Parameters:
        -----------
        Rt : float or array-like
            Solar radiation in [C_t(1-α-β), C_t(1-α)]
        Ct : float or array-like
            Clear sky radiation
        
        Returns:
        --------
        float or array-like
            Risk driver X_t in (α, α+β)
        """
        return 1 - Rt / Ct
    
    def iX(self, Xt, Ct):
        """
        Map risk driver X_t to solar radiation R_t.
        
        The function computes:
            X^(-1)(X_t) = C_t(1 - X_t)
        
        Parameters:
        -----------
        Xt : float or array-like
            Risk driver in (α, α+β)
        Ct : float or array-like
            Clear sky radiation
        
        Returns:
        --------
        float or array-like
            Solar radiation R_t in [C_t(1-α-β), C_t(1-α)]
        """
        return Ct * (1 - Xt)
    
    def eta(self, Rt, Ct):
        """
        Map solar radiation R_t to normalized variable X'_t.
        
        The function computes:
            η(R_t) = (1/β)(1 - α - R_t/C_t)
        
        Parameters:
        -----------
        Rt : float or array-like
            Solar radiation in [C_t(1-α-β), C_t(1-α)]
        Ct : float or array-like
            Clear sky radiation
        
        Returns:
        --------
        float or array-like
            Normalized risk driver X'_t in (0, 1)
        """
        return self.X_prime(self.X(Rt, Ct))
    
    def ieta(self, Xt_prime, Ct):
        """
        Map normalized variable X'_t to solar radiation R_t.
        
        The function computes:
            η^(-1)(X'_t) = C_t(1 - α - β · X'_t)
        
        Parameters:
        -----------
        Xt_prime : float or array-like
            Normalized risk driver in (0, 1)
        Ct : float or array-like
            Clear sky radiation
        
        Returns:
        --------
        float or array-like
            Solar radiation R_t in [C_t(1-α-β), C_t(1-α)]
        """
        return self.iX(self.iX_prime(Xt_prime), Ct)
    
    def RY(self, Rt, Ct):
        """
        Convert solar radiation R_t to transformed variable Y_t.
        
        The function computes:
            RY(R_t) = g((1/β)(1 - α - R_t/C_t))
        
        where g is the link function.
        
        Parameters:
        -----------
        Rt : float or array-like
            Solar radiation in [C_t(1-α-β), C_t(1-α)]
        Ct : float or array-like
            Clear sky radiation
        
        Returns:
        --------
        float or array-like
            Transformed variable Y_t in (-∞, ∞)
        """
        # Cloudiness index
        Xt = self.X(Rt, Ct)
        
        # Normalized risk driver
        Xt_prime = self.X_prime(Xt)
        
        # Transformed variable
        return self.Y(Xt_prime)
    
    def iRY(self, Yt, Ct):
        """
        Convert transformed variable Y_t to solar radiation R_t.
        
        The function computes:
            iRY(Y_t) = C_t(1 - α - β · g^(-1)(Y_t))
        
        where g^(-1) is the inverse link function.
        
        Parameters:
        -----------
        Yt : float or array-like
            Transformed variable in (-∞, ∞)
        Ct : float or array-like
            Clear sky radiation
        
        Returns:
        --------
        float or array-like
            Solar radiation R_t in [C_t(1-α-β), C_t(1-α)]
        """
        # Normalized risk driver
        Xt_prime = self.iY(Yt)
        
        # Cloudiness index
        Xt = self.iX_prime(Xt_prime)
        
        # Solar radiation
        return self.iX(Xt, Ct)
    
    def __repr__(self):
        """String representation."""
        return (f"SolarTransform(alpha={self.alpha:.6f}, beta={self.beta:.6f}, "
                f"link='{self.link}', version={self.VERSION})")


# Example usage
if __name__ == "__main__":
    # Create a solar transform
    st = SolarTransform(alpha=0.05, beta=0.90, link='invgumbel')
    
    # Example values
    Rt = 5.0  # Solar radiation
    Ct = 6.0  # Clear sky radiation
    
    # Forward transformations
    Xt = st.X(Rt, Ct)
    print(f"Risk driver X_t: {Xt:.6f}")
    
    Xt_prime = st.eta(Rt, Ct)
    print(f"Normalized X'_t: {Xt_prime:.6f}")
    
    Yt = st.RY(Rt, Ct)
    print(f"Transformed Y_t: {Yt:.6f}")
    
    # Backward transformation
    Rt_back = st.iRY(Yt, Ct)
    print(f"Recovered R_t: {Rt_back:.6f}")
    print(f"Original R_t: {Rt:.6f}")
    print(f"Difference: {abs(Rt - Rt_back):.10f}")
