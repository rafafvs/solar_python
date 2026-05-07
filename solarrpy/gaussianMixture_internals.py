import numpy as np
import pandas as pd
from scipy.stats import norm
from typing import Dict
from sklearn.mixture import GaussianMixture
from scipy.optimize import minimize
from typing import Dict, Optional

def gm_moments(means: np.ndarray, sd: np.ndarray, alpha: np.ndarray) -> pd.DataFrame:
    """
    Compute the first four non-central and central moments, mean, variance, 
    skewness, and excess kurtosis for a Gaussian Mixture.
    """
    # Ensure inputs are numpy arrays for vectorized operations
    means = np.asarray(means, dtype=float)
    sd = np.asarray(sd, dtype=float)
    alpha = np.asarray(alpha, dtype=float)
    
    # First moment (Expectation)
    m1 = np.sum(means * alpha)
    
    # Non-central moments
    m2 = np.sum((means**2 + sd**2) * alpha)
    m3 = np.sum((3 * means * sd**2 + means**3) * alpha)
    m4 = np.sum((3 * sd**4 + 6 * means**2 * sd**2 + means**4) * alpha)
    
    # Central moments: shift means by the global expectation
    delta_k = means - m1
    mu2 = np.sum((delta_k**2 + sd**2) * alpha)
    mu3 = np.sum((3 * delta_k * sd**2 + delta_k**3) * alpha)
    mu4 = np.sum((3 * sd**4 + 6 * delta_k**2 * sd**2 + delta_k**4) * alpha)
    
    # Derived statistics
    v_x = mu2
    sk_x = mu3 / (v_x**(1.5))
    kt_x = (mu4 / (v_x**2)) - 3.0

    # Return as a 1-row DataFrame to mimic dplyr::tibble behavior
    return pd.DataFrame({
        'm1': [m1],
        'm2': [m2],
        'm3': [m3],
        'm4': [m4],
        'mu2': [mu2],
        'mu3': [mu3],
        'mu4': [mu4],
        'mean': [m1],
        'variance': [v_x],
        'skewness': [sk_x],
        'kurtosis': [kt_x]
    })

def gm_moments_match(d: float, m1: float = 0.0, m2: float = 1.0, m3: float = 0.0, p: float = 0.5) -> Dict[str, np.ndarray]:
    """
    Match the first three moments of a 2-component Gaussian Mixture.
    """
    # Mixing probabilities
    probs = np.array([p, 1 - p], dtype=float)
    
    # 1. Solve for means using the first moment and distance 'd'
    mu1 = m1 + (1 - p) * d
    mu2 = m1 - p * d
    means = np.array([mu1, mu2], dtype=float)
    
    # 2. Solve for variances using the second and third moments
    # delta represents the difference in variances: var1 - var2
    delta = (m3 - p * (1 - p) * ((1 - p)**2 - p**2) * d**3) / (3 * p * (1 - p) * d)
    
    var2 = m2 - p * (1 - p) * d**2 - p * delta
    var1 = var2 + delta
    
    # 3. Convert variances to standard deviations
    sigma = np.sqrt(np.array([var1, var2], dtype=float))
    
    return {
        'means': means,
        'sigma': sigma,
        'probs': probs
    }

def gm_loglik(means: np.ndarray, sd: np.ndarray, alpha: np.ndarray, x: np.ndarray = None) -> float:
    """
    Compute the log-likelihood of a Gaussian Mixture.
    """
    if x is not None:
        # Ensure inputs are numpy arrays
        x = np.asarray(x, dtype=float)
        means = np.asarray(means, dtype=float)
        sd = np.asarray(sd, dtype=float)
        alpha = np.asarray(alpha, dtype=float)
        
        # Calculate the density for each component for each observation.
        # x[:, np.newaxis] reshapes x to an (N, 1) column vector.
        # means, sd, and alpha are (K,) row vectors.
        # Broadcasting creates an (N, K) matrix of component likelihoods.
        component_densities = norm.pdf(x[:, np.newaxis], loc=means, scale=sd)
        
        # Multiply by mixture weights and sum across components (axis=1) 
        # to get the total mixture density per observation (N,)
        mix_density = np.sum(alpha * component_densities, axis=1)
        
        # Sum the log of the densities to get the total log-likelihood
        loss = np.sum(np.log(mix_density))
    else:
        loss = 0.0
        
    return loss

def gm_fit(x: np.ndarray, components: int = 2, maxit: int = 30000, maxrestarts: int = 10) -> Dict[str, np.ndarray]:
    """
    Fit a Gaussian Mixture model using scikit-learn's Expectation-Maximization.
    """
    # scikit-learn expects 2D arrays of shape (n_samples, n_features)
    # We force the 1D input into a 2D column vector.
    x_matrix = np.asarray(x, dtype=float).reshape(-1, 1)
    
    # Initialize and fit the standard EM algorithm.
    # 'full' covariance type allows each component to have its own variance,
    # mapping exactly to mclust's modelNames = "V".
    gmm = GaussianMixture(
        n_components=components,
        covariance_type='full',
        max_iter=maxit,
        n_init=maxrestarts, # Number of initializations to run
        init_params='kmeans' # Standard robust initialization
    )
    
    gmm.fit(x_matrix)

    if not gmm.converged_:
        raise ValueError("Gaussian Mixture did not converge!")
    
    # Extract parameters. 
    # gmm.means_ is shape (n_components, 1) -> flatten to 1D
    # gmm.covariances_ is shape (n_components, 1, 1) for 1D 'full' -> flatten and sqrt
    # gmm.weights_ is shape (n_components,)
    means = gmm.means_.flatten()
    sd = np.sqrt(gmm.covariances_.flatten())
    p = gmm.weights_.flatten()
    
    # Enforce identifiability: Sort components by mean (Label Switching Fix)
    idx = np.argsort(means)
    
    return {
        'means': means[idx],
        'sd': sd[idx],
        'p': p[idx]
    }

def gm_fit_moments_match(x: np.ndarray, 
                         start_params: np.ndarray, 
                         mu_target: Optional[float] = None, 
                         var_target: Optional[float] = None, 
                         maxit: int = 1000, 
                         abstol: float = 1e-5) -> Dict[str, np.ndarray]:
    """
    Fit a 2-component Gaussian Mixture matching target mean and/or variance.
    """
    # Pre-process data: drop NaNs once to avoid calculating na.rm=True on every iteration
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    
    def objective(par):
        mu1, mu2, sd1, sd2, p1 = par
        p2 = 1.0 - p1
        
        # Calculate likelihoods
        density = p1 * norm.pdf(x, loc=mu1, scale=sd1) + p2 * norm.pdf(x, loc=mu2, scale=sd2)
        
        # Prevent log(0) underflow
        density = np.clip(density, a_min=1e-300, a_max=None)
        
        loglik = np.sum(np.log(density))
        return -loglik

    def constraints(par):
        mu1, mu2, sd1, sd2, p1 = par
        p2 = 1.0 - p1
        
        # Moments
        e_x = p1 * mu1 + p2 * mu2
        e_x2 = p1 * (sd1**2 + mu1**2) + p2 * (sd2**2 + mu2**2)
        v_x = e_x2 - e_x**2
        
        loss = []
        if mu_target is not None:
            loss.append(e_x - mu_target)
        if var_target is not None:
            loss.append(v_x - var_target)
            
        return np.array(loss)

    # Hardcoded bounds from the R implementation
    bounds = [
        (-3.0, 3.0),   # mu1
        (-3.0, 3.0),   # mu2
        (0.01, 3.0),   # sd1
        (0.01, 3.0),   # sd2
        (0.01, 0.99)   # p1
    ]
    
    # Define SciPy equality constraints
    cons = []
    if mu_target is not None or var_target is not None:
        cons.append({'type': 'eq', 'fun': constraints})
        
    # Run SLSQP Optimizer
    res = minimize(
        fun=objective,
        x0=start_params,
        method='SLSQP',
        bounds=bounds,
        constraints=cons,
        options={'maxiter': maxit, 'ftol': abstol}
    )
    
    if not res.success:
        # Note: Depending on pipeline architecture, you might want to raise an Exception here
        print(f"Optimization Warning: {res.message}")

    return {
        'means': res.x[0:2],
        'sd': res.x[2:4],
        'p': np.array([res.x[4], 1.0 - res.x[4]])
    }