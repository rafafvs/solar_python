"""
ARMA helper functions — Python translation of the R source.

Functions included:
    ARMA_vector_b          — state-space selection vector
    ARMA_companion_matrix  — companion / transition matrix
    ARMA_expectation       — h-step conditional mean
    ARMA_variance          — h-step forecast variance
    ARMA_covariance        — h-step forecast covariance
    ARMA_next_step         — iterate the state forward h steps
    ARMA_filter            — one-step-ahead filter (replaces the C backend)
    ARMA_forecast          — full state trajectory
    AR_variance            — long-run AR variance (closed-form for p≤4)
    AR_params_to_zeta      — phi → unconstrained zeta  (step-down Durbin-Levinson + atanh)
    ARMA_params_to_zeta    — ARMA wrapper for AR_params_to_zeta
    AR_params_to_phi       — zeta → constrained phi   (tanh + step-up Durbin-Levinson + Jacobian)
    ARMA_params_to_phi     — ARMA wrapper for AR_params_to_phi
    ARMA_loss_logLik       — Gaussian log-likelihood (constrained parametrisation)
    ARMA_loss_CSS          — conditional sum-of-squares loss
    ARMA_fit               — full estimation pipeline

Version: 1.0.0
"""

from __future__ import annotations

import warnings
import numpy as np
import pandas as pd
import numdifftools as nd # Equivalent to R's numDeriv
from scipy.optimize import minimize
from scipy.linalg import solve   # solve(A, b) → A^{-1} b
from scipy.stats import norm

# ---------------------------------------------------------------------------
# ARMA_vector_b
# ---------------------------------------------------------------------------

def arma_vector_b(ar_order: int, ma_order: int) -> np.ndarray:
    """
    Construct the vector b of an ARMA model state-space representation.
    
    Args:
        ar_order: Order of the AR component (p).
        ma_order: Order of the MA component (q).
        
    Returns:
        A column vector (numpy array) of shape (p + q, 1).
    """
    # Initialize component lists
    e_p = []
    e_q = []

    # Construct AR part: [1, 0, ..., 0] of length p
    if ar_order > 0:
        e_p = np.zeros(ar_order)
        e_p[0] = 1

    # Construct MA part: [1, 0, ..., 0] of length q
    if ma_order > 0:
        e_q = np.zeros(ma_order)
        e_q[0] = 1

    # Concatenate and reshape to column vector (n x 1)
    # np.concatenate handles empty lists if one order is 0
    b = np.concatenate([e_p, e_q]).reshape(-1, 1)
    
    return b

# ---------------------------------------------------------------------------
# ARMA_companion_matrix
# ---------------------------------------------------------------------------

def arma_companion_matrix(phi=None, theta=None):
    """
    Construct the companion matrix of an ARMA model.
    
    Args:
        phi: Array-like, AR parameters.
        theta: Array-like, MA parameters.
        
    Returns:
        A square numpy array of dimension (p+q).
        The result includes 'arOrder' and 'maOrder' as custom attributes 
        stored in a dict (or simply return them alongside).
    """
    # Standardize inputs to numpy arrays
    phi = np.atleast_1d(phi) if phi is not None else np.array([])
    theta = np.atleast_1d(theta) if theta is not None else np.array([])
    
    p = len(phi)
    q = len(theta)
    dim = p + q
    
    if dim == 0:
        return np.array([[]])

    # 1. Construct the first row: [phi, theta]
    first_row = np.concatenate([phi, theta]).reshape(1, -1)
    
    blocks = [first_row]

    # 2. Construct L_p (AR shift block)
    if p > 1:
        # Identity (p-1) x (p-1)
        I_p = np.eye(p - 1)
        # Add column of zeros to make it (p-1) x p
        L_p_ar = np.hstack([I_p, np.zeros((p - 1, 1))])
        # Add columns of zeros for MA part to make it (p-1) x (p+q)
        L_p = np.hstack([L_p_ar, np.zeros((p - 1, q))])
        blocks.append(L_p)

    # 3. Add separator zero row if both p and q exist
    if p >= 1 and q >= 1:
        blocks.append(np.zeros((1, dim)))

    # 4. Construct L_q (MA shift block)
    if q > 1:
        # Identity (q-1) x (q-1)
        I_q = np.eye(q - 1)
        # Add column of zeros to make it (q-1) x q
        L_q_ma = np.hstack([I_q, np.zeros((q - 1, 1))])
        # Add columns of zeros for AR part (at the start) to make it (q-1) x (p+q)
        L_q = np.hstack([np.zeros((q - 1, p)), L_q_ma])
        blocks.append(L_q)

    # Combine all blocks
    A = np.vstack(blocks)
    
    # Store orders as metadata (Pythonic approach: use an object or dictionary)
    # For a direct translation, we attach it to the array object
    A.flags.writeable = True # Ensure we can add attributes
    A.ar_order = p
    A.ma_order = q
    
    return A

# ---------------------------------------------------------------------------
# ARMA_expectation
# ---------------------------------------------------------------------------

def arma_expectation(h: int = 10, x0: np.ndarray = None, a: np.ndarray = None, 
                     b: np.ndarray = None, intercept: float = 0.0):
    """
    Compute the conditional mean of an ARMA model.
    """
    # Dimension of the state space
    pq = a.shape[0]
    
    # Build intercept vector c_ (pq x 1)
    c_vec = np.zeros((pq, 1))
    c_vec[0, 0] = intercept
    
    # Ensure x0 is a column vector (pq x 1)
    x0 = np.asarray(x0).reshape(-1, 1)
    
    # Identity matrix
    identity = np.eye(pq)
    
    # Pre-compute (I - A)^-1 * c_ vector
    # Equivalent to R's solve(I - A) %*% c_
    # Note: This assumes (I - A) is non-singular (stationary process)
    i_a_inv_c = np.linalg.solve(identity - a, c_vec)
    
    forecasts = {}
    
    # Iteratively compute A^k and the resulting expectation
    a_pow_k = identity.copy()
    
    for k in range(1, h + 1):
        # Update to A^k
        a_pow_k = a_pow_k @ a
        
        # Compute expectation vector: (I - A^k) * (I - A)^-1 * c + A^k * X0
        term_intercept = (identity - a_pow_k) @ i_a_inv_c
        term_state = a_pow_k @ x0
        
        expectation_vector = term_intercept + term_state
        
        # Extract the first element (y_{t+k})
        forecasts[f"t+{k}"] = float(expectation_vector[0, 0])
        
    return forecasts

# ---------------------------------------------------------------------------
# ARMA_variance
# ---------------------------------------------------------------------------

def arma_variance(h: int = 1, a: np.ndarray = None, b: np.ndarray = None, sigma2: float = 1.0):
    """
    Compute the long-term forecast variance of an ARMA model.
    """
    # Dimension
    dim = a.shape[0]
    
    # Outer product of the shock selector b: (dim x dim)
    # Ensure b is handled as a column vector for the outer product
    b = np.asarray(b).reshape(-1, 1)
    bb = b @ b.T
    
    # Initialize the sum of matrices
    # The R code starts with Identity. We follow that logic for 1:1 translation.
    current_sum = np.eye(dim)
    
    # List to store the [1,1] element of the cumulative sum at each step
    variances = {}
    variances["t+1"] = float(current_sum[0, 0] * sigma2)
    
    # Transition matrix power A^j
    a_pow_j = np.eye(dim)
    
    if h > 1:
        for i in range(1, h):
            # Update A^j to A^{i}
            a_pow_j = a @ a_pow_j
            
            # Add the variance contributed by the shock at step i
            # Term: A^j * (b * b^T) * (A^j)^T
            term = a_pow_j @ bb @ a_pow_j.T
            current_sum = current_sum + term
            
            # Store the [1,1] element of the cumulative sum
            variances[f"t+{(i + 1)}"] = float(current_sum[0, 0] * sigma2)
            
    return variances

# ---------------------------------------------------------------------------
# ARMA_covariance
# ---------------------------------------------------------------------------

def arma_covariance(h: int, k: int, a: np.ndarray, b: np.ndarray, sigma2: float = 1.0) -> np.ndarray:
    """
    Compute the ARMA conditional variance / covariance matrix.
    
    Args:
        h: Steps ahead for the first state.
        k: Steps ahead for the second state.
        a: Companion matrix (p+q x p+q).
        b: Innovation impact vector (p+q x 1).
        sigma2: Variance of the residuals.
        
    Returns:
        The (p+q x p+q) covariance matrix.
    """
    # Standardize b to column vector
    b = np.asarray(b).reshape(-1, 1)
    bb = b @ b.T
    
    # Initialize the sum matrix
    cv_x = np.zeros_like(a, dtype=float)
    
    # The number of common shocks is min(h, k)
    min_hk = min(h, k)
    
    # Direct translation of the R loop logic
    for j in range(min_hk):
        # R logic: h-1-j and k-1-j
        # Note: np.linalg.matrix_power handles 0 as Identity
        ah_j = np.linalg.matrix_power(a, h - 1 - j)
        ak_j = np.linalg.matrix_power(a, k - 1 - j)
        
        term = ah_j @ bb @ ak_j.T
        cv_x += term
        
    return sigma2 * cv_x

# ---------------------------------------------------------------------------
# ARMA_next_step
# ---------------------------------------------------------------------------

def arma_next_step(h: int = 1, x0: np.ndarray = None, a: np.ndarray = None, 
                   b: np.ndarray = None, intercept: float = 0.0, eps: np.ndarray = 0.0):
    """
    Compute the next-step state vector of an ARMA process.
    
    Args:
        h: Number of steps ahead.
        x0: Initial state vector of length (p + q).
        a: Companion matrix of dimension (p+q) x (p+q).
        b: Shock selector vector of length (p + q).
        intercept: Scalar intercept parameter.
        eps: Scalar 0 or vector of residuals of length h.
        
    Returns:
        The state vector after h steps as a (p+q, 1) numpy array.
    """
    # Standardize inputs to numpy arrays and column vectors
    x_t = np.asarray(x0).reshape(-1, 1)
    a = np.asarray(a)
    b = np.asarray(b).reshape(-1, 1)
    
    # Handle epsilon (shocks)
    # If eps is a scalar 0, create a vector of zeros of length h
    if np.isscalar(eps) and eps == 0:
        eps_vec = np.zeros(h)
    else:
        eps_vec = np.atleast_1d(eps)
        
    # Validations
    if len(eps_vec) != h:
        raise ValueError(f"The length of `eps` ({len(eps_vec)}) must be equal to `h` ({h})!")
    
    if x_t.shape[0] != a.shape[1]:
        raise ValueError(f"The length of `X0` ({x_t.shape[0]}) must match the columns of `A` ({a.shape[1]})!")

    # Forecasting/Simulation loop
    for step_idx in range(h):
        # 1. Transition: A * x_{t-1}
        # 2. Add Shock: b * eps[step]
        # Using @ for matrix multiplication and * for scalar-vector multiplication
        x_t = (a @ x_t) + (b * eps_vec[step_idx])
        
        # 3. Add Intercept to the first element (the process observation)
        x_t[0, 0] += intercept
        
    return x_t

# ---------------------------------------------------------------------------
# ARMA_filter
# ---------------------------------------------------------------------------

def arma_filter(x, a, b, intercept=0.0):
    """
    Filter the time-series to compute fitted values.
    
    Args:
        x: Numeric array of observations.
        a: Companion matrix (must have ar_order and ma_order attributes).
        b: Shock selector vector.
        intercept: Scalar intercept.
        
    Returns:
        np.ndarray: Fitted values (x_hat).
    """
    # Extract orders from matrix attributes (set during arma_companion_matrix)
    p = getattr(a, 'ar_order', 0)
    q = getattr(a, 'ma_order', 0)
    
    n = len(x)
    dim = a.shape[0]
    
    # Initialization
    x_hat = np.zeros(n)
    # State vector X_t initialized to zeros
    x_state = np.zeros((dim, 1))
    
    # Standardize b to column vector
    b = np.asarray(b).reshape(-1, 1)
    
    # Filtering Loop
    for t in range(n):
        # 1. Prediction Step (One-step ahead)
        # Based on ARMA_next_step logic: Predicted state = A * X_{t-1}
        x_pred_state = a @ x_state
        
        # 2. Predicted Observation (Fitted Value)
        # The first element of the state vector + intercept
        fitted_value = x_pred_state[0, 0] + intercept
        x_hat[t] = fitted_value
        
        # 3. Calculate Innovation (Residual)
        # epsilon_t = observed - predicted
        epsilon_t = x[t] - fitted_value
        
        # 4. Update Step
        # The current state X_t is the predicted state + the impact of the shock
        # Note: We add the intercept to the first element to align with ARMA_next_step
        x_state = x_pred_state
        x_state[0, 0] += intercept
        x_state = x_state + (b * epsilon_t)
        
    return x_hat

# ---------------------------------------------------------------------------
# ARMA_forecast
# ---------------------------------------------------------------------------

def arma_forecast(h: int, x0: np.ndarray, a: np.ndarray, b: np.ndarray, intercept: float = 0.0):
    """
    Fast ARMA state-space h-step forecast and trajectory history.
    
    Args:
        h: Integer, steps ahead.
        x0: Initial state vector of length p+q.
        a: (p+q) x (p+q) companion matrix.
        b: Numeric vector length p+q (shock selector).
        intercept: Scalar intercept.
        
    Returns:
        pd.DataFrame: A single-row DataFrame containing the final step and forecast,
                     with the full trajectory nested in a 'weights' column.
    """
    # Ensure inputs are correct shapes
    a = np.asmatrix(a)
    x_t = np.array(x0).reshape(-1, 1)
    dim = a.shape[0]
    
    # Storage for the full trajectory (h rows, dim columns)
    # We use a list of dictionaries or a 2D array for speed
    history = np.zeros((h, dim))
    
    # Forecast loop (Expectation: shocks = 0)
    current_x = x_t.copy()
    for i in range(h):
        # State transition
        current_x = a @ current_x
        # Add intercept to the first element (y_t component)
        current_x[0, 0] += intercept
        # Record the state (transpose to store as a row)
        history[i, :] = current_x.flatten()
        
    # Convert trajectory to DataFrame (mimicking dplyr::bind_rows)
    df_history = pd.DataFrame(history)
    # Add step column (1-based to match R)
    df_history['step'] = np.arange(1, h + 1)
    
    # Prepare the final output (tail(df_tT, 1) in R)
    # R extracts columns [5, 1] -> step and the first state component
    last_row = df_history.iloc[[-1]].copy()
    
    # Reorganize to match R: df_T <- tail(...)[:, c(5, 1)]
    # Python indices: 'step' is the last column, 0 is the first column
    df_t = last_row[['step', 0]].rename(columns={0: 'forecast'})
    
    # Nest the full history as 'weights'
    df_t['weights'] = [df_history]
    
    return df_t


# ---------------------------------------------------------------------------
# AR_variance  (closed-form for p ≤ 4, companion matrix otherwise)
# ---------------------------------------------------------------------------

def ar_variance(phi: np.ndarray, sigma2: float = 1.0) -> float:
    """
    Compute the long-term stationary variance of an AR model.
    
    Args:
        phi: Numeric array of AR parameters [phi1, phi2, ...].
        sigma2: Variance of the residuals (white noise).
        
    Returns:
        float: The stationary variance of the process.
    """
    phi = np.atleast_1d(phi)
    ar_order = len(phi)
    
    if ar_order == 0:
        return float(sigma2)
        
    if ar_order == 1:
        # AR(1): sigma^2 / (1 - phi1^2)
        return float(sigma2 / (1 - phi[0]**2))
    
    elif ar_order == 2:
        # AR(2) closed form
        p1, p2 = phi[0], phi[1]
        num = 1 - p2
        den = (1 - p2) * (1 - p1**2 - p2**2) - 2 * (p1**2) * p2
        return float(sigma2 * (num / den))
    
    elif ar_order == 3:
        # AR(3) recursive substitution logic
        p1, p2, p3 = phi[0], phi[1], phi[2]
        phi_t1 = (p1 + p2*p3) / (1 - p2 - p3**2 - p1*p3)
        phi_t2 = (p1 + p3) * phi_t1 + p2
        phi_t3 = p1*phi_t2 + p2*phi_t1 + p3
        phi_t0 = 1 / (1 - p1*phi_t1 - p2*phi_t2 - p3*phi_t3)
        return float(phi_t0 * sigma2)
    
    elif ar_order == 4:
        # AR(4) recursive substitution logic
        p1, p2, p3, p4 = phi[0], phi[1], phi[2], phi[3]
        phi_1_alt = (p1 + p3) / (1 - p4)
        phi_0_alt = p2 / (1 - p4)
        
        psi_1 = (phi_1_alt*p3 + p1*phi_0_alt*p4 + p1 + p3*p4)
        psi_1 /= (1 - phi_1_alt*(p3 + p1*p4) - p2*(1 + p4) - p4**2)
        
        psi_2 = phi_1_alt * psi_1 + phi_0_alt
        psi_3 = p1*psi_2 + p2*psi_1 + p4*psi_1 + p3
        psi_4 = p1*psi_3 + p2*psi_2 + p3*psi_1 + p4
        
        var_factor = 1 / (1 - p1*psi_1 - p2*psi_2 - p3*psi_3 - p4*psi_4)
        return float(sigma2 * var_factor)
    
    else:
        # Fallback for p > 4: Use state-space companion matrix
        # Note: We fix the positional argument bug from the R source here.
        # We assume h=100 is sufficient for convergence to stationary variance,
        # or we should use a Lyapunov solver (e.g., scipy.linalg.solve_discrete_lyapunov).
        a_matrix = arma_companion_matrix(phi=phi)
        b_vec = arma_vector_b(ar_order, 0)
        
        # We call the dictionary-based arma_variance and take the last step
        # as a proxy for the long-term variance.
        h_steps = 100 
        var_dict = arma_variance(h=h_steps, a=a_matrix, b=b_vec, sigma2=sigma2)
        return var_dict[f"t+{h_steps}"]


# ---------------------------------------------------------------------------
# AR_params_to_zeta  (step-down Durbin-Levinson + atanh)
# ---------------------------------------------------------------------------

def ar_params_to_zeta(phi: np.ndarray) -> np.ndarray:
    """
    Transform constrained AR parameters (phi) to unconstrained parameters (zeta).
    
    Args:
        phi: Array of AR coefficients of length p.
        
    Returns:
        Array of unconstrained parameters of length p.
    """
    phi = np.atleast_1d(phi).astype(float)
    p = len(phi)
    
    if p == 0:
        return np.array([], dtype=float)
    
    # Storage for intermediate AR vectors phi^(k)
    # We only need the current and previous, but for clarity 
    # we can use a list or iterate in-place.
    current_phi = phi.copy()
    kappa = np.zeros(p)
    
    # Iterate backwards from p down to 1
    # Python range(start, stop, step): p-1 down to 0
    for k_idx in range(p - 1, -1, -1):
        # Extract the last reflection coefficient
        k_val = current_phi[k_idx]
        
        # Tiny clipping for numerical safety
        if abs(k_val) >= 1.0:
            k_val = np.sign(k_val) * (1.0 - 1e-12)
        
        kappa[k_idx] = k_val
        
        # Back out phi^(k-1) if k > 0
        if k_idx > 0:
            # Formula: (phi_i^(k) + kappa_k * phi_{k-i}^(k)) / (1 - kappa_k^2)
            # Use slicing for vectorization
            num = current_phi[:k_idx] + k_val * np.flip(current_phi[:k_idx])
            den = 1.0 - k_val**2
            current_phi = num / den
            
    # Map to (-inf, inf)
    zeta = np.arctanh(kappa)
    return zeta

# ---------------------------------------------------------------------------
# ARMA_params_to_zeta  (step-down Durbin-Levinson + atanh)
# ---------------------------------------------------------------------------

def arma_params_to_zeta(phi=None, theta=None):
    """
    Transform ARMA parameters (phi, theta) to unconstrained space (zeta).
    
    Args:
        phi: Array-like, AR parameters.
        theta: Array-like, MA parameters.
        
    Returns:
        dict: Containing 'zeta_phi' and 'zeta_theta' arrays.
    """
    results = {
        "zeta_phi": np.array([]),
        "zeta_theta": np.array([])
    }

    # Process AR parameters
    # Following R logic: skip if None or if the first element is 0
    if phi is not None:
        phi_arr = np.atleast_1d(phi)
        if len(phi_arr) > 0 and phi_arr[0] != 0:
            results["zeta_phi"] = ar_params_to_zeta(phi_arr)

    # Process MA parameters
    if theta is not None:
        theta_arr = np.atleast_1d(theta)
        if len(theta_arr) > 0 and theta_arr[0] != 0:
            results["zeta_theta"] = ar_params_to_zeta(theta_arr)

    return results

# ---------------------------------------------------------------------------
# AR_params_to_phi  (tanh + step-up Durbin-Levinson + Jacobian)
# ---------------------------------------------------------------------------

def ar_params_to_phi(zeta: np.ndarray):
    """
    Transform unconstrained parameters (zeta) to constrained AR parameters (phi).
    Also returns the Jacobian matrices for the transformation.
    
    Args:
        zeta: Array of unconstrained parameters.
        
    Returns:
        dict: Containing 'phi', 'kappa', 'J_phi_kappa', and 'J_phi_zeta'.
    """
    zeta = np.atleast_1d(zeta).astype(float)
    p = len(zeta)
    
    # 1. Map zeta to reflection coefficients kappa (partial autocorrelations)
    kappa = np.tanh(zeta)
    
    # Storage for intermediate phi vectors and sensitivity matrices
    phis = [None] * p
    G = [None] * p # G[k] is (k+1)x(k+1) matrix of d phi^(k+1) / d kappa[0:k+1]

    # Base case: k = 0 (order 1 in R)
    phis[0] = np.array([kappa[0]])
    G[0] = np.array([[1.0]])

    # Recursive build
    for k in range(1, p):
        phi_prev = phis[k-1]
        phi_new = np.zeros(k + 1)
        
        # Build new phi coefficients for order k+1
        # phi_new[j] = phi_prev[j] - kappa[k] * phi_prev[k-j-1]
        for j in range(k):
            phi_new[j] = phi_prev[j] - kappa[k] * phi_prev[k - j - 1]
        phi_new[k] = kappa[k]
        phis[k] = phi_new

        # Build sensitivity matrix G[k] (size (k+1) x (k+1))
        G_prev = G[k-1]
        Gk = np.zeros((k + 1, k + 1))

        # 1. Columns m = 0...k-1 (propagation of derivatives from previous order)
        for m in range(k):
            for j in range(k):
                # Gk[j, m] = d phi_new[j] / d kappa[m]
                Gk[j, m] = G_prev[j, m] - kappa[k] * G_prev[k - j - 1, m]
        
        # 2. Column m = k (derivative with respect to the newest kappa)
        for j in range(k):
            # Gk[j, k] = d phi_new[j] / d kappa[k] = -phi_prev[k-j-1]
            Gk[j, k] = -phi_prev[k - j - 1]
            
        # 3. Last row: d phi_new[k] / d kappa[...]
        # phi_new[k] is just kappa[k], so d/d kappa[k] = 1, others = 0
        Gk[k, k] = 1.0
        
        G[k] = Gk

    phi = phis[p - 1]
    
    # Jacobian from phi to kappa
    j_phi_kappa = G[p - 1]
    
    # Jacobian from phi to zeta (chain rule)
    # d phi / d zeta = (d phi / d kappa) * (d kappa / d zeta)
    # d kappa / d zeta = sech^2(zeta) = 1 - tanh^2(zeta)
    sech2 = 1.0 - kappa**2
    j_phi_zeta = j_phi_kappa @ np.diag(sech2)
    
    return {
        "phi": phi,
        "kappa": kappa,
        "J_phi_kappa": j_phi_kappa,
        "J_phi_zeta": j_phi_zeta
    }

# ---------------------------------------------------------------------------
# ARMA_params_to_phi  (tanh + step-up Durbin-Levinson + Jacobian)
# ---------------------------------------------------------------------------

def arma_params_to_phi(zeta_phi=None, zeta_theta=None):
    """
    Convert unconstrained ARMA parameters to phi/theta and assemble the 
    block-diagonal Jacobian matrix.
    """
    # Initialize variables
    p, q = 0, 0
    phi, theta = np.array([]), np.array([])
    j_phi_zeta, j_theta_zeta = None, None

    # 1. Process AR Parameters
    if zeta_phi is not None:
        zeta_phi = np.atleast_1d(zeta_phi)
        if zeta_phi.size > 0:
            p = zeta_phi.size
            res_phi = ar_params_to_phi(zeta_phi)
            phi = res_phi['phi']
            j_phi_zeta = res_phi['J_phi_zeta']

    # 2. Process MA Parameters
    if zeta_theta is not None:
        zeta_theta = np.atleast_1d(zeta_theta)
        if zeta_theta.size > 0:
            q = zeta_theta.size
            res_theta = ar_params_to_phi(zeta_theta)
            theta = res_theta['phi'] # MA uses same logic as AR
            j_theta_zeta = res_theta['J_phi_zeta']

    # 3. Assemble the Global Jacobian Matrix (p+q, p+q)
    j_total = np.zeros((p + q, p + q))
    
    if p > 0:
        j_total[:p, :p] = j_phi_zeta
        
    if q > 0:
        # Bottom-right block starts at index p
        j_total[p : p + q, p : p + q] = j_theta_zeta

    return {
        "phi": phi,
        "theta": theta,
        "J": j_total
    }

# ---------------------------------------------------------------------------
# Loss functions (constrained parametrisation)
# ---------------------------------------------------------------------------

def arma_loss_log_lik(params, p=0, q=0, y=None, per_obs=False):
    """
    Compute the Gaussian log-likelihood for an ARMA(p, q) model.
    
    Args:
        params: Array-like or dict containing unconstrained parameters.
        p (int): AR order.
        q (int): MA order.
        y (np.ndarray): The observed time series.
        per_obs (bool): If True, returns an array of log-likelihoods.
        
    Returns:
        float or np.ndarray: Total or pointwise log-likelihood.
    """
    # 1. Handle Intercept
    intercept = 0.0
    # In Python, we check if 'intercept' is a key in a dict or index in a named series
    if isinstance(params, dict):
        intercept = params.get("intercept", 0.0)
        zeta_phi = params.get("zeta_phi", np.array([]))
        zeta_theta = params.get("zeta_theta", np.array([]))
    else:
        # Assuming params is a flat array: [zeta_phi..., zeta_theta...]
        # Note: Intercept handling in the R source suggests it might be appended.
        zeta_phi = params[0:p] if p > 0 else np.array([])
        zeta_theta = params[p:p+q] if q > 0 else np.array([])
        # If there's an extra element for intercept
        if len(params) > (p + q):
            intercept = params[-1]

    # 2. Convert parameters to constrained space
    # Uses the previously translated arma_params_to_phi
    phi_theta = arma_params_to_phi(zeta_phi, zeta_theta)
    
    # 3. Setup state-space matrices
    # Uses previously translated arma_companion_matrix and arma_vector_b
    a_matrix = arma_companion_matrix(phi_theta['phi'], phi_theta['theta'])
    b_vec = arma_vector_b(p, q)
    
    # 4. Generate fitted values
    # Uses previously translated arma_filter
    y_hat = arma_filter(y, a_matrix, b_vec, intercept=intercept)
    
    # 5. Compute residuals and estimate sigma
    # R: idx_excluded <- 1:max(c(p, q))
    k = max(p, q)
    eps_hat = (y - y_hat)[k:] # Remove the first k elements
    
    # sigma_hat = sqrt(mean(eps_hat^2))
    sigma_hat = np.sqrt(np.mean(eps_hat**2))
    
    # 6. Calculate Log-Likelihoods
    # log(f(z) / sigma) = log(f(z)) - log(sigma)
    # where z = eps / sigma
    z_hat = eps_hat / sigma_hat
    
    # Using scipy.stats.norm.logpdf for numerical stability
    # logpdf(z) already includes the -0.5 * log(2*pi) part
    log_liks = norm.logpdf(z_hat) - np.log(sigma_hat)
    
    if per_obs:
        return log_liks
    else:
        return np.sum(log_liks)

# ---------------------------------------------------------------------------
# ARMA_loss_CSS  (Conditional Sum of Squares)
# ---------------------------------------------------------------------------

def arma_loss_css(params, p=0, q=0, y=None):
    """
    Compute the Conditional Sum of Squares (CSS) for an ARMA(p, q) model.
    
    Args:
        params: Array-like or dict containing unconstrained parameters.
        p (int): AR order.
        q (int): MA order.
        y (np.ndarray): The observed time series.
        
    Returns:
        float: The sum of squared residuals.
    """
    # 1. Handle Intercept and Parameter Extraction
    intercept = 0.0
    if isinstance(params, dict):
        intercept = params.get("intercept", 0.0)
        zeta_phi = params.get("zeta_phi", np.array([]))
        zeta_theta = params.get("zeta_theta", np.array([]))
    else:
        # Assuming params is a flat array: [zeta_phi..., zeta_theta..., intercept]
        # logic matches R's params[1:p] and params[(p+1):(p+q)]
        zeta_phi = params[0:p] if p > 0 else np.array([])
        zeta_theta = params[p:p+q] if q > 0 else np.array([])
        # Check for intercept at the end if the vector length matches
        if len(params) > (p + q):
            intercept = params[-1]

    # 2. Convert parameters to constrained space
    phi_theta = arma_params_to_phi(zeta_phi, zeta_theta)
    
    # 3. Setup state-space matrices
    a_matrix = arma_companion_matrix(phi_theta['phi'], phi_theta['theta'])
    b_vec = arma_vector_b(p, q)
    
    # 4. Generate fitted values using the recursive filter
    y_hat = arma_filter(y, a_matrix, b_vec, intercept=intercept)
    
    # 5. Compute residuals and truncate
    # R: idx_excluded <- 1:max(c(p, q))
    k = max(p, q)
    eps_hat = (y - y_hat)[k:] 
    
    # 6. Return Sum of Squared Errors (SSE)
    return float(np.sum(eps_hat**2))

# ---------------------------------------------------------------------------
# ARMA_fit
# ---------------------------------------------------------------------------

def arma_fit(y, ar_order=1, ma_order=0, method="ML"):
    """
    Fit ARMA model parameters and compute robust standard errors.
    """
    n = len(y)
    pq = ar_order + ma_order
    
    # 1. Initialize unconstrained parameters zeta
    # R: runif(arOrder, -0.1, 0.1)
    zeta_init = np.random.uniform(-0.1, 0.1, pq)
    
    # 2. Define Loss Function
    # Note: scipy.minimize always minimizes. 
    # ML already returns -log_lik in our previous translation's context if we want to minimize.
    if method == "CSS":
        def loss_func(params):
            return arma_loss_css(params, p=ar_order, q=ma_order, y=y)
    elif method == "ML":
        def loss_func(params):
            # Minimize negative log-likelihood
            return -arma_loss_log_lik(params, p=ar_order, q=ma_order, y=y, per_obs=False)
    else:
        raise ValueError("Method must be 'CSS' or 'ML'")

    # 3. Optimization
    # BFGS is a good choice for smooth likelihood surfaces
    opt = minimize(loss_func, zeta_init, method='BFGS')
    zeta_qmle = opt.x
    
    # 4. Parameter Conversion & Jacobian
    # Extract zeta_phi and zeta_theta for the conversion function
    z_phi = zeta_qmle[:ar_order] if ar_order > 0 else None
    z_theta = zeta_qmle[ar_order:] if ma_order > 0 else None
    
    res = arma_params_to_phi(z_phi, z_theta)
    
    # 5. Numerical Derivatives for Inference
    # We need the derivatives of the LOG-LIKELIHOOD (not the negative loss)
    def loglik_wrapper(p_vec, per_obs=False):
        return arma_loss_log_lik(p_vec, p=ar_order, q=ma_order, y=y, per_obs=per_obs)

    # Hessian of Log-Likelihood (at the optimum)
    hess_func = nd.Hessian(lambda p: loglik_wrapper(p, per_obs=False))
    H = hess_func(zeta_qmle)
    
    # Jacobian (Score) per observation
    jac_func = nd.Jacobian(lambda p: loglik_wrapper(p, per_obs=True))
    S = jac_func(zeta_qmle) # Matrix of (n-k) x pq
    
    # B matrix (Outer product of scores)
    # B = sum(score_t * score_t^T)
    B = S.T @ S
    
    # 6. Covariance Matrices (Unconstrained Space)
    # V_star = solve(-H)
    try:
        V_star = np.linalg.inv(-H)
    except np.linalg.LinAlgError:
        V_star = np.full((pq, pq), np.nan)

    # 7. Covariance Matrices (Constrained Space via Delta Method)
    # V_orig = J * V_star * J^T
    J = res['J']
    V_orig = J @ V_star @ J.T
    
    # 8. Robust (Sandwich) Covariance
    # V_rob = V_star * B * V_star
    V_rob_star = V_star @ B @ V_star
    V_rob_orig = J @ V_rob_star @ J.T
    
    # 9. Extract Standard Errors
    std_errors_zeta = np.sqrt(np.diag(V_star))
    std_errors_phi = np.sqrt(np.diag(V_orig))
    
    std_errors_rob_zeta = np.sqrt(np.diag(V_rob_star))
    std_errors_rob_phi = np.sqrt(np.diag(V_rob_orig))

    # Assemble Results
    res.update({
        'loss': opt.fun,
        'method': method,
        'H': H,
        'vcov': {
            'V_star': V_star,
            'V_orig': V_orig,
            'V_rob_star': V_rob_star,
            'V_rob_orig': V_rob_orig
        },
        'std_errors': {
            'zeta': std_errors_zeta,
            'phi': std_errors_phi
        },
        'std_errors_rob': {
            'zeta': std_errors_rob_zeta,
            'phi': std_errors_rob_phi
        }
    })
    
    return res