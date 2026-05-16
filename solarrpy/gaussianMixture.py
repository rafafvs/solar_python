import numpy as np
import pandas as pd
import warnings
from scipy.stats import norm
import numdifftools as nd
from .gaussianMixture_internals import gm_moments, gm_fit, gm_fit_moments_match

class GaussianMixtureModel:
    """
    Gaussian Mixture Model for time series classification and moment matching.
    """
    _version = "1.0.0"
    def __init__(self, components: int = 2, maxit: int = 5000, maxrestarts: int = 500, abstol: float = 1e-08):
        # Control parameters
        self.components = components
        self.maxit = maxit
        self.maxrestarts = maxrestarts
        self.abstol = abstol
        
        # Primary parameters
        self._means = np.linspace(-3.0, 3.0, num=components)
        self._sd = np.ones(components, dtype=float)
        init_p = np.full(shape=components, fill_value=1.0 / components)
        self._p = init_p / np.sum(init_p)
        
        # State placeholders (Mapped from R6 private list)
        self._x = None
        self._loglik = np.nan
        self._fitted = None
        
        self._vcov = None
        self._std_means = np.full(components, np.nan)
        self._std_sd = np.full(components, np.nan)
        self._std_p = np.full(components, np.nan)
        
        self._use_empiric = False
        self._empiric_params = {}
        self._responsibilities = None
        
        # Initializing moments dictionary structure
        self._moments = {
            'm1': np.nan, 'm2': np.nan, 'm3': np.nan, 'm4': np.nan,
            'mean': np.nan, 'variance': np.nan, 'skewness': np.nan, 
            'kurtosis': np.nan, 'nobs': np.nan
        }

    def logLik(self, x: np.ndarray, params: dict = None, per_obs: bool = False) -> np.ndarray | float:
        """
        Compute the log-likelihood of the dataset.
        """
        if params is None:
            params = self.coefficients
            
        # Ensure array types
        x = np.asarray(x, dtype=float)
        means = np.asarray(params['means'], dtype=float)
        sd = np.asarray(params['sd'], dtype=float)
        p = np.asarray(params['p'], dtype=float)
        
        # Vectorized density calculation using broadcasting
        # x is reshaped to (N, 1) to broadcast against (K,) parameter arrays
        component_densities = norm.pdf(x[:, np.newaxis], loc=means, scale=sd)
        
        # Weight by probabilities and sum across components (axis=1)
        mix_density = np.sum(p * component_densities, axis=1)
        
        # Protect against log(0) mathematically underflowing to -inf
        # Using 1e-300 as a safe lower bound for 64-bit floats
        mix_density = np.clip(mix_density, a_min=1e-300, a_max=None)
        
        # Compute log-likelihood per observation
        log_likelihood = np.log(mix_density)
        
        if not per_obs:
            # np.nansum mirrors R's sum(..., na.rm = TRUE)
            log_likelihood = np.nansum(log_likelihood)
            
        return log_likelihood

    def E_step(self, x: np.ndarray = None, params: dict = None) -> pd.DataFrame:
        """
        Compute the posterior probabilities (responsibilities).
        """
        # State fallback
        if x is None:
            x = self._x
        if params is None:
            params = self.coefficients
            
        x = np.asarray(x, dtype=float)
        means = np.asarray(params['means'], dtype=float)
        sd = np.asarray(params['sd'], dtype=float)
        p = np.asarray(params['p'], dtype=float)
        
        # 1. Numerator: Weighted likelihoods via broadcasting (N x K)
        component_densities = norm.pdf(x[:, np.newaxis], loc=means, scale=sd)
        weighted_densities = p * component_densities
        
        # 2. NaN Coercion: Replicate R's ifelse(is.na(x)|is.nan(x), 0, x)
        weighted_densities = np.nan_to_num(weighted_densities, nan=0.0)
        
        # 3. Denominator: Row sums
        # keepdims=True ensures the shape is (N, 1) to allow proper broadcasting during division
        row_sums = np.sum(weighted_densities, axis=1, keepdims=True)
        
        # 4. Normalization: Safely divide to get posterior probabilities
        # Prevent 0/0 division errors if an observation is extremely far from all components
        with np.errstate(divide='ignore', invalid='ignore'):
            responsibilities = np.where(
                row_sums > 0, 
                weighted_densities / row_sums, 
                0.0  # Assign 0 responsibility if the probability space collapsed
            )
            
        # 5. Output formatting
        cols = [f"B{k}" for k in range(self.components)]
        return pd.DataFrame(responsibilities, columns=cols)

    def classify(self, x: np.ndarray = None) -> pd.DataFrame:
        """
        Classify the time series into the components with highest likelihood.
        """
        if x is None:
            x = self._x
            
        x = np.asarray(x, dtype=float)
        params = self.coefficients
        n = len(x)
        
        # E-step: posterior probabilities (Returns a DataFrame)
        resp_df = self.E_step(x)
        resp = resp_df.values
        
        # 1. Vectorized Classification
        # np.argmax returns 0-based indices. 
        # We add 1 to match R's 1-based indexing for class labels (1, 2, ..., K).
        class_idx = np.argmax(resp, axis=1)
        classification = class_idx
        
        # 2. Vectorized Uncertainty
        uncertainty = 1.0 - np.max(resp, axis=1)
        
        # 3. Vectorized B_hat (One-Hot Encoding)
        B_hat = np.zeros((n, self.components), dtype=float)
        B_hat[np.arange(n), class_idx] = 1.0
        
        # 4. Vectorized x_hat
        # x[:, np.newaxis] broadcasts the (N,) array to multiply cleanly against the (N, K) B_hat
        x_hat = B_hat * x[:, np.newaxis]
        
        # 5. Vectorized z_hat
        means = np.asarray(params['means'], dtype=float)
        sd = np.asarray(params['sd'], dtype=float)
        z_hat = (x_hat - means) / sd * B_hat
        
        # 6. Construct output DataFrame dynamically
        out_dict = {'classification': classification}
        
        for k in range(self.components):
            out_dict[f'B{k}'] = B_hat[:, k]
            
        for k in range(self.components):
            out_dict[f'x{k}'] = x_hat[:, k]
            
        for k in range(self.components):
            out_dict[f'z{k}'] = z_hat[:, k]
            
        out_dict['uncertainty'] = uncertainty
        
        return pd.DataFrame(out_dict)

    def fit(self, x: np.ndarray, mu_target: float = None, var_target: float = None) -> None:
        """
        Fit the mixture parameters using EM, optionally matching target moments.
        """
        x = np.asarray(x, dtype=float)
        n = len(x)
        
        # Store full series in state
        self._x = x
        x_fit = x[~np.isnan(x)]

        if len(x_fit) < self.components:
            raise ValueError("Not enough valid observations to fit the requested number of components.")
            
        # 2. Unconstrained Fit
        # Assuming gm_fit is accessible in the module namespace
        clust = gm_fit(x_fit, components=self.components, maxit=self.maxit, maxrestarts=self.maxrestarts)
        
        # 3. Moments Matching (Constrained Fit)
        if self.components == 2 and (mu_target is not None or var_target is not None):
            # Flatten parameters and drop p2 (index 1 of p array) to match expected start_params signature
            # R: unlist(purrr::flatten(clust)[-6])
            start_params = np.array([
                clust['means'][0], clust['means'][1],
                clust['sd'][0], clust['sd'][1],
                clust['p'][0]
            ], dtype=float)
            
            clust = gm_fit_moments_match(
                x=x_fit, 
                start_params=start_params, 
                mu_target=mu_target, 
                var_target=var_target, 
                maxit=self.maxit, 
                abstol=self.abstol
            )
            
        # 4. Update Internal State
        self._means = clust['means']
        self._sd = clust['sd']
        self._p = clust['p']
        
        # 5. Trigger Update Cascade
        # These methods must be implemented on the class
        self._responsibilities = self.E_step(x)
        self.update_logLik()
        self._fitted = self.classify(x)
        self.update_empiric_parameters()
        self.Hessian()

    def update(self, means: np.ndarray = None, sd: np.ndarray = None, p: np.ndarray = None) -> None:
        """
        Update the parameters inside the object.
        """
        # Update mean parameters
        if means is not None:
            self._means = np.asarray(means, dtype=float).flatten()
            
        # Update standard deviation parameters
        if sd is not None:
            self._sd = np.asarray(sd, dtype=float).flatten()
            
        # Update probability parameters
        if p is not None:
            self._p = np.asarray(p, dtype=float).flatten()
            
        # Note: Parameter naming (e.g., 'mu1', 'mu2') is omitted here because 
        # standard numpy arrays do not hold element-wise string names. 
        # Name binding is handled dynamically in the `summary` property.
        
        # Constraint check:
        if len(self._means) != self.components:
            raise ValueError("Length mismatch between means and components")

    def update_logLik(self) -> None:
        """
        Update the log-likelihood with the current parameters.
        """
        # Failsafe: Cannot compute likelihood if data hasn't been ingested
        if self._x is None:
            self._loglik = np.nan
            return
            
        x = self._x
    
        # Update log-likelihood slot
        self._loglik = self.logLik(x[~np.isnan(x)])

    def update_empiric_parameters(self) -> None:
        """
        Compute the parameters on the classified time series.
        """
        # Failsafe: Ensure data exists before operating
        if self._fitted is None or self._x is None:
            return
            
        df = self._fitted.copy()
        df['x'] = self._x
        
        # Group by the hard classification labels.
        # pandas .groupby() automatically sorts by the group keys by default,
        # which matches dplyr::arrange(classification).
        grouped = df.groupby('classification')['x']
        
        # Compute empirical parameters
        # Pandas mean() and std() automatically drop NaNs (acting like na.rm=TRUE)
        mu = grouped.mean()
        sd = grouped.std(ddof=1)
        p = grouped.count() / len(df)
        
        # Check for collapsed classes
        if len(mu) < self.components:
            warnings.warn(
                "Classification resulted in fewer classes than components! "
                "Falling back to theoretical parameters for empirical slots.",
                RuntimeWarning
            )
            self._empiric_params['means'] = self._means.copy()
            self._empiric_params['sd'] = self._sd.copy()
            self._empiric_params['p'] = self._p.copy()
        else:
            # Extract the computed values as numpy arrays
            self._empiric_params['means'] = mu.values
            self._empiric_params['sd'] = sd.values
            self._empiric_params['p'] = p.values

        # Note: Naming logic (names(private$empiric_params$means) <- ...) is omitted 
        # as we are handling vector naming dynamically in the summary method.

    def filter(self) -> None:
        """
        Update the responsibilities, log-likelihood, classification, 
        empirical parameters, and Hessian based on current parameters.
        """
        # 1. State Guard: Prevent execution if no data is loaded
        if self._x is None:
            raise RuntimeError("Cannot filter: No data (x) has been ingested into the model.")
            
        # 2. Update posterior probabilities
        self._responsibilities = self.E_step()
        
        # 3. Update log-likelihood
        self.update_logLik()
        
        # 4. Final classification of each component
        self._fitted = self.classify()
        
        # 5. Update empiric parameters
        self.update_empiric_parameters()
        
        # 6. Compute hessian matrix
        self.Hessian()

    def Hessian(self) -> None:
        """
        Compute the robust variance-covariance matrix and standard errors
        using numerical differentiation.
        """
        # Ensure we have data to evaluate
        if self._x is None:
            return

        params = self.coefficients
        K = self.components
        
        m = np.asarray(params['means'], dtype=float)
        s = np.asarray(params['sd'], dtype=float)
        p = np.asarray(params['p'], dtype=float)
        
        # FIX: Drop the strictly *last* probability rather than hardcoding index 1
        p_red = p[:-1]
        
        # Flatten parameters into a 1D array for numdifftools
        flat_params = np.concatenate([m, s, p_red])
        
        def loglik_wrapper(par_array, per_obs=False):
            # Reconstruct parameters by slicing (length is known deterministically)
            m_tmp = par_array[0:K]
            s_tmp = par_array[K:2*K]
            p_red_tmp = par_array[2*K:]
            
            # Recalculate the dropped K-th probability
            p_last = 1.0 - np.sum(p_red_tmp)
            p_tmp = np.append(p_red_tmp, p_last)
            
            par_dict = {'means': m_tmp, 'sd': s_tmp, 'p': p_tmp}
            
            # Filter data and compute log-likelihood
            x_filtered = self._x[~np.isnan(self._x)]
            return self.logLik(x_filtered, params=par_dict, per_obs=per_obs)
            
        # 1. Numerical Differentiation
        # nd.Hessian and nd.Jacobian return callable functions
        H_fun = nd.Hessian(lambda p: loglik_wrapper(p, per_obs=False))
        J_fun = nd.Jacobian(lambda p: loglik_wrapper(p, per_obs=True))
        
        H = H_fun(flat_params)
        J = J_fun(flat_params)
        
        # 2. Robust Variance-Covariance (Sandwich Estimator)
        # Vectorized equivalent of the R for-loop: S = J^T @ J
        S = J.T @ J
        
        # Safely invert the Hessian
        try:
            H_inv = np.linalg.inv(H)
        except np.linalg.LinAlgError:
            warnings.warn("Hessian is singular; using pseudo-inverse.", RuntimeWarning)
            H_inv = np.linalg.pinv(H)
            
        Cv = H_inv @ S @ H_inv
        
        # 3. Delta Method Transformation Matrix (Jac)
        n_red = len(flat_params)
        Jac = np.eye(n_red)
        last_row = np.zeros(n_red)
        
        # The derivatives of p_K with respect to p_1 ... p_{K-1} are all -1
        if K > 1:
            last_row[-(K-1):] = -1.0
            
        Jac = np.vstack((Jac, last_row))
        
        # 4. Full Covariance Matrix
        self._vcov = Jac @ Cv @ Jac.T
        
        # 5. Safely Extract Standard Errors
        try:
            std_errors = np.sqrt(np.diag(self._vcov))
        except np.linalg.LinAlgError:
            warnings.warn(
                f"Variance-covariance matrix has negative diagonal elements."
                "This may indicate the model did not converge to a proper local maximum, "
                "near-degenerate components, or a constrained optimum.",
                RuntimeWarning
            )
            std_errors = np.sqrt(np.abs(np.diag(self._vcov)))
        
        self._std_means = std_errors[0:K]
        self._std_sd = std_errors[K:2*K]
        self._std_p = std_errors[2*K:]

    def use_empiric_parameters(self) -> None:
        """
        Substitute the empiric parameters with EM parameters. 
        If evaluated again the EM parameters will be substituted back.
        """
        self._use_empiric = not self._use_empiric

    def __str__(self) -> str:
        """
        Readable string representation for print() and console output.
        Replicates the tabular console formatting of the R version.
        """
        # Helper function to format values and standard errors
        def format_params(vals, errs):
            formatted = []
            for v, e in zip(vals, errs):
                # \033[1;31m adds red ANSI color to the standard error
                formatted.append(f"{v:.3f} (\033[1;31m{e:.3f}\033[0m)")
            return " ".join(formatted)

        # Retrieve standard errors via the property
        se = self.std_errors
        
        means_str = format_params(self.means, se['means'])
        sd_str = format_params(self.sd, se['sd'])
        p_str = format_params(self.p, se['p'])

        # Build headers
        lbl_clean = f"{self.components} components"
        
        msg_0_clean = f"#################### GaussianMixture ({lbl_clean}) ####################"
        msg_0 = msg_0_clean + "\n"
        
        msg_1 = f"Means: {means_str} \n"
        msg_2 = f"Stddev: {sd_str} \n"
        msg_3 = f"Probs: {p_str} \n"
        
        # Safe Nobs calculation
        nobs = len(self._x) if self._x is not None else 0
        msg_4 = f"Nobs: {nobs}\n"
        msg_5 = f"Use empirical moments: {self.use_empiric}\n"
        
        # Safe log-likelihood formatting
        loglik_val = self.loglik if self.loglik is not None and not np.isnan(self.loglik) else 0.0
        msg_6 = f"Log-Likelihood: {loglik_val:.3f}\n"
        msg_7 = f"Version: {self._version}\n"

        # Dynamically generate the dashed line based on the header length
        line = "-" * len(msg_0_clean) + "\n"

        # Return the concatenated string instead of calling print()
        return f"{line}{msg_0}{line}{msg_1}{msg_2}{msg_3}{line}{msg_4}{msg_5}{line}{msg_6}{msg_7}"
        
    def __repr__(self) -> str:
        """
        Unambiguous developer representation of the object.
        """
        return (f"GaussianMixtureModel(components={self.components}, "
                f"maxit={self.maxit}, maxrestarts={self.maxrestarts}, "
                f"abstol={self.abstol})")


    @property
    def means(self) -> np.ndarray:
        if self._use_empiric and 'means' in self._empiric_params:
            return self._empiric_params['means']
        return self._means

    @property
    def sd(self) -> np.ndarray:
        if self._use_empiric and 'sd' in self._empiric_params:
            return self._empiric_params['sd']
        return self._sd

    @property
    def p(self) -> np.ndarray:
        if self._use_empiric and 'p' in self._empiric_params:
            return self._empiric_params['p']
        return self._p

    @property
    def coefficients(self) -> dict:
        return {'means': self.means, 'sd': self.sd, 'p': self.p}

    @property
    def use_empiric(self) -> bool:
        return self._use_empiric

    @property
    def std_errors(self) -> dict:
        return {'means': self._std_means, 'sd': self._std_sd, 'p': self._std_p}

    @property
    def model(self) -> pd.DataFrame:
        return pd.DataFrame({
            'means': self.means, 
            'sd': self.sd, 
            'p': self.p
        })

    @property
    def loglik(self) -> float:
        return self._loglik

    @property
    def fitted(self) -> pd.DataFrame:
        if self._fitted is None:
            return pd.DataFrame()
            
        df = self._fitted.copy()
        if 'uncertainty' in df.columns:
            # Pandas median() automatically drops NaNs, unlike R's median()
            limit_uncertainty = df['uncertainty'].median()
            df['uncertainty'] = np.where(df['uncertainty'] < limit_uncertainty, 
                                        df['uncertainty'], 
                                        0.5)
        return df

    @property
    def moments(self) -> pd.DataFrame:
        # Assuming gm_moments is accessible in the scope
        moms_df = gm_moments(self.means, self.sd, self.p)
        
        if self._x is not None and self._w is not None:
            nobs = len(self._x[self._w != 0])
        else:
            nobs = np.nan
            
        moms_df['nobs'] = nobs
        self._moments = moms_df.iloc[0].to_dict() # Update internal state
        return moms_df

    @property
    def summary(self) -> pd.DataFrame:
        # Generate parameter names dynamically since we dropped them from __init__
        components = len(self.means)
        terms = [f"mu{i+1}" for i in range(components)] + \
                [f"sd{i+1}" for i in range(components)] + \
                [f"p{i+1}" for i in range(components)]
                
        estimates = np.concatenate([self.means, self.sd, self.p])
        stderrs = np.concatenate([self.std_errors['means'], self.std_errors['sd'], self.std_errors['p']])
        
        df = pd.DataFrame({
            'term': terms,
            'estimate': estimates,
            'std.error': stderrs
        })
        
        df['statistic'] = df['estimate'] / df['std.error']
        
        # Corrected p-value calculation
        df['p_value'] = 2.0 * (1.0 - norm.cdf(np.abs(df['statistic'])))
        
        return df