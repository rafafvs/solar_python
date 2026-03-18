import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from patsy import build_design_matrices
import warnings

class FormulaString(str):
    """
    Helper class to replicate R's ability to attach attributes (like 'coef_names') 
    to a formula string.
    """
    def __new__(cls, value, coef_names=None):
        obj = str.__new__(cls, value)
        obj.coef_names = coef_names if coef_names is not None else []
        return obj

def seasonalModel_formula(formula, order=1, period=365, sin=True, cos=True, t_idx="n"):
    """
    Create a Fourier formula
    """
    if order == 0:
        return formula
    
    # Extract existing coefficient names if the attribute exists
    coefs_names_attr = getattr(formula, "coef_names", [])
    coefs_names = []
    base_formula = str(formula)
    
    # Add sin terms
    if sin:
        new_coef_name = f"sin_{order}_{period}"
        if new_coef_name not in coefs_names_attr:
            # Translated `base::pi` to `np.pi` for Python compatibility (e.g., Patsy/Statsmodels)
            base_formula += f" + I(np.sin(2 * np.pi / {period} * {t_idx} * {order}))"
            coefs_names.append(new_coef_name)
            
    # Add cos terms
    if cos:
        new_coef_name = f"cos_{order}_{period}"
        if new_coef_name not in coefs_names_attr:
            base_formula += f" + I(np.cos(2 * np.pi / {period} * {t_idx} * {order}))"
            coefs_names.append(new_coef_name)

    # Return new string with updated attributes
    return FormulaString(base_formula, coef_names=coefs_names_attr + coefs_names)


def seasonalModel_formula_dt(formula, order=1, period=365, sin=True, cos=True, t_idx="n"):
    """
    Create a Fourier formula for the differential
    """
    if order == 0:
        return formula
    
    coefs_names_attr = getattr(formula, "coef_names", [])
    coefs_names = []
    base_formula = str(formula)
    
    # Add sin terms
    if sin:
        new_coef_name = f"sin_{order}_{period}"
        if new_coef_name not in coefs_names_attr:
            base_formula += f" + I({order} * 2 * np.pi / {period} * np.cos(2 * np.pi / {period} * {t_idx} * {order}))"
            coefs_names.append(new_coef_name)
            
    # Add cos terms
    if cos:
        new_coef_name = f"cos_{order}_{period}"
        if new_coef_name not in coefs_names_attr:
            # NOTE: Preserved exact original logic, which computes the derivative of cos without the negative sign
            base_formula += f" + I({order} * 2 * np.pi / {period} * np.sin(2 * np.pi / {period} * {t_idx} * {order}))"
            coefs_names.append(new_coef_name)

    return FormulaString(base_formula, coef_names=coefs_names_attr + coefs_names)


def seasonalModel_params_to_zeta(b):
    """
    From constraint to unconstrained parameters
    """
    # Accommodate pandas Series to preserve names like R's named vectors
    b_vals = b.values if isinstance(b, pd.Series) else b
    
    # Reparametrized coefficients (0-based indexing)
    b0_star = np.log(b_vals[0])
    b1_star = np.arctanh(np.sqrt(b_vals[1]**2 + b_vals[2]**2) / b_vals[0])
    b2_star = np.arctan2(b_vals[1], b_vals[2]) # np.arctan2(y, x) matches R's atan2(y, x)
    
    b_star_vals = [b0_star, b1_star, b2_star]
    
    if isinstance(b, pd.Series):
        names = [f"{name}_star" for name in b.index]
        return pd.Series(b_star_vals, index=names)
    else:
        return np.array(b_star_vals)


def seasonalModel_params_to_phi(b_star):
    """
    From unconstrained to constraint parameters
    """
    b_star_vals = b_star.values if isinstance(b_star, pd.Series) else b_star
    
    # Reparametrized coefficients
    b0 = np.exp(b_star_vals[0])
    b1 = b0 * np.tanh(b_star_vals[1]) * np.sin(b_star_vals[2])
    b2 = b0 * np.tanh(b_star_vals[1]) * np.cos(b_star_vals[2])
    
    b_vals = [b0, b1, b2]
    
    if isinstance(b_star, pd.Series):
        # stringr::str_remove_all equivalent
        names = [str(name).replace("_star", "") for name in b_star.index]
        return pd.Series(b_vals, index=names)
    else:
        return np.array(b_vals)


def seasonalModel_params_to_zeta_jacobian(b_star):
    """
    Jacobian from unconstrained to constraint parameters
    """
    b_star_vals = b_star.values if isinstance(b_star, pd.Series) else b_star
    
    # Extract the coefficients
    b0_star = b_star_vals[0]
    b1_star = b_star_vals[1]
    b2_star = b_star_vals[2]
    
    # Initialize 3x3 matrix
    J = np.zeros((3, 3))
    
    # Derivatives of b0 wrt b0, b1, b2
    J[0, 0] = np.exp(b0_star)                                  # d_b0_d_b0 = b0
    J[0, 1] = 0                                                # d_b0_d_b1
    J[0, 2] = 0                                                # d_b0_d_b2
    
    # Derivatives of b1 wrt b0, b1, b2
    J[1, 0] = J[0, 0] * np.tanh(b1_star) * np.sin(b2_star)     # d_b1_d_b0 = b1
    J[1, 1] = J[0, 0] * (1 - np.tanh(b1_star)**2) * np.sin(b2_star) # d_b1_d_b1
    J[1, 2] = J[0, 0] * np.tanh(b1_star) * np.cos(b2_star)     # d_b1_d_b2 = b2
    
    # Derivatives of b2 wrt b0, b1, b2
    J[2, 0] = J[1, 2]                                          # d_b2_d_b0 = b2
    J[2, 1] = J[0, 0] * (1 - np.tanh(b1_star)**2) * np.cos(b2_star) # d_b2_d_b1
    J[2, 2] = -J[1, 0]                                         # d_b2_d_b2 = -b1
    
    colnames = ["d_b0_star", "d_b1_star", "d_b2_star"]
    rownames = ["b0", "b1", "b2"]
    
    return pd.DataFrame(J, index=rownames, columns=colnames)

class SeasonalModel:
    """
    Seasonal Model class
    
    The `SeasonalModel` class implements a seasonal regression model as a linear combination of sine and cosine functions.
    This model is designed to capture periodic effects in time series data, particularly for applications involving seasonal trends.
    
    Version 1.0.1
    """

    def __init__(self, order=1, period=365):
        """
        Initialize a `SeasonalModel` object.
        :param order: int or list, number of combinations of sines and cosines.
        :param period: int or list, seasonality period. The default is 365.
        """
        # Ensure order and period are lists for the iteration logic in fit
        self.__period = [period] if np.isscalar(period) else list(period)
        self.__order = [order] if np.isscalar(order) else list(order)
        
        self.extra_params = {}
        
        # Private fields equivalent
        self.__version = "1.0.1"
        self.__model = None
        self.__dmodel = None
        self._mformula = None
        self._dformula = None
        self.__std_errors = pd.Series(dtype=float)
        self._external_regressors = []
        
        # Mutable parameters for prediction (since statsmodels Results are immutable)
        self.__model_params = pd.Series(dtype=float)
        self.__dmodel_params = pd.Series(dtype=float)
        self.__coef_names = []
        self.__dcoef_names = []

    def fit(self, formula, data, **kwargs):
        """
        Fit a seasonal model as a linear combination of sine and cosine functions and
        eventual external regressors specified in the formula.
        """
        if not hasattr(formula, "coef_names"):
            formula = FormulaString(str(formula), coef_names=[])
        
        base_formula = formula
        base_formula_dt = formula

        ## Handle string inputs gracefully by wrapping them in the custom FormulaString (from previous script)
        #if not hasattr(base_formula, 'coef_names'):
        #    from types import SimpleNamespace
        #    base_formula = SimpleNamespace(coef_names=[], value=str(formula))
        #    # Mocking the string conversion for simplicity if FormulaString isn't imported
        #    setattr(base_formula, '__str__', lambda self: self.value)
        #    
        #if not hasattr(base_formula_dt, 'coef_names'):
        #    from types import SimpleNamespace
        #    base_formula_dt = SimpleNamespace(coef_names=[], value=str(formula))
        #    setattr(base_formula_dt, '__str__', lambda self: self.value)

        # Loop through periods and orders to build formulas
        if len(self.__order) == 1:
            for p in self.__period:
                for o in self.__order:
                    for i in range(1, o + 1):
                        base_formula = seasonalModel_formula(base_formula, order=i, period=p)
                        base_formula_dt = seasonalModel_formula_dt(base_formula_dt, order=i, period=p)
        elif len(self.__order) == len(self.__period):
            for j in range(len(self.__period)):
                p = self.__period[j]
                for i in range(1, self.__order[j] + 1):
                    base_formula = seasonalModel_formula(base_formula, order=i, period=p)
                    base_formula_dt = seasonalModel_formula_dt(base_formula_dt, order=i, period=p)

        self._mformula = str(base_formula)
        self._dformula = str(base_formula_dt)
        
        # Fit seasonal model
        self.__model = smf.ols(self._mformula, data=data, **kwargs).fit()
        self.__dmodel = smf.ols(self._dformula, data=data, **kwargs).fit()
        
        # Store coefficients locally to allow updating
        self.__model_params = self.__model.params.copy()
        dcoefs = self.__dmodel.params.copy()
        
        # Map differential coefficients (Sin to -Sin, Intercept to 0)
        sin_mask = dcoefs.index.str.contains("sin")
        dcoefs.loc[sin_mask] = -dcoefs.loc[sin_mask]
        
        int_mask = dcoefs.index.str.contains("Intercept")
        if int_mask.any():
            dcoefs.loc[int_mask] = 0
            
        self.__dmodel_params = dcoefs

        # Detect regressors from statsmodels exog_names
        exog_names = self.__model.model.exog_names
        
        idx_seasonal = [i for i, x in enumerate(exog_names) if "sin" in x or "cos" in x]
        idx_external = [i for i, x in enumerate(exog_names) if "sin" not in x and "cos" not in x and "Intercept" not in x]
        
        self._external_regressors = [exog_names[i] for i in idx_external]
        
        n_seasonal_reg = len(idx_seasonal)
        n_external_regressors = len(idx_external)
        
        # Map standard coefficient names
        coefs_names = exog_names.copy()
        seasonal_names = getattr(base_formula, "coef_names", [])
        
        for k, idx in enumerate(idx_seasonal):
            if k < len(seasonal_names):
                coefs_names[idx] = seasonal_names[k]
                
        if len(exog_names) - n_seasonal_reg - n_external_regressors == 1:
            int_idx = exog_names.index("Intercept") if "Intercept" in exog_names else 0
            coefs_names[int_idx] = "Intercept"
            
        self.__coef_names = coefs_names

        # Build equivalent "standard names" for the differential model too.
        d_exog_names = self.__dmodel.model.exog_names
        d_idx_seasonal = [i for i, x in enumerate(d_exog_names) if "np.sin" in x or "np.cos" in x]
        d_coefs_names = d_exog_names.copy()
        for k, idx in enumerate(d_idx_seasonal):
            if k < len(seasonal_names):
                d_coefs_names[idx] = seasonal_names[k]
        if "Intercept" in d_exog_names:
            d_coefs_names[d_exog_names.index("Intercept")] = "Intercept"
        self.__dcoef_names = d_coefs_names
        
        # Extract Standard Errors
        self.__std_errors = self.__model.bse.copy()
        self.__std_errors.index = self.__coef_names

    def predict(self, n=None, newdata=None):
        """
        Predict method for the class `seasonalModel`.
        """
        if newdata is None:
            if n is None:
                # Predict on training data
                X_exog = self.__model.model.exog
                return pd.Series(np.dot(X_exog, self.__model_params.values))
            else:
                n_val = n if isinstance(n, (list, tuple, np.ndarray, pd.Series)) else [n]
                newdata = pd.DataFrame({'n': n_val})
                
        # Use patsy to recreate the exact design matrix format from newdata
        X_df = build_design_matrices([self.__model.model.data.design_info], newdata, return_type='dataframe')[0]
        
        # Dot product with updated parameters
        return X_df.dot(self.__model_params)

    def differential(self, n=None, newdata=None):
        """
        Compute the differential of the sinusoidal function.
        """
        if newdata is None:
            if n is None:
                # Predict on training data
                X_exog = self.__dmodel.model.exog
                return pd.Series(np.dot(X_exog, self.__dmodel_params.values))
            else:
                n_val = n if isinstance(n, (list, tuple, np.ndarray, pd.Series)) else [n]
                newdata = pd.DataFrame({'n': n_val})
                
        # Build design matrix for differential model
        X_df = build_design_matrices([self.__dmodel.model.data.design_info], newdata, return_type='dataframe')[0]
        return X_df.dot(self.__dmodel_params)

    def update(self, coefficients):
        """
        Update the model's parameters.
        :param coefficients: dict or pd.Series, new parameters.
        """
        new_coefs = self.coefficients.copy()
        names_old = new_coefs.index.tolist()
        
        coef_series = coefficients if isinstance(coefficients, pd.Series) else pd.Series(coefficients)
        names_new = coef_series.index.tolist()
        
        if len(names_new) != len(names_old):
            warnings.warn("In seasonalModel$update(): The length of new `coefficients` do not match the length of the old coefficients.")
            
        # Update values safely matching the R logic
        for i, name in enumerate(names_new):
            if name in names_old:
                new_coefs[name] = coef_series.iloc[i]
                self.__std_errors[name] = np.nan
                
        # Remap standard names back to statsmodels exog_names for matrix operations
        mapped_params = pd.Series(new_coefs.values, index=self.__model.model.exog_names)
        self.__model_params = mapped_params
        
        # Update coefficients of the differential, aligned to differential design matrix columns
        # Map by the stored "standard names" for the differential model.
        d_vals = new_coefs.reindex(self.__dcoef_names).values
        dcoefs = pd.Series(d_vals, index=self.__dmodel.model.exog_names)

        # Apply sign/intercept adjustments matching fit()
        sin_mask = dcoefs.index.str.contains("np.sin", regex=False)
        dcoefs.loc[sin_mask] = -dcoefs.loc[sin_mask]

        if "Intercept" in dcoefs.index:
            dcoefs.loc["Intercept"] = 0

        self.__dmodel_params = dcoefs

    def update_std_errors(self, std_errors):
        """
        Update the parameter's std. errors.
        """
        new_std_errors = self.std_errors.copy()
        names_old = new_std_errors.index.tolist()
        
        err_series = std_errors if isinstance(std_errors, pd.Series) else pd.Series(std_errors)
        names_new = err_series.index.tolist()
        
        if len(names_new) != len(names_old):
            warnings.warn("In seasonalModel$update_std.errors(): The length of new `std.errors` do not match the length of the old std. errors!")
            
        for i, name in enumerate(names_new):
            if name in names_old:
                new_std_errors[name] = err_series.iloc[i]
                
        self.__std_errors = new_std_errors

    def __str__(self):
        """Print method for the class."""
        msg = "----------------------- seasonalModel ----------------------- \n"
        msg += f" - Order: {self.order}\n - Period: {self.period}\n"
        
        if not self._external_regressors:
            msg += "- External regressors: 0 \n"
        else:
            n_ext = len(self._external_regressors)
            ext_names = ", ".join(self._external_regressors)
            msg += f"- External regressors: {n_ext} ({ext_names})\n"
            
        msg += f"- Version: {self.__version}\n"
        msg += "--------------------------------------------------------------\n"
        if self.__model is not None:
            msg += str(self.__model.summary().tables[1]) # Only print coef table to mimic R print
        return msg
        
    def __repr__(self):
        return self.__str__()

    # Active bindings / Properties
    
    @property
    def coefficients(self):
        """Named vector, fitted coefficients (using standard names)."""
        coefs = self.__model_params.copy()
        coefs.index = self.__coef_names
        return coefs

    @property
    def model(self):
        """A slot with the fitted object."""
        return self.__model

    @property
    def period(self):
        """Integer scalar, periodicity of the seasonality."""
        return self.__period[0] if len(self.__period) == 1 else self.__period

    @property
    def order(self):
        """Integer scalar, number of combinations of sines and cosines."""
        return self.__order[0] if len(self.__order) == 1 else self.__order

    @property
    def omega(self):
        """Integer, periodicity."""
        p = np.array(self.__period)
        return 2 * np.pi / p

    @property
    def std_errors(self): # snake_case to match Python conventions
        """Named vector, with the parameters' std. errors."""
        return self.__std_errors

    @property
    def tidy(self):
        """A dataframe with estimated parameters and std. errors."""
        return pd.DataFrame({
            'term': self.coefficients.index,
            'estimate': self.coefficients.values,
            'std.error': self.__std_errors.values
        })