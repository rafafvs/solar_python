import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
import re
import warnings

# --- External Helper Functions (Translating the implied R external functions) ---

def seasonal_model_formula(formula, order, period):
    """
    Literal translation of seasonalModel_formula.
    Appends sin and cos terms to the formula string.
    """
    # R formula concatenation: formula + sin(...) + cos(...)
    term = f" + np.sin(2 * np.pi * {order} * n / {period}) + np.cos(2 * np.pi * {order} * n / {period})"
    return formula + term

def seasonal_model_formula_dt(formula, order, period):
    """
    Literal translation of seasonalModel_formula_dt.
    Assumed to be identical structure for the formula definition step.
    """
    return seasonal_model_formula(formula, order, period)

# --- Main Class Translation ---

class SeasonalModel:
    """
    Seasonal Model
    
    Version 1.0.1
    """

    def __init__(self, order=1, period=365):
        """
        Initialize a `seasonalModel` object.
        """
        self.extra_params = {}
        
        # Private attributes (mapping private$..var)
        self._period = period
        self._order = order
        self._model = None
        self._dmodel = None
        self._mformula = None
        self._dformula = None
        self._std_errors = None
        self._external_regressors = None
        self._version = "1.0.1"

    def fit(self, formula, data, **kwargs):
        """
        Fit a seasonal model.
        """
        # Formula with standard names
        base_formula = formula
        base_formula_dt = formula
        
        # Handle scalar vs vector order/period logic from R
        # In Python, we ensure they are lists to iterate
        orders = [self._order] if np.isscalar(self._order) else self._order
        periods = [self._period] if np.isscalar(self._period) else self._period
        
        if len(orders) == 1:
            for period in periods:
                for order in orders:
                    for i in range(1, order + 1):
                        base_formula = seasonal_model_formula(base_formula, order=i, period=period)
                        base_formula_dt = seasonal_model_formula_dt(base_formula_dt, order=i, period=period)
                        
        elif len(orders) == len(periods):
            for j in range(len(periods)):
                period = periods[j]
                for i in range(1, orders[j] + 1):
                    base_formula = seasonal_model_formula(base_formula, order=i, period=period)
                    base_formula_dt = seasonal_model_formula_dt(base_formula_dt, order=i, period=period)
        
        # Store the main formula
        self._mformula = base_formula
        # Note: Python formulas don't typically carry 'attr' like R. 
        # We assume the variable names are handled by statsmodels internals.

        # Fit seasonal model
        # R: lm(private$mformula, data = data)
        self._model = smf.ols(formula=self._mformula, data=data).fit(**kwargs)

        # Fit the differential w.r.t. t
        self._dformula = base_formula_dt
        # R: lm(private$dformula, data = data)
        # Note: We fit the structure here, coefficients are overwritten later.
        self._dmodel = smf.ols(formula=self._dformula, data=data).fit(**kwargs)

        # Update coefficients
        # R: dcoefs <- private$..model$coefficients
        dcoefs = self._model.params.copy()
        
        # Align names (In Python statsmodels, names are indices)
        # We assume the structure is identical so direct assignment works, 
        # but we use the dmodel's index to be safe.
        dcoefs.index = self._dmodel.params.index
        
        # Flip signs for sine terms
        # R: dcoefs[stringr::str_detect(names(dcoefs), "sin")] <- -dcoefs[...]
        for name in dcoefs.index:
            if "sin" in name:
                dcoefs[name] = -dcoefs[name]
            # Zero out Intercept
            # R: dcoefs[stringr::str_detect(names(dcoefs), "(Intercept)")] <- 0
            if "Intercept" in name:
                dcoefs[name] = 0.0
        
        self._dmodel.params = dcoefs

        # Detect number of regressors
        n_regressors = len(self._model.params)
        
        # Extract regressors from the formula excluding target variable
        # R: formula.tools::get.vars(formula(private$..model))[-c(1)]
        # Python equivalent: extracting exog names excluding Intercept
        regressors = list(self._model.params.index)
        if "Intercept" in regressors:
            regressors.remove("Intercept")
            
        # Seasonal regressors
        # R: which(stringr::str_detect(regressors, "sin|cos"))
        seasonal_regressors = [r for r in regressors if "sin" in r or "cos" in r]
        n_seasonal_reg = len(seasonal_regressors)
        
        # External regressors
        # R: which(!stringr::str_detect(regressors, "sin|cos"))
        self._external_regressors = [r for r in regressors if "sin" not in r and "cos" not in r]
        n_external_regressors = len(self._external_regressors)
        
        # Store standard names (Logic from R script regarding attributes skipped as it's R-specific,
        # but we preserve the logic of identifying them)
        coefs_names = regressors # Simplified placeholder for the complex R naming logic
        
        # Update std.errors
        # R: broom::tidy(private$..model)$std.error
        self._std_errors = self._model.bse.copy()

    def predict(self, n=None, newdata=None):
        """
        Predict method.
        """
        if newdata is None:
            if n is None:
                return self._model.predict()
            else:
                return self._model.predict(exog=pd.DataFrame({'n': n}))
        else:
            return self._model.predict(exog=newdata)

    def differential(self, n=None, newdata=None):
        """
        Compute the differential.
        """
        if newdata is None:
            if n is None:
                return self._dmodel.predict()
            else:
                return self._dmodel.predict(exog=pd.DataFrame({'n': n}))
        else:
            return self._dmodel.predict(exog=newdata)

    def update(self, coefficients):
        """
        Update the model's parameters.
        """
        # Extract old coefficients
        new_coefs = self.coefficients.copy()
        
        # Extract names
        names_old = new_coefs.index
        
        # Handle Dictionary or Series input for coefficients
        if isinstance(coefficients, dict):
            names_new = list(coefficients.keys())
            values_new = list(coefficients.values())
        elif isinstance(coefficients, pd.Series):
            names_new = coefficients.index
            values_new = coefficients.values
        else:
            # Fallback for list/array (assuming order matches)
            warnings.warn("Coefficients provided as list/array. Assuming order matches.")
            names_new = names_old[:len(coefficients)]
            values_new = coefficients

        # Warning
        if len(names_new) != len(names_old):
            warnings.warn("In seasonalModel.update(): The length of new `coefficients` does not match the length of the old coefficients.")
            
        # Update only if they are present
        for i in range(len(names_new)):
            name = names_new[i]
            val = values_new[i]
            if name in names_old:
                new_coefs[name] = val
                # R: private$..std.errors[names_new[i]] <- NA_integer_
                if self._std_errors is not None and name in self._std_errors:
                    self._std_errors[name] = np.nan

        # Update the parameters inside the lm object
        self._model.params = new_coefs

        # Update coefficients of the differential
        dcoefs = new_coefs.copy()
        # R: names(dcoefs) <- names(private$..dmodel$coefficients) (implied alignment)
        
        for name in dcoefs.index:
            if "sin" in name:
                dcoefs[name] = -dcoefs[name]
            if "Intercept" in name:
                dcoefs[name] = 0.0
                
        self._dmodel.params = dcoefs

    def update_std_errors(self, std_errors):
        """
        Update the parameter's std. errors.
        """
        new_std_errors = self.std_errors.copy()
        
        if isinstance(std_errors, dict):
            names_new = list(std_errors.keys())
            values_new = list(std_errors.values())
        elif isinstance(std_errors, pd.Series):
            names_new = std_errors.index
            values_new = std_errors.values
        else:
            names_new = new_std_errors.index[:len(std_errors)]
            values_new = std_errors

        names_old = new_std_errors.index

        if len(names_new) != len(names_old):
            warnings.warn("In seasonalModel.update_std_errors(): The length of new `std.errors` do not match the length of the old std. errors!")

        for i in range(len(names_new)):
            name = names_new[i]
            val = values_new[i]
            if name in names_old:
                new_std_errors[name] = val
        
        self._std_errors = new_std_errors

    def print(self):
        """
        Print method.
        """
        print("----------------------- seasonalModel -----------------------")
        msg_1 = f" - Order: {self.order}\n - Period: {self.period}\n"
        
        if not self._external_regressors:
            msg_2 = "- External regressors: 0 \n"
        else:
            n_ext = len(self._external_regressors)
            msg_2 = f"- External regressors: {n_ext} ({self._external_regressors})\n"
            
        msg_3 = f"- Version: {self._version}\n"
        print(msg_1 + msg_2 + msg_3, end="")
        print("--------------------------------------------------------------")
        if self._model is not None:
            print(self._model.params)
        else:
            print("Model not fitted")

    # --- Active Bindings (Properties) ---

    @property
    def coefficients(self):
        if self._model is None:
            return None
        return self._model.params

    @property
    def model(self):
        return self._model

    @property
    def period(self):
        return self._period

    @property
    def order(self):
        return self._order

    @property
    def omega(self):
        # R: 2 * base::pi / private$..period
        # Handle potential vector period
        p = np.array(self._period)
        return 2 * np.pi / p

    @property
    def std_errors(self):
        return self._std_errors

    @property
    def tidy(self):
        if self._model is None:
            return None
        return pd.DataFrame({
            'term': self.coefficients.index,
            'estimate': self.coefficients.values,
            'std.error': self._std_errors.values if self._std_errors is not None else [np.nan]*len(self.coefficients)
        })

    def __repr__(self):
        """Python's equivalent to default print behavior"""
        self.print()
        return ""