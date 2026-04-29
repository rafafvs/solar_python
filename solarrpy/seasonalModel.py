import numpy as np
import pandas as pd
import statsmodels.api as sm
import warnings

def add_seasonal_features(data, time_col="n", orders=[1], periods=[365], sin = True, cos = True):
    """Returns a new dataframe with Fourier columns added."""
    # Work on a copy to avoid mutating the original data
    X = data.copy()
    
    # Standardize inputs to lists so we can iterate cleanly
    periods = periods if isinstance(periods, list) else [periods]
    orders = orders if isinstance(orders, list) else [orders]
    
    # If only one order is provided, stretch it to match the number of periods
    orders_to_use = [orders[0]] * len(periods) if len(orders) == 1 else orders
    
    # Generate the features using native vectorized math
    for j, period in enumerate(periods):
        order = orders_to_use[j]
        for i in range(1, order + 1):
            if sin:
                X[f'sin_{i}_{period}'] = np.sin(2 * np.pi / period * X[time_col] * i)
            if cos:
                X[f'cos_{i}_{period}'] = np.cos(2 * np.pi / period * X[time_col] * i)

    # Return the new dataframe with the seasonal features
    return X

def add_seasonal_features_dt(data, time_col="n", orders=[1], periods=[365], sin = True, cos = True):
    """Returns a new dataframe with Fourier columns added."""
    # Work on a copy to avoid mutating the original data
    X = data.copy()

    # Standardize inputs to lists so we can iterate cleanly
    periods = periods if isinstance(periods, list) else [periods]
    orders = orders if isinstance(orders, list) else [orders]
    
    # If only one order is provided, stretch it to match the number of periods
    orders_to_use = [orders[0]] * len(periods) if len(orders) == 1 else orders
    
    # Generate the features using native vectorized math
    for j, period in enumerate(periods):
        order = orders_to_use[j]
        for i in range(1, order + 1):
            if sin:
                X[f'sin_{i}_{period}'] = 2 * np.pi * i / period * np.cos(2 * np.pi * i / period * X[time_col])
            if cos:
                X[f'cos_{i}_{period}'] = 2 * np.pi * i / period * np.sin(2 * np.pi * i / period * X[time_col])
    
    return X

def seasonalModel_params_to_zeta(b):
    """Constrained (b0, b1, b2) to unconstrained (b0_star, b1_star, b2_star)."""

    b0_star = np.log(b[0])
    b1_star = np.arctanh(np.sqrt(b[1]**2 + b[2]**2) / b[0])
    b2_star = np.arctan2(b[1], b[2])

    return np.array([b0_star, b1_star, b2_star])

def seasonalModel_params_to_phi(b_star):
    """Unconstrained (b0_star, b1_star, b2_star) to constrained (b0, b1, b2)."""
    b0 = np.exp(b_star[0])
    b1 = b0 * np.tanh(b_star[1]) * np.sin(b_star[2])
    b2 = b0 * np.tanh(b_star[1]) * np.cos(b_star[2])

    return np.array([b0, b1, b2])

def seasonalModel_params_to_zeta_jacobian(b_star):
    """Jacobian from unconstrained (z0, z1, z2) to constraint (b0, b1, b2)."""
    # Extract the coefficients
    b0_star = b_star[0]
    b1_star = b_star[1]
    b2_star = b_star[2]
    
    # Initialize
    J = np.zeros((3, 3))
    
    # Derivatives of b0 wrt b0, b1, b2
    J[1,1] = np.exp(b0_star)
    J[1,2] = 0
    J[1,3] = 0

    # Derivatives of b1 wrt b0, b1, b2
    J[2,1] = J[1,1] * np.tanh(b1_star) * np.sin(b2_star)
    J[2,2] = J[1,1] * (1 - np.tanh(b1_star)**2) * np.sin(b2_star)
    J[2,3] = J[1,1] * np.tanh(b1_star) * np.cos(b2_star)

    # Derivatives of b2 wrt b0, b1, b2
    J[3,1] = J[2,3]
    J[3,2] = J[1,1] * (1 - np.tanh(b1_star)**2) * np.cos(b2_star)
    J[3,3] = -J[2,1]

    J_df = pd.DataFrame(
        J, 
        columns=['d_b0_star', 'd_b1_star', 'd_b2_star'],
        index=['b0', 'b1', 'b2']
    )

    return J_df

class SeasonalModel:

    def __init__(self, orders = 1, periods = 365):
        """
        Initializes the SeasonalModel
        """

        self.version = "1.0.1"
        self.extra_params = {}
        
        # Standardize inputs to lists to support multiple seasonalities
        self._periods = periods if isinstance(periods, list) else [periods]
        self._orders = orders if isinstance(orders, list) else [orders]
        
        # State attributes (Initialized as None, populated during .fit() )
        self._model = None       # Will hold the statsmodels object
        self._dcoefs = None      # Will hold the derivative coefficients
        self._std_errors = None  # Will hold the standard errors
        self.external_regressors = [] # Will hold the external regressors

    def fit(self, data, target_col, time_col="n", external_regressors=[], include_intercept=True):
        """Fit the seasonal model."""

        # Store the external regressors
        self.external_regressors = external_regressors

        # 1. Initialize our matrices with the base regressors
        X = data[self.external_regressors + [time_col]].copy()
        y = data[target_col]

        # Add seasonal featuress
        X = add_seasonal_features(X, time_col=time_col, orders=self._orders, periods=self._periods)

        # Drop the raw time column because it's not a regressor itself
        X = X.drop(columns=[time_col])
        
        # 4. Fit the main model
        if include_intercept:
            X = sm.add_constant(X)

        self._model = sm.OLS(y, X).fit()

        self._std_errors = self._model.bse.copy() # Save the standard errors to use in update()
        
        # 5. Create the derivative coefficients (No second fit)
        self.dcoefs = self._model.params.copy()
        
        if 'const' in self.dcoefs:
            self.dcoefs['const'] = 0.0
        
        # Negates cos coefficients because of the derivative of cos(x) = -sin(x)
        # Different from R because R code used derivatives expressions 
        # and not the names of the coefficients
        cos_cols = [col for col in self.dcoefs.index if 'cos' in col]
        self.dcoefs[cos_cols] *= -1
        
        return self

    def predict(self, data, time_col="n"):
        """Predict method for the class `seasonalModel`."""

        # Handle missing data (Return training predictions)
        if data is None:
            return self._model.fittedvalues

        # Extract base columns needed for feature engineering
        self._training_data = data[self.external_regressors + [time_col]].copy()

        X = self._training_data.copy()
        # Add seasonal features
        X = add_seasonal_features(X, time_col=time_col, orders=self._orders, periods=self._periods)

        # Drop the raw time column because it's not a regressor itself
        X = X.drop(columns=[time_col])

        # Make sure the model includes the intercept
        if 'const' in self._model.params.index and 'const' not in X.columns:
            X = sm.add_constant(X, has_constant='add')
        
        # Select the columns that are in the model
        X = X[self._model.params.index]

        # Generate the predictions
        return self._model.predict(X)

    def differential(self, data=None, time_col="n"):
        """Compute the differential of the sinusoidal function."""

        # Handle missing data (Return training differential)
        if data is None:
            data = self._training_data  # saved during fit()
            
        # Extract base columns needed for feature engineering
        X = data[self.external_regressors + [time_col]].copy()
        
        # Add seasonal features
        X = add_seasonal_features_dt(X, time_col=time_col, orders=self._orders, periods=self._periods)
        
        # Drop the raw time column because it's not a regressor itself
        X = X.drop(columns=[time_col])

        # Make sure the model includes the intercept
        if 'const' in self._model.params.index and 'const' not in X.columns:
            X = sm.add_constant(X, has_constant='add')
        
        # Select the columns that are in the model
        X = X[self._model.params.index]
        
        # Generate the differential
        return X.dot(self.dcoefs)

    def update(self, coefficients):
        """Update model coefficients by name."""
        # 1. Warn if lengths differ
        if len(coefficients) != len(self._model.params):
            warnings.warn("In seasonalModel.update(): The length of new coefficients does not match the length of the old coefficients.")

        # 2. Update only coefficients that exist, invalidate their std errors
        for name, value in coefficients.items():
            if name in self._model.params.index:
                self._model.params[name] = value
                self._std_errors[name] = np.nan

        # 3. Recompute derivative coefficients from updated params
        self.dcoefs = self._model.params.copy()
        if 'const' in self.dcoefs:
            self.dcoefs['const'] = 0.0
        cos_cols = [col for col in self.dcoefs.index if 'cos' in col]
        self.dcoefs[cos_cols] *= -1
        return self
    
    def update_std_errors(self, std_errors):
        """Update standard errors of the model parameters."""
        
        # 1. Safety check for lengths
        if len(std_errors) != len(self._std_errors):
            warnings.warn(
                "In SeasonalModel.update_std_errors(): The length of `std_errors` "
                "does not match the length of the old standard errors."
            )

        # 2. Safe update (only overwrite keys that already exist in the model)
        for name, value in std_errors.items():
            if name in self._std_errors.index:
                self._std_errors[name] = value

        return self

    def __str__(self):
        """Print method for the class `seasonalModel`."""
        msg = f"----------------------- seasonalModel -----------------------\n"
        msg += f" - Order: {self._orders}\n"
        msg += f" - Period: {self._periods}\n"
        if len(self.external_regressors) == 0:
            msg += f" - External regressors: 0\n"
        else:
            msg += f" - External regressors: {len(self.external_regressors)} ({self.external_regressors})\n"
        msg += f" - Version: {self.version}\n"
        msg += f"--------------------------------------------------------------\n"
        msg += str(self._model.summary2(float_format="%.2f").tables[1])
        return msg

    def __repr__(self):
        return self.__str__()

    @property
    def coefficients(self):
        """Named Series of fitted coefficients."""
        # Because we saved this as a Pandas Series during fit(), 
        # it natively contains both the names (index) and the values.
        return self._model.params

    @property
    def model(self):
        """A slot with the fitted statsmodels object."""
        return self._model 

    @property
    def periods(self):
        """List of periodicities of the seasonality."""
        return self._periods

    @property
    def orders(self):
        """List of the number of combinations of sines and cosines."""
        return self._orders

    @property
    def omega(self):
        """List of calculated angular frequencies."""
        # Since we upgraded to support multiple periods, 
        # this uses a list comprehension to calculate omega for all of them.
        return [2 * np.pi / p for p in self._periods]

    @property
    def std_errors(self):
        """Series with the parameters' standard errors."""
        return self._std_errors

    @property
    def tidy(self):
        """A DataFrame with estimated parameters and standard errors."""
        return pd.DataFrame({
            "term": self._model.params.index,
            "estimate": self._model.params.values,
            "std.error": self._std_errors.values
        })