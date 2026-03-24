import numpy as np
import pandas as pd
from scipy.optimize import minimize
import copy
import warnings
from .seasonalModel import SeasonalModel
from .seasonalSolarFunctions import SeasonalSolarFunctions

def control_seasonalClearsky(order=1, order_H0=1, period=365, include_intercept=True, 
                             include_trend=False, delta0=1.4, lower=0, upper=3, 
                             by=0.001, ntol=0, quiet=False):
    """
    Control parameters for a `seasonalClearsky` object
    Returns a dictionary of control parameters.
    """
    # Equivalent to R's structure(..., class = c("control", "list"))
    return {
        "order": order,
        "order_H0": order_H0,
        "period": period,
        "include.intercept": include_intercept,  # Preserved dot notation in key for data mapping
        "include.trend": include_trend,
        "delta0": delta0,
        "lower": lower,
        "upper": upper,
        "by": by,
        "ntol": ntol,
        "quiet": quiet
    }

def clearsky_delta_optimizer(x, Ct, lower=0, upper=3, by=0.01, ntol=0):
    """
    Optimizer for Solar Clear sky.
    Find the best parameter delta for fitting clear sky radiation.
    """
    x = np.asarray(x)
    Ct = np.asarray(Ct)
    
    # Grid of points (R's seq includes upper bound)
    grid = np.arange(lower, upper + by/2, by)
    
    # Loss: Equivalent to purrr::map_dbl
    loss = np.array([np.sum((delta * Ct) - x < 0) for delta in grid])
    
    # opt dataframe and filtering
    opt = pd.DataFrame({'delta': grid, 'loss': loss})
    
    # Return minimum delta satisfying the constraint (R's [1] translates to .iloc[0])
    valid_opt = opt[opt['loss'] <= ntol]
    if not valid_opt.empty:
        delta = valid_opt.iloc[0]['delta']
    else:
        delta = np.nan # Fallback if no delta satisfies the condition
        
    return delta

def clearsky_optimizer(seasonal_model_Ct, newdata, ntol=0):
    """
    Optimizer for clear sky model with restricted least squares (RLS).
    """
    # Clone the model. Using deepcopy or R6 equivalent .clone() 
    if hasattr(seasonal_model_Ct, 'clone'):
        sm = seasonal_model_Ct.clone(True)
    else:
        sm = copy.deepcopy(seasonal_model_Ct)
        
    def loss_function(params, sm_obj, data):
        # Update the parameters
        # If params is a numpy array but the class expects a pd.Series, we wrap it
        if hasattr(sm_obj, 'coefficients') and isinstance(sm_obj.coefficients, pd.Series):
            params_to_pass = pd.Series(params, index=sm_obj.coefficients.index)
        else:
            params_to_pass = params
            
        sm_obj.update(params_to_pass)
        
        # Prediction
        pred = sm_obj.predict(newdata=data)
        
        # Violations: ifelse(data$GHI - pred > 0, 1, 0)
        violations = np.where(data['GHI'] - pred > 0, 1, 0)
        
        # Check number of violations lower than ntol
        mse = np.sum((data['GHI'] - pred)**2) + 1000000 * (np.sum(violations) - ntol)
        return mse

    # Initial parameters
    init_params = sm.coefficients
    init_vals = init_params.values if isinstance(init_params, pd.Series) else init_params

    # Optimal parameters (R's optim defaults to Nelder-Mead for multi-variable without bounds)
    opt = minimize(loss_function, init_vals, args=(sm, newdata), method='Nelder-Mead')
    
    # Update the parameters with the optimized array
    if isinstance(init_params, pd.Series):
        final_params = pd.Series(opt.x, index=init_params.index)
    else:
        final_params = opt.x
        
    sm.update(final_params)
    
    return sm

def clearsky_outliers(x, Ct, date=None, threshold=0.0001, quiet=False):
    """
    Impute clear sky outliers.
    Detect and impute outliers with respect to a maximum level of radiation (Ct)
    """
    # Initialize a dataset
    data = pd.DataFrame({'Ct': Ct, 'x': x})
    
    # Eventually add a date for non-stationary data
    if date is not None:
        data['date'] = pd.to_datetime(date)
        data['Month'] = data['date'].dt.month
        data['Day'] = data['date'].dt.day

    # Number of observations
    nobs = len(data)
    
    # Detect problems and violations
    outliers_na = data.index[data['x'].isna()].tolist()
    outliers_lo = data.index[(data['x'] <= 0) & (~data['x'].isna())].tolist()
    outliers_hi = data.index[(data['x'] >= data['Ct']) & (~data['x'].isna())].tolist()
    
    # Complete outliers index (unique elements)
    idx_outliers = list(set(outliers_na + outliers_lo + outliers_hi))
    
    # Check presence of outliers and impute them
    if not idx_outliers:
        if not quiet:
            print("No outliers!")
        data_clean = data.copy()
    else:
        # Verbose message
        if not quiet:
            percentage = (len(idx_outliers) / nobs) * 100
            print(f"Outliers: {len(idx_outliers)} ({percentage:.2f} %)")
            
        # Dataset without outliers
        data_no_outliers = data.drop(index=idx_outliers)
        
        # Imputed dataset
        data_clean = data.copy()
        
        # Impute outliers
        for i in idx_outliers:
            df_n = data.loc[i]
            
            if date is not None:
                # Data for the same day and month (without outliers)
                mask = (data_no_outliers['Month'] == df_n['Month']) & \
                       (data_no_outliers['Day'] == df_n['Day']) & \
                       (data_no_outliers['date'] != df_n['date'])
                df_day = data_no_outliers[mask]
            else:
                df_day = data_no_outliers

            # Imputed data depending on outliers type
            if i in outliers_na:
                # Outlier is an NA
                data_clean.at[i, 'x'] = df_day['x'].mean()
            elif i in outliers_lo:
                # Outlier is under minimum value for the day
                data_clean.at[i, 'x'] = df_day['x'].min()
            elif i in outliers_hi:
                # Outlier is above maximum value for the day
                data_clean.at[i, 'x'] = data_clean.at[i, 'Ct'] * (1 - threshold)

    if date is None:
        out_date = None
    else:
        out_date = data.loc[idx_outliers, 'date'].tolist()

    # Calculate errors identically to R
    # mean(abs((data$x - data_clean$x)/data$x))*100
    mape = np.mean(np.abs((data['x'] - data_clean['x']) / data['x'])) * 100
    # sd(data$x - data_clean$x)
    mse = np.std(data['x'] - data_clean['x'], ddof=1) # ddof=1 for sample standard deviation like R's sd()

    # Structure output data
    return {
        "x": data_clean['x'].values,
        "original": data.loc[idx_outliers, 'x'].values,
        "imputed": data_clean.loc[idx_outliers, 'x'].values,
        "index": idx_outliers,  # 0-based indexing in Python
        "index_type": {"na": outliers_na, "lo": outliers_lo, "hi": outliers_hi},
        "date": out_date,
        "n": len(idx_outliers),
        "MAPE": mape,
        "MSE": mse,
        "threshold": threshold
    }

class SeasonalClearsky(SeasonalModel):
    """
    R6 implementation for a clear sky seasonal model
    Version 1.0.1
    """

    def __init__(self, control=None):
        """
        Initialize a `seasonalClearsky` object.
        :param control: dict, control parameters. See `control_seasonalClearsky` for more details.
        """
        if control is None:
            control = control_seasonalClearsky()
            
        # Call parent initializer to set up the inherited attributes
        super().__init__(order=control['order'], period=control['period'])
        
        self.lat = np.nan
        
        # Private fields
        self.__version = "1.0.1"
        self.__coefficients_orig = None
        self.__delta = np.nan
        self.__ssf = None
        self.__control = control
        
        # R code directly overwrote these inherited private properties, 
        # but in Python they are already set via super().__init__()
        # self._SeasonalModel__order = control['order']
        # self._SeasonalModel__period = control['period']

    @property
    def control(self):
        """Named list (dict in Python), control parameters."""
        return self.__control

    @property
    def ssf(self):
        """Solar Seasonal Functions"""
        return self.__ssf

    def fit(self, x, date, lat, clearsky=None, method="paper"):
        """
        Fit the seasonal model for clear sky radiation.
        :param x: Numeric vector, time series of actual solar radiation (GHI).
        :param date: Time index/dates.
        :param lat: Reference latitude.
        :param clearsky: CAMS clear sky data.
        :param method: 'repo' (2-step uniform scaling heuristic) or 'paper' (1-step Constrained Least Squares).
        """

        if method not in ["repo", "paper"]:
            raise ValueError("method must be either 'repo' or 'paper'")

        control = self.control
        include_intercept = control['include.intercept']
        include_trend = control['include.trend']
        
        if clearsky is None:
            raise ValueError("`clearsky` time series must be specified.")
            
        # Add the function to compute extraterrestrial radiation
        self.__ssf = SeasonalSolarFunctions("spencer")
        
        # Store reference latitude (handle scalar or array-like)
        self.lat = lat[0] if isinstance(lat, (list, tuple, np.ndarray, pd.Series)) else lat
        
        # Initialize the dataset
        date_series = pd.to_datetime(date)
        data = pd.DataFrame({'date': date_series})
        data['Year'] = data['date'].dt.year
        data['Month'] = data['date'].dt.month
        data['Day'] = data['date'].dt.day
        data['t'] = data['Year'] - data['Year'].max()
        data['n'] = data['date'].dt.dayofyear
        data['Rt'] = np.asarray(x)
        
        # Note: Hon requires alt. We pass alt=None as the R code implicitly omitted it.
        data['H0'] = self.ssf.Hon(data['n'], self.lat, alt=None)
        data['clearsky'] = np.asarray(clearsky)
        
        # ========================================================================
        # 1. Daily maximum clearsky
        # ========================================================================
        base_formula = "clearsky ~ H0"
        
        if control['order_H0'] > 1:
            for i in range(2, control['order_H0'] + 1):
                col_name = f"H0_{i}"
                data[col_name] = data['H0'] ** i
                base_formula += f" + {col_name}"
                
        if include_trend:
            base_formula += " + t"
            
        if not include_intercept:
            base_formula += " - 1"
            
        formula_to_fit = getattr(base_formula, 'value', str(base_formula))
        
        # Fit the coefficients of the clear sky max model
        # We always run the unconstrained OLS first to build the design matrix and get a starting point
        super().fit(formula=formula_to_fit, data=data)
        
        # ========================================================================
        # 2. Optimization Branching
        # ========================================================================
        if method == "repo":
            # --- The Package's 2-Step Heuristic ---
            data['Ct_hat'] = self.predict(newdata=data)
            
            delta_val = clearsky_delta_optimizer(
                data['Rt'], data['Ct_hat'] * control['delta0'], 
                control['lower'], control['upper'], control['by'], control['ntol']
            )
            
            # Store delta parameter
            self.__delta = delta_val * control['delta0']
            
            # Scale OLS coefficients and std errors uniformly
            current_coefs = self.coefficients.copy() * self.__delta
            current_std_errs = self.std_errors.copy() * self.__delta

        elif method == "paper":
            # --- The Paper's 1-Step Constrained Least Squares (CLS) ---
            from scipy.optimize import minimize, LinearConstraint
            
            # Extract Design Matrix (X) and Target (y) from the statsmodels object
            X = self.model.model.exog
            y_target = data['clearsky'].values
            ghi = data['Rt'].values
            
            def objective_function(beta):
                return np.sum((X.dot(beta) - y_target)**2)

            def objective_jacobian(beta):
                return 2 * X.T.dot(X.dot(beta) - y_target)

            # Strict constraint: C_t >= R_t  =>  X @ beta >= GHI
            linear_constraint = LinearConstraint(X, lb=ghi, ub=np.inf)
            ols_start = self.coefficients.values
            
            if not control.get('quiet', False):
                print("Running Constrained Least Squares (Paper method)...")
                
            cls_result = minimize(
                objective_function, 
                x0=ols_start, 
                method='trust-constr', 
                jac=objective_jacobian,
                constraints=[linear_constraint],
                options={'maxiter': 10000, 'verbose': 0}
            )
            
            if not cls_result.success:
                warnings.warn(f"Constrained optimization failed or reached maxiter: {cls_result.message}")
            
            self.__delta = 1.0 # No scalar multiplier used in this method
            current_coefs = pd.Series(cls_result.x, index=self.coefficients.index)
            # The paper keeps the OLS standard errors (dependent purely on X'X) 
            current_std_errs = self.std_errors.copy() 
        
        # ========================================================================
        # 3. Naming Convention & Property Updating (Shared Logic)
        # ========================================================================

        # Standard names for coefficients
        coefs_names = []
        # Get original names from SeasonalModel (using property rather than mangled private attribute directly)
        orig_names = list(self.coefficients.index)
        
        if include_intercept:
            coefs_names.append("delta_0")
            orig_names.pop(0)
            for i in range(1, control['order_H0'] + 1):
                coefs_names.append(f"delta_extra{i}")
                orig_names.pop(0)
        else:
            for i in range(1, control['order_H0'] + 1):
                coefs_names.append(f"delta_extra{i}")
                orig_names.pop(0)
                
        if include_trend:
            coefs_names.append("t")
            orig_names.pop(0)
            
        # Ensure self.order evaluates correctly if it's a list or scalar
        order_val = self.order[0] if isinstance(self.order, list) else self.order
        if order_val > 0:
            for name in orig_names:
                coefs_names.append(f"delta_{name}")

        # Store original coefficients
        self.__coefficients_orig = self.coefficients.copy()
        
        # Update coefficients values
        super().update(current_coefs)
        
        # Directly update the coefficient names to mimic R's env manipulation
        self._SeasonalModel__coef_names = coefs_names
        
        # Update std errors values and names
        super().update_std_errors(current_std_errs)
        
        # Re-index std_errors to reflect the new names
        new_std_errs = self.std_errors
        new_std_errs.index = coefs_names
        self._SeasonalModel__std_errors = new_std_errs

    def predict(self, n=None, newdata=None):
        """
        Predict method for `seasonalClearsky` object.
        """
        if newdata is None:
            if n is None:
                # Assuming super() predict maps correctly without args to fitted values
                return super().predict()
            else:
                n_arr = np.atleast_1d(n)
                H0 = self.ssf.Hon(n_arr, self.lat, alt=None)
                newdata_df = pd.DataFrame({'n': n_arr, 'H0': H0})
                
                if self.control['order_H0'] > 1:
                    for i in range(2, self.control['order_H0'] + 1):
                        newdata_df[f"H0_{i}"] = newdata_df['H0'] ** i
                        
                return super().predict(newdata=newdata_df)
        else:
            newdata = newdata.copy()
            newdata['H0'] = self.ssf.Hon(newdata['n'], self.lat, alt=None)
            
            if self.control['order_H0'] > 1:
                for i in range(2, self.control['order_H0'] + 1):
                    newdata[f"H0_{i}"] = newdata['H0'] ** i
                    
            return super().predict(newdata=newdata)

    def differential(self, n=None, newdata=None):
        """
        Differential method for `seasonalClearsky` object.
        """
        if newdata is None:
            if n is None:
                return super().differential()
            else:
                n_arr = np.atleast_1d(n)
                H0 = self.ssf.Hon(n_arr, self.lat, alt=None, deriv=True)
                newdata_df = pd.DataFrame({'n': n_arr, 'H0': H0})
                
                if self.control['order_H0'] > 1:
                    for i in range(2, self.control['order_H0'] + 1):
                        # NOTE: Translating exactly as written in R.
                        newdata_df[f"H0_{i}"] = i * (newdata_df['H0'] ** (i - 1))
                        
                return super().differential(newdata=newdata_df)
        else:
            newdata = newdata.copy()
            newdata['H0'] = self.ssf.Hon(newdata['n'], self.lat, alt=None, deriv=True)
            
            if self.control['order_H0'] > 1:
                for i in range(2, self.control['order_H0'] + 1):
                    # NOTE: Translating exactly as written in R.
                    newdata[f"H0_{i}"] = i * (newdata['H0'] ** (i - 1))
                    
            return super().differential(newdata=newdata)

    def __str__(self):
        """Print method for `seasonalClearsky` object."""
        msg = "----------------------- seasonalClearsky ----------------------- \n"
        msg += f" - Order: {self.order}\n - Period: {self.period}\n"
        msg += "- External regressors: 1 (H0) \n"
        msg += f"- Version: {self.__version}\n"
        msg += "--------------------------------------------------------------\n"
        if self.model is not None:
            msg += str(self.model.summary().tables[1])
        return msg
        
    def __repr__(self):
        return self.__str__()

