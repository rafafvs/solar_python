import numpy as np
import copy
import pandas as pd
from scipy.optimize import minimize

from .seasonalModel import SeasonalModel
from .seasonalSolarFunctions import SeasonalSolarFunctions
from .zzz import number_of_day

def control_seasonalClearsky(orders=1, order_H0=1, periods=365, include_intercept=True,
                             include_trend=False, delta0=1.4, lower=0, upper=3,
                             by=0.001, ntol=0, quiet=False):
    """Control parameters for a ``SeasonalClearsky`` object.

    Parameters
    ----------
    orders : int or list of int, default 1
        Fourier expansion order(s). Scalar matches the R API; a list enables
        multi-period extensions and is accepted transparently downstream
        (``SeasonalModel`` wraps non-list arguments to lists internally).
    periods : int or list of int, default 365
        Seasonal period(s) in days. Same scalar/list contract as ``orders``.
    """

    return {
        "orders": orders,
        "order_H0": order_H0,
        "periods": periods,
        "include_intercept": include_intercept,
        "include_trend": include_trend,
        "delta0": delta0,
        "lower": lower,
        "upper": upper,
        "by": by,
        "ntol": ntol,
        "quiet": quiet
    }

def clearsky_delta_optimizer(x, Ct, lower=0, upper=3, by=0.01, ntol=0):
    """Optimizer for the delta parameter in the clear sky model."""
    
    # Grid of candidate delta values
    grid = np.arange(lower, upper + by, by)
    
    # Loss function: count the number of violations where delta * Ct < GHI
    losses = np.array([np.sum(delta * Ct - x < 0) for delta in grid])
    
    # Return the first delta where violations <= ntol
    valid = np.where(losses <= ntol)[0]
    return grid[valid[0]] if len(valid) > 0 else None

def clearsky_optimizer(seasonal_model_Ct, data, ntol=0):
    # Deep copy to avoid mutating the original model
    sm = copy.deepcopy(seasonal_model_Ct)

    def loss_function(params, sm, data):
        # Update model with candidate parameters
        param_dict = dict(zip(sm.model.params.index, params))
        sm.update(param_dict)

        # Predict clear-sky values
        pred = sm.predict(newdata=data)

        # Count violations where GHI exceeds prediction
        violations = np.sum(data['GHI'] - pred > 0)

        # Penalized MSE — large penalty per violation beyond ntol
        mse = np.sum((data['GHI'] - pred) ** 2) + 1_000_000 * (violations - ntol)
        return mse

    # Optimize starting from current model coefficients
    opt = minimize(
        loss_function,
        x0=sm.model.params.values,
        args=(sm, data),
        method='Nelder-Mead' # gradient-free, mirrors R's default optim()
    )

    # Update model with optimal parameters
    sm.update(dict(zip(sm.model.params.index, opt.x)))
    
    return sm

def clearsky_outliers(x, Ct, date, threshold=0.0001, quiet=False):
    """
    Detect outliers in the clear sky model. 
    If x<0, imputed to be equal to min(x) for that day.
    If x>max(Ct), imputed to be Ct*(1-threshold).
    If x is NA, imputed to be the mean(x) for that day.
    """
    # Initialize a dataset
    data = pd.DataFrame({'Ct': Ct, 'x': x})
    
    # Eventually add a date for non-stationary data
    if date is not None:
        data['date'] = pd.to_datetime(date)
        data['Month'] = data['date'].dt.month
        data['Day'] = data['date'].dt.day

    # Detect problems and violations
    outliers_na = data.index[data['x'].isna()].tolist()
    outliers_lo = data.index[data['x'] <= 0].tolist()
    outliers_hi = data.index[data['x'] >= data['Ct']].tolist()
    
    # Complete outliers index
    idx_outliers = list(set(outliers_na + outliers_lo + outliers_hi))

    if len(idx_outliers) == 0:
        if not quiet:
            print("No outliers!")
        data_clean = data.copy()
    else:
        if not quiet:
            pct = (len(idx_outliers) / len(data)) * 100
            print(f"Outliers: {len(idx_outliers)} ({pct:.2f} %)")
            
        # Dataset without outliers
        data_no_outliers = data.drop(index=idx_outliers)
        
        # Imputed dataset
        data_clean = data.copy()
        
        # Convert outlier lists to sets for O(1) lightning-fast lookups in the loop
        set_na = set(outliers_na)
        set_lo = set(outliers_lo)
        set_hi = set(outliers_hi)
        
        # Impute outliers
        for i in idx_outliers:
            df_n = data.loc[i]
            
            # equivalent to R's !missing(date)
            if date is not None: 
                # Filter for the same day and month, but a different date/year
                mask = (
                    (data_no_outliers['Month'] == df_n['Month']) & 
                    (data_no_outliers['Day'] == df_n['Day']) & 
                    (data_no_outliers['date'] != df_n['date'])
                )
                df_day = data_no_outliers[mask]
            else:
                df_day = data_no_outliers
                
            # Impute data depending on outlier type
            if i in set_na:
                data_clean.at[i, 'x'] = df_day['x'].mean()
                
            elif i in set_lo:
                data_clean.at[i, 'x'] = df_day['x'].min()
                
            elif i in set_hi:
                data_clean.at[i, 'x'] = data_clean.at[i, 'Ct'] * (1 - threshold)

    if date is None:
        out_date = None
    else:
        out_date = data.loc[idx_outliers, 'date'].tolist()

    # Fixed error calculations
    mape = (np.abs((data['x'] - data_clean['x']) / data['x'])).mean() * 100
    mse = (data['x'] - data_clean['x']).std(ddof=1)

    return {
        "x": data_clean['x'].values,
        "original": data.loc[idx_outliers, 'x'].values,
        "imputed": data_clean.loc[idx_outliers, 'x'].values,
        "index": idx_outliers,
        "index_type": {"na": outliers_na, "lo": outliers_lo, "hi": outliers_hi},
        "date": out_date,
        "n": len(idx_outliers),
        "MAPE": mape,
        "MSE": mse,
        "threshold": threshold
    }

class SeasonalClearsky(SeasonalModel):
    """Seasonal model for clearsky radiation."""
    def __init__(self, control=None):
        if control is None:
            control = control_seasonalClearsky()

        # Call parent initializer to set up the inherited attributes
        super().__init__(orders=control['orders'], periods=control['periods'])
        
        self.lat = None

        self.version = "1.0.1"
        self.coefficients_orig = None
        self.delta = None
        self._ssf = None
        self._control = control

    def fit(self, x, date, lat, clearsky, H0=None, alt=None,
        optimiser="delta_optimiser"):
        """
        Fit the seasonal model for clearsky radiation (GHI ~ H0 + ...).
    
        Parameters
        ----------
        optimiser : {"delta_optimiser", "constrained"}
            "delta_optimiser" : grid-search for a scalar delta that scales
                the OLS-fitted coefficients to form an upper envelope (fast,
                interpretable, assumes correct model shape).
            "constrained" : refit all model parameters via constrained /
                penalized optimization to enforce GHI <= pred directly
                (more flexible, slower, no post-hoc coefficient scaling).
        """
        _VALID_OPTIMISERS = {"delta_optimiser", "constrained"}
        if optimiser not in _VALID_OPTIMISERS:
            raise ValueError(
                f"`optimiser` must be one of {_VALID_OPTIMISERS}, got '{optimiser}'."
            )
    
        control = self._control
        include_intercept = control["include_intercept"]
        include_trend     = control["include_trend"]
        order_H0          = control["order_H0"]
    
        if clearsky is None:
            raise ValueError("`clearsky` time series must be specified.")
    
        self._ssf = SeasonalSolarFunctions(method='spencer')
        self.lat  = float(np.atleast_1d(lat)[0])
    
        # ------------------------------------------------------------------ #
        # Build dataset
        # ------------------------------------------------------------------ #
        date_series = pd.to_datetime(date)
        data = pd.DataFrame({'date': date_series})
        data['Year']  = data['date'].dt.year
        data['Month'] = data['date'].dt.month
        data['Day']   = data['date'].dt.day
        data['t']     = data['Year'] - data['Year'].max()
        data['n']     = number_of_day(data["date"])
        data['Rt']    = np.asarray(x)
        data['GHI']   = np.asarray(x)   # alias used by clearsky_optimizer
    
        if H0 is None:
            data["H0"] = self._ssf.Hon(data["n"].values, self.lat)
        else:
            data["H0"] = H0
    
        data["clearsky"] = np.asarray(clearsky)
    
        # ------------------------------------------------------------------ #
        # Step 1 — OLS fit: clearsky ~ H0 + harmonics (shared by both paths)
        # ------------------------------------------------------------------ #
        for i in range(2, order_H0 + 1):
            data[f"H0_{i}"] = data["H0"] ** i
    
        ext_regressors = ["H0"] + [f"H0_{i}" for i in range(2, order_H0 + 1)]
        if include_trend:
            ext_regressors.append("t")
    
        super().fit(
            data                = data,
            target_col          = "clearsky",
            time_col            = "n",
            external_regressors = ext_regressors,
            include_intercept   = include_intercept
        )
    
        data["Ct_hat"] = self.predict(newdata=data)
    
        # ------------------------------------------------------------------ #
        # Step 2 — Optimisation (path diverges here)
        # ------------------------------------------------------------------ #
        if optimiser == "delta_optimiser":
            # --- Trim dataset to only what's needed downstream ---
            selected_columns = ["n", "H0", "Rt", "Ct_hat"]
            if include_trend:
                selected_columns.insert(1, "t")
            data = data[selected_columns]
    
            delta_val = clearsky_delta_optimizer(
                data['Rt'],
                data['Ct_hat'] * control['delta0'],
                control['lower'],
                control['upper'],
                control['by'],
                control['ntol']
            )
            if delta_val is None:
                raise ValueError(
                    "delta_optimiser failed to find a valid delta. "
                    "Consider widening [lower, upper] or increasing ntol."
                )
    
            self.coefficients_orig = self._model.params.copy()
            self.delta = delta_val * control['delta0']
    
            std_errors = self._std_errors.copy()
            super().update(self.coefficients_orig * self.delta)
            super().update_std_errors(std_errors * self.delta)
    
        elif optimiser == "constrained":
            # --- Refit all parameters; no scalar delta scaling needed ---
            refitted_model = clearsky_optimizer(
                seasonal_model_Ct = self,
                data              = data,          # must contain 'GHI' column
                ntol              = control['ntol']
            )
            # Propagate refitted parameters back into self
            new_params = refitted_model._model.params
            super().update(dict(zip(new_params.index, new_params.values)))
    
            # delta is not meaningful here but set to 1 for API consistency
            self.coefficients_orig = self._model.params.copy()
            self.delta = 1.0
    
        return self

    def predict(self, n=None, newdata=None, alt=None):
        """Predict method for `seasonalClearsky` object."""
        # Case 1: training fitted values
        if newdata is None and n is None:
            return super().predict()
            
        # Case 2: scalar/array of day-of-year integers
        if newdata is None and n is not None:
            newdata = pd.DataFrame({"n": np.atleast_1d(n)})
        
        # 1. Calculate Extraterrestrial radiation (H0) for the given day(s)
        if 'H0' not in newdata.columns:
            newdata['H0'] = self._ssf.Hon(n=newdata['n'], lat=self.lat, alt=alt)
        
        # 2. Add polynomial terms for H0 if specified in the control dictionary
        order_H0 = self.control["order_H0"]
        if order_H0 > 1:
            for i in range(2, order_H0 + 1):
                newdata[f'H0_{i}'] = newdata['H0'] ** i

        # --- Parent Predict ---
        # We pass the dataframe with H0 and H0_i terms to the parent class.
        # The parent will generate the seasonal sines/cosines based on 'n'
        # and execute the matrix dot product with the shifted 'delta' coefficients.
        return super().predict(newdata, "n")

    def differential(self, n=None, newdata=None):
        """Compute the first derivative of the clear sky model with respect to time (n)."""
        # --- 1. Handle Inputs ---
        if newdata is None and n is None:
            # Return the fitted derivatives from the training data
            return super().differential()
            
        if newdata is None and n is not None:
            newdata = pd.DataFrame({"n": np.atleast_1d(n)})

        # --- 2. Extraterrestrial Radiation ---
        # To compute the derivative of a polynomial, we need base_H0 and dH0/dn.
        H0_base = self._ssf.Hon(n=newdata['n'], lat=self.lat, deriv=False)
        H0_deriv = self._ssf.Hon(n=newdata['n'], lat=self.lat, deriv=True)

        # The 'H0' column passed to the parent needs to be dH0/dn
        newdata['H0'] = H0_deriv

        # --- 3. Polynomial Terms ---
        order_H0 = self.control["order_H0"]
        if order_H0 > 1:
            for i in range(2, order_H0 + 1):
                # Chain Rule: d/dn(H0^i) = i * H0^(i-1) * (dH0/dn)
                newdata[f'H0_{i}'] = i * (H0_base ** (i - 1)) * H0_deriv

        # --- 4. Parent Delegation ---
        # The parent model takes differentiated external regressors, 
        # computes derivatives of internal seasonal sines/cosines, 
        # and multiplies the entire matrix by the standard model coefficients.
        return super().differential(newdata, "n")

    def __str__(self):
        """Print method for `seasonalClearsky` object."""
        msg = "==============================================================\n"
        msg += "--------------------- seasonalClearsky -----------------------\n"
        msg += "==============================================================\n"
        msg += f" - Order:  {self.orders} \n"
        msg += f" - Period: {self.periods} \n"
        msg += " - External regressors: 1 (H0) \n"
        msg += f" - Version: {self.version} \n"
        msg += "-----------------------------------------------------------\n"
        
        if self.model is not None:
            msg += "\n===========================================================================\n"
            msg += str(self.model.summary2(float_format="%.2f").tables[1])
            msg += "\n===========================================================================\n"
        
        return msg

    def __repr__(self):
        return self.__str__()
        
    @property
    def control(self):
        return self._control

    @property
    def std_errors(self):
        """Series with the parameters' standard errors."""
        return self._std_errors

    @property
    def ssf(self):
        return self._ssf
