import numpy as np
import pandas as pd
import warnings
from solarrpy import seasonalModel

# Assuming SeasonalModel is defined in a file named seasonal_model.py
# from seasonal_model import SeasonalModel 
# If running in a single script, ensure the SeasonalModel class is defined above this.

# --- 1. Placeholder Helpers for External R Dependencies ---
# These functions replace the external R packages/functions used in the original code.

def control_seasonal_clearsky():
    """
    Default control parameters dictionary.
    """
    return {
        'order': 1,
        'period': 365,
        'order_H0': 1,
        'include_intercept': True,
        'include_trend': False,
        'delta0': 1.0,
        'lower': 0.0,
        'upper': 2.0,
        'by': 0.01,
        'ntol': 5
    }

def number_of_day(dates):
    """Returns the day of the year (1-366)."""
    return pd.to_datetime(dates).dayofyear

class SeasonalSolarFunctions:
    """
    Mock implementation of seasonalSolarFunctions.
    Calculates Extraterrestrial Radiation (H0).
    """
    def __init__(self, algorithm="spencer"):
        self.algorithm = algorithm
    
    def Hon(self, n, lat, deriv=False):
        """
        Calculate extraterrestrial radiation or its derivative.
        
        Args:
            n (array-like): Day of year.
            lat (float): Latitude.
            deriv (bool): If True, return derivative.
        """
        # 
        # Simplified placeholder logic:
        rads = 2 * np.pi * n / 365
        if deriv:
            return -np.sin(rads) # Dummy derivative
        return 10 + 5 * np.cos(rads) # Dummy H0 value

def clearsky_delta_optimizer(Rt, Ct_hat_scaled, lower, upper, by, ntol):
    """
    Placeholder for the optimization function.
    Returns a scalar scaling factor (delta).
    """
    # Logic to find optimal delta would go here
    return 1.05 

# --- 2. Main Class Conversion ---

class SeasonalClearsky(SeasonalModel):
    """
    R6 implementation for a clear sky seasonal model.
    
    Inherits from SeasonalModel.
    Version: 1.0.1
    """

    def __init__(self, control=None):
        """
        Initialize a `SeasonalClearsky` object.
        
        Args:
            control (dict): Control parameters.
        """
        if control is None:
            control = control_seasonal_clearsky()
            
        # Initialize parent with order and period from control
        super().__init__(order=control.get('order', 1), period=control.get('period', 365))
        
        self.lat = np.nan
        self._control = control
        self._ssf = None  # seasonalSolarFunctions instance
        self._coefficients_orig = None
        self._delta = None
        self._version = "1.0.1"

    @property
    def control(self):
        return self._control

    @property
    def ssf(self):
        return self._ssf

    def fit(self, x, date, lat, clearsky=None):
        """
        Fit the seasonal model for clear sky radiation.

        Args:
            x (array-like): Time series of solar radiation (Rt).
            date (array-like): Time series of dates.
            lat (float): Reference latitude.
            clearsky (array-like): Time series of target clear sky radiation.
        """
        if clearsky is None:
            raise ValueError("`clearsky` time series must be specified.")
            
        # Extract control parameters
        control = self.control
        include_intercept = control.get('include_intercept', True)
        include_trend = control.get('include_trend', False)
        order_H0 = control.get('order_H0', 1)

        # Initialize Solar Functions
        self._ssf = SeasonalSolarFunctions("spencer")
        
        # Store latitude (handle scalar)
        self.lat = lat if np.isscalar(lat) else lat[0]
        
        # Initialize the dataset
        dates = pd.to_datetime(date)
        data = pd.DataFrame({'date': dates})
        
        data['Year'] = dates.year
        data['Month'] = dates.month
        data['Day'] = dates.day
        data['t'] = data['Year'] - data['Year'].max()
        data['n'] = number_of_day(dates)
        data['Rt'] = x
        
        # Compute Extraterrestrial Radiation (H0)
        data['H0'] = self._ssf.Hon(data['n'], self.lat)
        data['clearsky'] = clearsky
        
        # ========================================================================
        # 1. Daily maximum clearsky fit
        # Formula: clearsky ~ H0 + H0^2 + ... + trend + seasonal
        # ========================================================================
        
        base_formula = "clearsky ~ H0"
        
        # Add polynomial terms for H0
        if order_H0 > 1:
            for i in range(2, order_H0 + 1):
                col_name = f"H0_{i}"
                data[col_name] = data['H0'] ** i
                base_formula += f" + {col_name}"
                
        # Add trend
        if include_trend:
            base_formula += " + t"
            
        # Handle Intercept
        if not include_intercept:
            base_formula += " - 1"
            
        # Fit coefficients using parent method
        super().fit(formula=base_formula, data=data)
        
        # Initial fit average clear sky
        data['Ct_hat'] = self.predict(newdata=data)
        
        # ========================================================================
        # 2. Optimization
        # ========================================================================
        
        # Optimize the fit using the helper function
        delta_res = clearsky_delta_optimizer(
            Rt=data['Rt'].values, 
            Ct_hat_scaled=data['Ct_hat'].values * control.get('delta0', 1.0), 
            lower=control.get('lower'), 
            upper=control.get('upper'), 
            by=control.get('by'), 
            ntol=control.get('ntol')
        )
        
        # --- Rename Coefficients (Mirroring R logic) ---
        # The R code renames coefficients to specific domain names (delta_0, delta_extra, etc.)
        
        orig_names = list(self._model.params.index)
        coefs_names = []
        
        # Note: We iterate specifically to match the order expected by the R script logic.
        
        # 1. Intercept -> delta_0
        if include_intercept:
            if "Intercept" in orig_names:
                coefs_names.append("delta_0")
                orig_names.remove("Intercept")
        
        # 2. H0 terms -> delta_extra
        if "H0" in orig_names:
             coefs_names.append("delta_extra1")
             orig_names.remove("H0")
             
        for i in range(2, order_H0 + 1):
            name = f"H0_{i}"
            if name in orig_names:
                coefs_names.append(f"delta_extra{i}")
                orig_names.remove(name)
        
        # 3. Trend -> t
        if include_trend and "t" in orig_names:
            coefs_names.append("t")
            orig_names.remove("t")
            
        # 4. Seasonal Terms -> delta_{original_name}
        # Any remaining terms are assumed to be seasonal sine/cosine terms
        for name in orig_names:
            coefs_names.append(f"delta_{name}")
            
        # Store original coefficients
        self._coefficients_orig = self.coefficients.copy()
        self._delta = delta_res * control.get('delta0', 1.0)
        
        # Update coefficients values by scaling with delta
        new_coefs = self.coefficients * self._delta
        
        # Apply new names if dimensions match
        if len(new_coefs) == len(coefs_names):
            new_coefs.index = coefs_names
            self._model.params = new_coefs
            
            # Update standard errors
            new_se = self.std_errors * self._delta
            new_se.index = coefs_names
            self.update_std_errors(new_se)
        else:
            warnings.warn("Coefficient name mapping mismatch. Updating values without renaming.")
            self.update(new_coefs)

    def predict(self, n=None, newdata=None):
        """
        Predict method for `seasonalClearsky` object.
        """
        control = self.control
        order_H0 = control.get('order_H0', 1)
        
        if newdata is None:
            if n is None:
                # Use fitted values from internal model
                return self._model.predict()
            else:
                # Construct newdata from n, calculating H0
                H0 = self._ssf.Hon(n, self.lat)
                newdata_df = pd.DataFrame({'n': n, 'H0': H0})
                
                # Add polynomial terms
                if order_H0 > 1:
                    for i in range(2, order_H0 + 1):
                        newdata_df[f"H0_{i}"] = newdata_df['H0'] ** i
                
                return self._model.predict(exog=newdata_df)
        else:
            # If newdata is provided
            # If 'n' is present but 'H0' is missing, calculate H0
            if 'H0' not in newdata.columns and 'n' in newdata.columns:
                 newdata['H0'] = self._ssf.Hon(newdata['n'], self.lat)
            
            # Add polynomial terms if missing
            if order_H0 > 1:
                for i in range(2, order_H0 + 1):
                    col = f"H0_{i}"
                    if col not in newdata.columns:
                        newdata[col] = newdata['H0'] ** i
                        
            return self._model.predict(exog=newdata)

    def differential(self, n=None, newdata=None):
        """
        Differential method for `seasonalClearsky` object.
        """
        control = self.control
        order_H0 = control.get('order_H0', 1)
        
        if newdata is None:
            if n is None:
                return self._dmodel.predict()
            else:
                # Calculate derivative of H0
                H0_deriv = self._ssf.Hon(n, self.lat, deriv=True)
                newdata_df = pd.DataFrame({'n': n, 'H0': H0_deriv})
                
                # Apply polynomial chain rule logic from R script
                # R: newdata[[paste0("H0_", i)]] <- i * newdata$H0^(i-1)
                if order_H0 > 1:
                    for i in range(2, order_H0 + 1):
                         newdata_df[f"H0_{i}"] = i * (newdata_df['H0'] ** (i-1))
                
                return self._dmodel.predict(exog=newdata_df)
        else:
            if 'H0' not in newdata.columns and 'n' in newdata.columns:
                newdata['H0'] = self._ssf.Hon(newdata['n'], self.lat, deriv=True)
            
            if order_H0 > 1:
                for i in range(2, order_H0 + 1):
                    col = f"H0_{i}"
                    if col not in newdata.columns:
                        # R logic assumes H0 column holds the base for the power rule here
                        newdata[col] = i * (newdata['H0'] ** (i-1))
                        
            return self._dmodel.predict(exog=newdata)

    def __repr__(self):
        msg = []
        msg.append("----------------------- SeasonalClearsky -----------------------")
        msg.append(f" - Order: {self.order}")
        msg.append(f" - Period: {self.period}")
        msg.append("- External regressors: 1 (H0)")
        msg.append(f"- Version: {self._version}")
        msg.append("--------------------------------------------------------------")
        if self._model is not None:
             msg.append(str(self.coefficients))
        else:
             msg.append("Model not fitted.")
        return "\n".join(msg)