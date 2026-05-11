import pandas as pd
import numpy as np
from scipy.stats import norm
import copy

# --- Imports from existing solarrpy modules (Phase 0 stabilization) ---
# Class names corrected to PascalCase to match the actual definitions in
# solarrpy. The lowercase R-style aliases were never defined.
from .solarTransform import SolarTransform
from .seasonalClearsky import SeasonalClearsky, clearsky_outliers
from .seasonalModel import SeasonalModel
from .ARMA_model import ARMAModel
from .zzz import number_of_day
from .solarModel_internals import solarModel_match_params

# --- Phase 1 calibration (paper Appendix A.1) — re-export the module-level
# entrypoint here so callers can use ``from solarrpy.SolarModel import
# calibrate_window`` per the SoRad reproduction roadmap.  The actual
# implementation lives in :mod:`solarrpy.calibration`.
from .calibration import (  # noqa: E402,F401
    calibrate_window,
    calibrate_all_windows,
    CalibrationWindowResult,
    CalibrationAllWindowsResult,
)

# --- Phase-0 stubs for modules that have not been ported yet ---
# Each name below is referenced inside SolarModel methods. Phase 0 keeps the
# file importable by replacing the missing imports with placeholders that
# raise NotImplementedError if/when the relevant code path runs. Phase 2-5
# will provide real implementations and remove these stubs.
_PHASE_NOT_DONE = (
    "{name} is not implemented in this Phase 0 stabilization snapshot. "
    "It will land in Phase {phase} of the SoRad reproduction roadmap; "
    "see .cursor/plans/sorad_r-to-python_roadmap_*.plan.md."
)


def sGARCH(*args, **kwargs):  # noqa: N802 (preserve R name)
    raise NotImplementedError(
        _PHASE_NOT_DONE.format(name="sGARCH", phase="3 (deferred)")
    )


def solarMixture(*args, **kwargs):  # noqa: N802
    raise NotImplementedError(
        _PHASE_NOT_DONE.format(name="solarMixture", phase="1 (Step 6 EM)")
    )


def solarMoments(*args, **kwargs):  # noqa: N802
    raise NotImplementedError(
        _PHASE_NOT_DONE.format(name="solarMoments", phase="2 (conditional moments)")
    )


def solarMoments_conditional(*args, **kwargs):  # noqa: N802
    raise NotImplementedError(
        _PHASE_NOT_DONE.format(name="solarMoments_conditional", phase="2")
    )


def solarMoments_unconditional(*args, **kwargs):  # noqa: N802
    raise NotImplementedError(
        _PHASE_NOT_DONE.format(name="solarMoments_unconditional", phase="2")
    )


def dsolarGHI(*args, **kwargs):  # noqa: N802
    raise NotImplementedError(
        _PHASE_NOT_DONE.format(name="dsolarGHI", phase="2 (forecast density)")
    )


def qsolar_ghi(*args, **kwargs):
    raise NotImplementedError(
        _PHASE_NOT_DONE.format(name="qsolar_ghi", phase="2 (forecast density)")
    )

"""
SolarModel Class

The `SolarModel` class in Python implements a comprehensive workflow for modeling solar radiation data.
It provides methods to fit seasonal models, detect and handle outliers, perform data transformation,
and apply time-series models such as ARMA and GARCH for residual analysis. The class can use both
seasonal and Gaussian Mixture Models (GMMs) to capture the distributional properties of the data.

Overview:
---------
This class orchestrates the end-to-end fitting process for a solar radiation model:

- **Seasonal Modeling:** Fits seasonal cycles to the input data, optionally on both clear-sky and actual GHI.
- **Outlier Detection:** Identifies and treats outliers in the solar data using robust procedures.
- **Transformation:** Transforms the data (e.g., via power/Box–Cox or log) for improved statistical properties.
- **Time Series Modeling:** Fits ARMA models to the mean and GARCH models to the residual variance.
- **Mixture Modeling:** Fits a Gaussian Mixture Model to capture non-Gaussian features in the standardized residuals.
- **Likelihood Methods:** Computes (log-)likelihoods and related probability quantities.

Usage Example:
--------------
```python
# Model specification (typically via a SolarModelSpec object)
spec = SolarModelSpec()
spec.set_mean_model(arOrder=1, maOrder=1)
spec.specification("Bologna")

# Model initialization and fit
model = SolarModel(spec)
model.fit()

# Accessing results
summary = model.summary()
loglikelihood = model.loglik
fitted_data = model.data

# Updating coefficients and filtering new data
params = model.coefficients
model.update(params)
model.filter()

# Fit a model with realized clear sky
spec.control['stochastic_clearsky'] = True
model = SolarModel(spec)
model.fit()
```

Notes:
------
- The class is designed to be flexible and modular, supporting different fitting steps controlled by the specification.
- It is versioned (see `self.version`).
- Most core components are set up in the constructor and filled out during fitting.

Version: 1.0.1
"""

class SolarModel:
    def __init__(self, spec):
        """
        Initialize the SolarModel.
        :param spec: A SolarModelSpec object.
        """
        self.version = "1.0.1"
        
        # 1. Create seasonal skeleton (366 days to handle leap years)
        dates = pd.date_range(start="2020-01-01", end="2020-12-31", freq='D')
        self._seasonal_data = pd.DataFrame({
            'date': dates,
            'Month': dates.month,
            'Day': dates.day,
            'n': number_of_day(dates) # 1-based day of year
        }).drop(columns=['date']).sort_values(['Month', 'Day'])

        # 2. Setup internal data
        # spec.data is assumed to be a Pandas DataFrame from the previous translation
        self._data = spec.data.copy()
        if 'H0' in self._data.columns:
            self._data = self._data.drop(columns=['H0'])
        
        self._data['loglik'] = 0.0
        self._monthly_data = pd.DataFrame({
            'Month': range(1, 13),
            'Yt_tilde_uncond': 0.0,
            'sigma_uncond': 1.0
        })
        
        self._loglik = None
        self._spec = copy.deepcopy(spec)
        
        # 3. Component placeholders
        # These would be instances of the classes we translated previously
        self._transform = SolarTransform(alpha=0, beta=1, link=self._spec.transform['link'])
        self._seasonal_model_ct = None
        self._seasonal_model_yt = None
        self._arma = None
        self._seasonal_variance = None
        self._garch = None
        self._nm_model = None
        
        self._outliers = None
        self._moments = {'conditional': None, 'unconditional': None}
        self.interpolated = False

        self._hessian = None
        self._jacobian = None

    # =================================================================
    # The Fit Pipeline
    # =================================================================

    def fit(self):
        """Execute the full estimation pipeline."""
        self.fit_seasonal_model_ct()
        self.compute_risk_drivers()
        self.fit_transform()
        self.fit_seasonal_model_yt()
        self.fit_monthly_mean()
        self.fit_arma()
        self.fit_seasonal_variance()
        self.fit_garch()
        self.fit_nm_model()
        self.update_nm_classification()
        self.update_moments()
        self.update_logLik()

    def fit_seasonal_model_ct(self):
        """Fit clear sky model and predict Ct envelope."""
        control = self._spec.clearsky
        # Train filter
        train_mask = self._data['isTrain'] & (self._data['weights'] != 0)
        train_data = self._data[train_mask]
        
        # Initialize and fit
        sm_ct = SeasonalClearsky(control=control)
        sm_ct.fit(
            x=train_data[self._spec.target],
            date=train_data['date'],
            lat=self._spec.coords['lat'],
            clearsky=train_data['clearsky'],
            optimiser="constrained",
        )
        
        # Predict for whole dataset
        self._data['Ct'] = sm_ct.predict(newdata=self._data)
        self._seasonal_model_ct = sm_ct

    def compute_risk_drivers(self):
        """Calculate risk drivers (Xt) and handle outliers."""
        target = self._spec.target
        data = self._data
        
        # Outlier detection (Using your previously translated function)
        env_limit = data['clearsky'] if self._spec.stochastic_clearsky else data['Ct']
        self._outliers = clearsky_outliers(data['GHI'], env_limit, date=data['date'], quiet=self._spec.quiet)
        
        # Update GHI with imputed values
        data['GHI'] = self._outliers['x']
        
        # Compute Transformed Variable (Xt)
        data['Xt'] = self._transform.X(data[target], env_limit)
    
        # Update data
        self._data['Xt'] = data['Xt']

    def fit_transform(self):
        """Optimize alpha/beta for the solar transform."""
        control = self._spec.transform
        train_data = self._data[self._data['isTrain'] & (self._data['weights'] != 0)]
        
        # Fit parameters. BoundTransform.fit only consumes ``Xt`` and an
        # ``epsilon`` scaling factor — the legacy ``min_pos`` / ``max_pos``
        # controls have no consumer in the current Python translation. The
        # paper's prescription (App. A.1.2) is implemented inside
        # BoundTransform.fit itself; ``threshold`` here maps to its epsilon.
        epsilon = control.get('threshold', 1e-3)
        params = self._transform.fit(train_data['Xt'], epsilon=epsilon)
        self._transform.update(params['alpha'], params['beta'])
        
        # Post-process Xt to avoid boundary issues in ARMA-GARCH
        self._data['Xt'] = np.where(self._data['Xt'] <= params['Xt_min'], params['Xt_min'] * (1 + control['delta']), self._data['Xt'])
        self._data['Xt'] = np.where(self._data['Xt'] >= params['Xt_max'], params['Xt_max'] * (1 - control['delta']), self._data['Xt'])
        
        # Calculate Yt (The Gaussianized variable)
        self._data['Yt'] = self._transform.Y(self._transform.X_prime(self._data['Xt']))

    def fit_seasonal_model_yt(self):
        """Fit seasonal mean to Yt."""
        control = self._spec.seasonal_mean
        target = self._spec.target
        train_data = self._data[self._data['isTrain'] & (self._data['weights'] != 0)]
        
        sm_yt = SeasonalModel(orders=control['order'], periods=control['period'])
        # Note: SeasonalModel.fit signature is fit(data, target_col, time_col,
        # external_regressors, include_intercept). It does not accept
        # include_trend; that flag is handled at the seasonalClearsky layer.
        sm_yt.fit(train_data, target_col='Yt', time_col='n', include_intercept=control['include_intercept'])
        
        self._data['Yt_bar'] = sm_yt.predict(self._data)
        self._data['Yt_tilde'] = self._data['Yt'] - self._data['Yt_bar']
        self._data[f'{target}_bar'] = self._transform.iRY(self._data['Yt_bar'], self._data['Ct'])
        
        # Store seasonal model for Yt
        self._seasonal_model_yt = sm_yt

    def fit_monthly_mean(self):
        """Center Yt_tilde by month."""
        if self._spec.seasonal_mean['monthly_mean']:
            train_data = self._data[self._data['isTrain'] & (self._data['weights'] != 0)]
            monthly_mean = train_data.groupby('Month')['Yt_tilde'].mean().reset_index()
            monthly_mean.columns = ['Month', 'Yt_tilde_uncond']
            
            # Merge back
            self._data = self._data.merge(monthly_mean, on='Month', how='left')
            self._data['Yt_tilde'] = self._data['Yt_tilde'] - self._data['Yt_tilde_uncond']
            
            # Update monthly tracking
            self._monthly_data['Yt_tilde_uncond'] = monthly_mean['Yt_tilde_uncond']

    def fit_arma(self):
        """Fit ARMA model to deseasonalized series."""
        control = self._spec.mean_model
        train_data = self._data[self._data['isTrain']]
        
        arma = ARMAModel(
            ar_order=control['arOrder'],
            ma_order=control['maOrder'],
            include_intercept=control['include_intercept'],
        )
        arma.fit(train_data['Yt_tilde'])
        
        self._data['Yt_tilde_hat'] = arma.filter(self._data['Yt_tilde'])
        self._data['eps'] = self._data['Yt_tilde'] - self._data['Yt_tilde_hat']
        self._arma = arma

    def fit_seasonal_variance(self):
        """Fit seasonal variance to ARMA residuals."""
        control = self._spec.seasonal_variance
        train_data = self._data[self._data['isTrain'] & (self._data['weights'] != 0)].copy()
        train_data['eps2'] = train_data['eps']**2
        
        sv = SeasonalModel(orders=control['order'], periods=control['period'])
        # See note in fit_seasonal_model_yt about include_trend.
        sv.fit(train_data, target_col='eps2', time_col='n', include_intercept=True)
        
        self._data['sigma_bar'] = np.sqrt(sv.predict(self._data))
        # Compute standardized residuals
        self._data['eps_tilde'] = self._data['eps'] / self._data['sigma_bar']
        self._seasonal_variance = sv
        
        # Compute monthly corrective variance
        self.fit_monthly_variance()
        # Correction of parameters to ensure unitary variance
        self.correct_seasonal_variance()

    def fit_monthly_variance(self):
        if self._spec.seasonal_variance['monthly_mean']:
            train_data = self._data[self._data['isTrain'] & (self._data['weights'] != 0)]
            monthly_sd = train_data.groupby('Month')['eps_tilde'].std().reset_index()
            monthly_sd.columns = ['Month', 'sigma_uncond']
            
            self._data = self._data.merge(monthly_sd, on='Month', how='left', suffixes=('', '_new'))
            self._data['sigma_uncond'] = self._data['sigma_uncond_new']
            self._data = self._data.drop(columns=['sigma_uncond_new'])
            
            self._data['eps_tilde'] = self._data['eps'] / (self._data['sigma_bar'] * self._data['sigma_uncond'])
            self._monthly_data['sigma_uncond'] = monthly_sd['sigma_uncond']

    def correct_seasonal_variance(self):
        """Ensure standardized residuals have exactly unitary variance."""
        if self._spec.seasonal_variance['correction']:
            train_data = self._data[self._data['isTrain'] & (self._data['weights'] != 0)]
            correction_factor = train_data['eps_tilde'].var()
            
            # Update the seasonal variance model coefficients
            new_coefs = self._seasonal_variance.coefficients * correction_factor
            self._seasonal_variance.update(new_coefs)
            
            # Recalculate
            self._data['sigma_bar'] = np.sqrt(self._seasonal_variance.predict(self._data))
            self._data['eps_tilde'] = self._data['eps'] / (self._data['sigma_bar'] * self._data['sigma_uncond'])

    def fit_garch(self):
        """Fit GARCH model to standardize residuals further.

        Phase 0 note
        ------------
        The SoRad scope of the reproduction roadmap explicitly excludes
        sGARCH (see plan section "Scope exclusions"), so the only supported
        path here is ``garch_variance == False``: ``sigma`` is left at 1.0
        and ``u_tilde`` is the seasonal-variance-standardized ``eps_tilde``.
        Activating ``garch_variance`` would require Phase 3 of the prior
        plan (sGARCH translation), which has been deferred.
        """
        control = self._spec
        if control.garch_variance:
            # Building the sGARCH stub raises NotImplementedError; only do
            # so if the user actually opted into GARCH.
            train_data = self._data[self._data['isTrain'] & (self._data['weights'] != 0)]
            garch = sGARCH(arch_order=control.variance_model['archOrder'],
                           garch_order=control.variance_model['garchOrder'])
            garch.fit(train_data['eps_tilde'], weights=train_data['weights'])
            self._data['sigma'] = np.sqrt(garch.filter(self._data['eps_tilde']))
            self._data['u_tilde'] = self._data['eps_tilde'] / self._data['sigma']
            self._garch = garch
        else:
            self._data['sigma'] = 1.0
            self._data['u_tilde'] = self._data['eps_tilde']
            self._garch = None

    def fit_nm_model(self):
        """Fit Gaussian Mixture Model to standardized residuals."""
        control = self._spec.mixture_model
        train_data = self._data[self._data['isTrain'] & (self._data['weights'] != 0)]
        
        # Gaussian Mixture fitted only on train data
        nm = solarMixture(components=2, abstol=control['abstol'], maxit=control['maxit'], maxrestarts=control['maxrestarts'])

        # Logic for matching expectation or variance (Duffie & Beckman style constraint)
        if control['match_expectation'] or control['match_variance']:
            # 1. Calculate monthly empirical moments
            # We aggregate by Month to get a vector of 12 targets
            target = train_data.groupby('Month')['u_tilde'].agg(
                mu_target='mean', 
                var_target='var'
            ).reset_index()
            
            # 2. Handle non-empirical matching (target 0 mean, 1 variance)
            if not control['match_empiric']:
                target['mu_target'] = 0.0
                target['var_target'] = 1.0
            
            # 3. Prepare target vectors (Use None for R's NA to signal the fitter to ignore the constraint)
            mu_target = target['mu_target'].values if control['match_expectation'] else None
            var_target = target['var_target'].values if control['match_variance'] else None

            # 4. Fit with constraints
            nm.fit(
                x=train_data['u_tilde'], 
                date=train_data['date'], 
                weights=train_data['weights'],
                mu_target=mu_target, 
                var_target=var_target,
                method=control['method']
            )
        else:
            # Standard fit without moment constraints
            nm.fit(
                x=train_data['u_tilde'], 
                date=train_data['date'], 
                weights=train_data['weights'],
                method=control['method']
            )
        
        self._nm_model = nm

    def update(self, params=None):
        """
        Update the parameters of all sub-models within the SolarModel.

        Args:
            params (dict): A nested dictionary containing parameter updates.
                           Should follow the structure of self.coefficients.
        """
        if params is None:
            print("`params` is missing, nothing to update!")
            return

        if not isinstance(params, dict):
            params = solarModel_match_params(params, self.coefficients)

        if 'params' in params:
            self._transform.update(params['alpha'], params['beta'])

        # Update clear sky model
        self._seasonal_model_ct.update(params['seasonal_model_ct'])
        # Update seasonal mean model
        self._seasonal_model_yt.update(params['seasonal_model_yt'])
        # Update ARMA model
        self._arma.update(params['ARMA'])
        # Update GARCH model
        self._garch.update(params['GARCH'])
        # Update seasonal variance model
        self._seasonal_variance.update(params['seasonal_variance'])

        # Update Gaussian Mixture Model (NM_model)
        means = np.column_stack((params['NM_mu_up'], params['NM_mu_dw']))
        sds = np.column_stack((params['NM_sd_up'], params['NM_sd_dw']))
        probs = np.column_stack((params['NM_p_up'], 1 - params['NM_p_up']))

        # Update the NM sub-model
        self._nm_model.update(means=means, sd=sds, p=probs)

        # Invalidate the stored total log-likelihood
        self._loglik = None

    def update_moments(self):
        """
        Update the internal state for conditional and unconditional moments.
        Calculates expected values and variances across the time series.
        """
        # 1. Update conditional moments
        # Uses the 'data' property which includes seasonal joins
        # and the 'spec' property for model constraints.
        self._moments['conditional'] = solarMoments_conditional(
            self._data, 
            control_model=self._spec
        )

        # 2. Update unconditional moments
        # Passes the ARMA and GARCH sub-models explicitly as required by the R logic.
        self._moments['unconditional'] = solarMoments_unconditional(
            self.data,
            ARMA=self._arma,
            GARCH=self._garch
        )

    def update_logLik(self):
        """
        Update individual and total log-likelihood values.
        Synchronizes the observation-level densities with the aggregate model score.
        """
        # 1. Update individual log-likelihoods
        # self.log_lik() is assumed to return a pandas Series or numpy array 
        # aligned with the internal data.
        self._data['loglik'] = self.logLik()

        # 2. Update the total log-likelihood attribute
        # We use skipna=False to match R's default behavior where 
        # a single NA makes the whole sum NA. 
        # If your pipeline handles NaNs elsewhere, you can use the default skipna=True.
        self._loglik = self.data['loglik'].sum() 

    def update_risk_drivers(self):
        """
        Synchronize the clear sky envelope, risk drivers, and solar transform.
        Ensures that the Gaussianized series Yt is updated following any 
        parameter changes in the astronomical components.
        """
        # 1. Update clear sky (Ct)
        # Using the data property to ensure seasonal features are joined for prediction.
        # self._seasonal_model_ct is an instance of SeasonalClearsky.
        self._data['Ct'] = self._seasonal_model_ct.predict(newdata=self.data)

        # 2. Update risk driver (Xt)
        # This method handles outlier detection and computes the raw risk driver.
        # It is assumed to be defined as part of the SolarModel class.
        self.compute_risk_drivers()

        # 3. Update solar transform and compute Yt
        # This method re-fits alpha/beta and computes the Gaussianized variable.
        # It is assumed to be defined as part of the SolarModel class.
        self.fit_transform()

    def update_nm_classification(self, run_filter=False):
        """
        Update the classification of the Bernoulli random variable and Gaussian components.
        
        Args:
            run_filter (bool): If True, updates the mixture classification 
                               via the NM model before re-classifying.
        """
        if run_filter:
            # Update mixture classification via the latent model
            self._nm_model.filter()

        # 1. Prepare data and clean existing classification columns
        data = self._data.copy()
        cols_to_drop = [c for c in ['B', 'z1', 'z2'] if c in data.columns]
        if cols_to_drop:
            data = data.drop(columns=cols_to_drop)

        # 2. Classify the series month-by-month
        monthly_results = []
        
        for month in range(1, 13):
            # Filter data for the specific month
            month_mask = data['Month'] == month
            df_month = data[month_mask].copy()
            
            if df_month.empty:
                continue
                
            # Classify using the sub-model for this month (0-indexed list)
            # Assuming self._nm_model.models is the list of 12 sub-models
            nm_results = self._nm_model.models[month - 1].classify(df_month['u_tilde'])
            
            # Extract Bernoulli variable and Gaussian components
            df_month['B'] = nm_results['B1']
            df_month['z1'] = nm_results['z1']
            df_month['z2'] = nm_results['z2']
            
            # Select only the key and new columns for the join
            monthly_results.append(df_month[['date', 'B', 'z1', 'z2']])

        # 3. Consolidate results
        all_classifications = pd.concat(monthly_results)
        
        # Merge back to the internal data storage using 'date' as the key
        self._data = self._data.merge(all_classifications, on='date', how='left')

    def filter(self, fit=True):
        """
        Re-propagate values through the modeling pipeline to synchronize
        residuals with the current model parameters.
        
        Args:
            fit (bool): If True, re-estimates monthly corrections and 
                        applies seasonal variance scaling.
        """
        target = self._spec.target
        data = self._data

        # 1. Update Seasonal Mean of Yt
        # Uses the seasonalModel_Yt sub-model
        data['Yt_bar'] = self._seasonal_model_yt.predict(newdata=self.data)
        
        # 2. Update Seasonal Mean of the Target variable (GHI/clearsky)
        # Using the inverse transform iRY(Yt_bar, Ct)
        data[f'{target}_bar'] = self._transform.iRY(data['Yt_bar'], data['Ct'])
        
        # 3. Compute Deseasonalized Series
        data['Yt_tilde'] = data['Yt'] - data['Yt_bar']
        
        # Re-estimate monthly mean correction if specified
        if fit:
            self.fit_monthly_mean()
            
        # 4. Update ARMA Filter and Residuals
        data['Yt_tilde_hat'] = self._arma.filter(data['Yt_tilde'])
        data['eps'] = data['Yt_tilde'] - data['Yt_tilde_hat']
        
        # 5. Update Seasonal Volatility and Standardized Residuals
        # Predict sigma_bar (seasonal standard deviation)
        data['sigma_bar'] = np.sqrt(self._seasonal_variance.predict(data))
        data['eps_tilde'] = data['eps'] / data['sigma_bar']
        
        # Re-estimate monthly variance and apply parameter correction if specified
        if fit:
            self.fit_monthly_variance()
            self.correct_seasonal_variance()
            
        # 6. Update GARCH Conditional Volatility
        if self._spec.garch_variance:
            data['sigma'] = np.sqrt(self._garch.filter(data['eps_tilde']))
        else:
            data['sigma'] = 1.0
            
        # 7. Final Standardized Residuals
        data['u_tilde'] = data['eps_tilde'] / data['sigma']

    def compute_moments(self, t_now, t_hor, theta=0, quiet=False):
        """
        Compute conditional moments for a sequence of future dates.

        Renamed from ``moments`` in Phase 0 to resolve a method/property name
        collision: ``moments`` is also a ``@property`` returning the cached
        ``{'conditional': ..., 'unconditional': ...}`` dict.

        Args:
            t_now (str/datetime): Reference date for the forecast.
            t_hor (list/pd.DatetimeIndex): Sequence of future dates to forecast.
            theta (float): Shift parameter for the mixture model.
            quiet (bool): If True, suppresses diagnostic messages.

        Returns:
            pd.DataFrame: A combined DataFrame of moments for each date in t_hor.
        """
        # 1. Ensure t_hor is iterable (handles single date or list)
        if isinstance(t_hor, (str, pd.Timestamp)):
            t_hor_list = [t_hor]
        else:
            t_hor_list = t_hor
            
        # Map over the horizon dates and call solarMoments
        results = []
        for t in t_hor_list:
            moment_row = solarMoments(
                t_now, t, self._data, self._arma, self._garch, self._nm_model,
                self._transform, theta=theta, quiet=quiet
            )
            results.append(moment_row)
        
        return pd.concat(results, ignore_index=True) if results else pd.DataFrame()

    def VaR(self, moments=None, t_now=None, t_hor=None, theta=0, ci=0.05):
        """
        Compute Value at Risk for the solar model.
        """
        # 1. Generate moments if not provided
        if moments is None:
            # Generate sequence from t_now + 1 day to t_hor
            start_date = pd.to_datetime(t_now) + pd.Timedelta(days=1)
            t_seq = pd.date_range(start=start_date, end=t_hor, freq='D')
            
            # Call the internal moments method (renamed in Phase 0).
            moments = self.compute_moments(t_now, t_seq, theta=theta, quiet=False)

        # 2. Add realized GHI from internal data
        # self.data is the property that joins seasonal features
        moments = pd.merge(
            moments, 
            self.data[['date', 'GHI']], 
            on='date', 
            how='left'
        )

        # 3. Initialize VaR column
        moments['VaR'] = np.nan
        
        # 4. Loop through each time step to calculate the quantile
        for i in range(len(moments)):
            row = moments.iloc[i]
            
            # Define the CDF of the Gaussian Mixture in the Yt space
            # F(x) = p1 * Phi((x - mu1)/sigma1) + (1 - p1) * Phi((x - mu0)/sigma0)
            def cdf_y(x):
                term1 = norm.cdf(x, row['M_Y1'], row['S_Y1']) * row['p1']
                term2 = norm.cdf(x, row['M_Y0'], row['S_Y0']) * (1 - row['p1'])
                return term1 + term2
                
            # Compute Value at Risk using the inverse transform logic
            moments.at[i, 'VaR'] = qsolar_ghi(
                ci, 
                row['Ct'], 
                row['alpha'], 
                row['beta'], 
                cdf_y, 
                link=self.spec.transform['link']
            )

        # 5. Identify Violations
        # et = 1 if GHI < VaR, else 0
        moments['et'] = (moments['GHI'] < moments['VaR']).astype(int)

        # 6. Ensure date components are present for the final selection
        moments['Year'] = moments['date'].dt.year
        moments['Month'] = moments['date'].dt.month
        moments['Day'] = moments['date'].dt.day

        # Select and return relevant variables
        return moments[['date', 'Year', 'Month', 'Day', 'GHI', 'VaR', 'et']]

    def logLik(self, moments=None, target="Yt", quasi=False):
        """
        Compute the log-likelihood of the model.
        
        Args:
            moments (pd.DataFrame): Optional pre-computed moments. Defaults to conditional moments.
            target (str): The variable space for the likelihood ("Yt" or "GHI").
            quasi (bool): If True, uses a single Gaussian density (pseudo-likelihood).
                        If False, uses the full Gaussian Mixture density.
                        
        Returns:
            pd.Series: Vector of log-likelihood values for each observation.
        """
        # 1. Handle default moments
        if moments is None:
            moments = self.moments['conditional']
        
        # 2. Join with internal data to get target values, weights, and training flags
        # Equivalent to dplyr::left_join(moments, private$..data, by = "date")
        m = moments.merge(
            self._data[['date', target, 'weights', 'isTrain']], 
            on='date', 
            how='left'
        )
        
        # Initialize log-likelihood array
        loglik = np.zeros(len(m))
        
        # 3. Density calculation in physical space (GHI)
        if target == "GHI":
            # This involves the change-of-variables (Jacobian) handled by dsolar_ghi
            for i in range(len(m)):
                row = m.iloc[i]
                
                # Skip zero-weighted training points (R logic: moments$loglik[i] <- 0)
                if row['isTrain'] and row['weights'] == 0:
                    continue
                    
                if quasi:
                    # Normal density function for the quasi-likelihood
                    pdf_func = lambda x: norm.pdf(x, row['e_Yt'], row['sd_Yt'])
                else:
                    # Mixture density function: p1*N(mu1, sd1) + (1-p1)*N(mu0, sd0)
                    pdf_func = lambda x: (norm.pdf(x, row['M_Y1'], row['S_Y1']) * row['p1'] + 
                                        norm.pdf(x, row['M_Y0'], row['S_Y0']) * (1 - row['p1']))
                
                # dsolarGHI handles the transformation and Jacobian adjustment
                val = dsolarGHI(
                    row['GHI'], 
                    row['Ct'], 
                    row['alpha'], 
                    row['beta'], 
                    pdf_func, 
                    link=self.spec.transform['link']
                )
                loglik[i] = np.log(val)
                
        # 4. Density calculation in standardized space (Yt)
        else:
            if quasi:
                # Standardize and compute log-pdf: log(phi(z) / sigma)
                z = (m[target] - m['e_Yt']) / m['sd_Yt']
                loglik = norm.logpdf(z) - np.log(m['sd_Yt'])
            else:
                # Mixture log-likelihood in Yt space
                z1 = (m[target] - m['M_Y1']) / m['S_Y1']
                z0 = (m[target] - m['M_Y0']) / m['S_Y0']
                
                # Density = (p1 * phi(z1)/s1) + ((1-p1) * phi(z0)/s0)
                density = (norm.pdf(z1) / m['S_Y1'] * m['p1'] + 
                        norm.pdf(z0) / m['S_Y0'] * (1 - m['p1']))
                loglik = np.log(density)
                
        return pd.Series(loglik, index=m['date'])

    def __str__(self):
        """String representation of the SolarModel."""
        # Get specifications
        data = self._spec.dates['data']
        train = self._spec.dates['train']
        test = self._spec.dates['test']
        
        train_perc = f"{train['perc']*100:.2f}"
        test_perc = f"{test['perc']*100:.2f}"
        
        lines = [
            f"--------------------- solarModel ({self.place}) ---------------------",
            f"Model: ARMA({self._ARMA.order[0]}, {self._ARMA.order[1]})-GARCH({self._GARCH.order[0]}, {self._GARCH.order[1]})",
            f"Target: {self._spec.target}",
            f" Coordinates: (Lat: {self._spec.coords['lat']}, Lon: {self._spec.coords['lon']}, Alt: {self._spec.coords['alt']})",
            f" Dates: {data['from']} - {data['to']}",
            f" Observations: {data['nobs']}",
            "---------------------------------------------------------------",
            f"Train dates: {train['from']} - {train['to']} ({train['nobs']} points ~ {train_perc}%)",
            f" Test dates: {test['from']} - {test['to']} ({test['nobs']} points ~ {test_perc}%)",
            "---------------------------------------------------------------",
            f"Log-Likelihood: {self._loglik:.8f}" if self._loglik else "Log-Likelihood: Not computed",
            f"Interpolated: {self._interpolated}",
            f"Version: {self.VERSION}"
        ]
        
        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    # =================================================================
    # Active Properties (Getters)
    # =================================================================

    @property
    def place(self):
        """Character, optional name of the location considered."""
        return self._spec.place

    @property
    def model_name(self):
        """Character, model's name combining link function and model orders."""
        link = self._spec.transform['link']
        ar, ma = self._arma.order
        arch, garch = self._garch.order
        return f"{link}-ARMA({ar}, {ma})GARCH({arch}, {garch})"

    @property
    def monthly_data(self):
        """A data frame that contains monthly parameters, including GMM coefficients."""
        if self._nm_model is not None and not self._nm_model.is_empty:
            return pd.merge(self._monthly_data, self._nm_model.coefficients, on='Month', how='left')
        return self._monthly_data

    @property
    def seasonal_data(self):
        """A data frame containing seasonal and monthly parameters."""
        return pd.merge(self._seasonal_data, self.monthly_data, on='Month', how='left')

    @property
    def data(self):
        """The full fitted data frame joined with seasonal and monthly parameters."""
        # Equivalent to dplyr::select(self$seasonal_data, -n)
        seasonal = self.seasonal_data.drop(columns=['n'])
        # Equivalent to dplyr::left_join(private$..data, seasonal_data, by = c("Month", "Day"))
        return pd.merge(self._data, seasonal, on=['Month', 'Day'], how='left')

    @property
    def loglik(self):
        """The log-likelihood computed on train data."""
        return self._loglik

    @property
    def spec(self):
        """The specification object governing the model."""
        return self._spec

    @property
    def location(self):
        """A data frame with coordinates and location metadata."""
        loc_dict = {
            'place': [self.place],
            'target': [self._spec.target]
        }
        # Merge the coordinates dictionary (lat, lon, alt)
        loc_dict.update({k: [v] for k, v in self._spec.coords.items()})
        return pd.DataFrame(loc_dict)

    @property
    def transform(self): return self._transform

    @property
    def seasonal_model_Ct(self): return self._seasonal_model_ct

    @property
    def seasonal_model_Yt(self): return self._seasonal_model_yt

    @property
    def ARMA(self): return self._arma

    @property
    def seasonal_variance(self): return self._seasonal_variance

    @property
    def GARCH(self): return self._garch

    @property
    def NM_model(self): return self._nm_model

    @property
    def moments(self):
        """Returns conditional and unconditional moments with extra metadata columns."""
        moments_dict = self._moments.copy()
        
        # Metadata to add to moments for downstream forecasting/plotting
        target_bar_col = f"{self._spec.target}_bar"
        extra_cols = ['date', target_bar_col, 'Ct', 'p1']
        
        # Extract from self.data (which handles the seasonal joins)
        data_extra = self.data[extra_cols].copy()
        data_extra['alpha'] = self._transform.alpha
        data_extra['beta'] = self._transform.beta
        
        # Merge metadata into both moment types if they are populated
        if moments_dict['conditional'] is not None:
            moments_dict['conditional'] = pd.merge(
                moments_dict['conditional'], data_extra, on='date', how='left'
            )
            
        if moments_dict['unconditional'] is not None:
            moments_dict['unconditional'] = pd.merge(
                moments_dict['unconditional'], data_extra, on='date', how='left'
            )
            
        return moments_dict

    @property
    def coefficients(self):
        """Aggregates all sub-model parameters into a structured nested dictionary."""
        
        # Gaussian Mixture Model (NM) parameter formatting (12 months)
        nm_coefs = self._nm_model.coefficients
        
        def format_nm(values, prefix):
            # Creates a single-row DataFrame with columns like mu_up_1, mu_up_2...
            cols = [f"{prefix}_{i}" for i in range(1, 13)]
            return pd.DataFrame([values.values], columns=cols)

        params = {
            "location": self.location,
            "params": pd.DataFrame({'alpha': [self._transform.alpha], 'beta': [self._transform.beta]}),
            "seasonal_model_Ct": pd.DataFrame([self._seasonal_model_ct.coefficients]),
            "seasonal_model_Yt": pd.DataFrame([self._seasonal_model_yt.coefficients]),
            "ARMA": pd.DataFrame([self._arma.coefficients]),
            "seasonal_variance": pd.DataFrame([self._seasonal_variance.coefficients]),
            "GARCH": pd.DataFrame([self._garch.coefficients]),
            "NM_mu_up": format_nm(nm_coefs['mu1'], "mu_up"),
            "NM_mu_dw": format_nm(nm_coefs['mu2'], "mu_dw"),
            "NM_sd_up": format_nm(nm_coefs['sd1'], "sd_up"),
            "NM_sd_dw": format_nm(nm_coefs['sd2'], "sd_dw"),
            "NM_p_up": format_nm(nm_coefs['p1'], "p")
        }

        # Add stochastic clearsky correlations if they exist
        if hasattr(self._nm_model, 'rho_up') and self._nm_model.rho_up is not None:
            params['rho_up'] = format_nm(nm_coefs['rho_up'], "rho_up")
            
        if hasattr(self._nm_model, 'rho_dw') and self._nm_model.rho_dw is not None:
            params['rho_dw'] = format_nm(nm_coefs['rho_dw'], "rho_dw")

        return params

    @property
    def var_theta(self):
        """
        Computes the Variance-Covariance matrix using the Sandwich Estimator.
        H^-1 * (J' * J) * H^-1
        """
        H = self._hessian
        J = self._jacobian
        
        # solve(H) in R is equivalent to np.linalg.inv(H)
        # crossprod(J, J) is J.T @ J
        H_inv = np.linalg.inv(H)
        sandwich = H_inv @ (J.T @ J) @ H_inv
        
        return sandwich