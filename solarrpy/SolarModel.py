import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Tuple
import copy
from solarrpy import solarMixture, sGARCH, solarMoments


class SolarModel:
    """
    Solar Model class for comprehensive solar radiation modeling.
    
    This class implements a complete solar model that includes fitting seasonal models,
    detecting outliers, performing transformations, and applying time-series models
    such as ARMA and GARCH.
    
    Example:
    --------
    # Model specification
    spec = SolarModelSpec()
    spec.set_mean_model(ar_order=1, ma_order=1)
    spec.specification("Bologna")
    
    # Model fit
    model = SolarModel(spec)
    model.fit()
    print(model)
    """
    
    VERSION = "1.0.1"
    
    def __init__(self, spec):
        """
        Initialize a SolarModel.
        
        Parameters:
        -----------
        spec : SolarModelSpec
            Specification object containing model setup
        """
        # Create seasonal data by month and day for a leap year (366 days)
        dates = pd.date_range(start='2020-01-01', end='2020-12-31', freq='D')
        seasonal_data = pd.DataFrame({
            'date': dates,
            'Month': dates.month,
            'Day': dates.day,
            'n': self._number_of_day(dates)
        })
        seasonal_data = seasonal_data.drop('date', axis=1)
        seasonal_data = seasonal_data.sort_values(['Month', 'Day'])
        
        # Initialize private attributes
        self._data = spec.data.drop('H0', axis=1, errors='ignore').copy()
        self._data['loglik'] = 0.0
        self._seasonal_data = seasonal_data
        self._monthly_data = pd.DataFrame({
            'Month': range(1, 13),
            'Yt_tilde_uncond': 0.0,
            'sigma_uncond': 1.0
        })
        self._loglik = None
        self._spec = copy.deepcopy(spec)
        self._transform = SolarTransform(alpha=0, beta=1, link=self._spec.transform['link'])
        self._seasonal_model_Ct = None
        self._seasonal_model_Yt = None
        self._ARMA = None
        self._seasonal_variance = None
        self._GARCH = None
        self._NM_model = None
        self._outliers = None
        self._hessian = None
        self._jacobian = None
        self._moments = {'conditional': None, 'unconditional': None}
        self._interpolated = False
    
    @staticmethod
    def _number_of_day(dates):
        """Calculate day of year for given dates."""
        if isinstance(dates, pd.DatetimeIndex):
            return dates.dayofyear
        return pd.to_datetime(dates).dayofyear
    
    # ==================== Public Fitting Methods ====================
    
    def fit(self):
        """
        Initialize and fit a SolarModel object given the specification.
        
        This method sequentially fits all model components:
        1. Clear sky seasonal model
        2. Risk drivers
        3. Solar transform
        4. Seasonal mean
        5. Monthly mean correction
        6. ARMA model
        7. Seasonal variance
        8. GARCH variance
        9. Mixture model
        """
        print("Fitting solar model...")
        
        # 1) Clear sky
        self.fit_seasonal_model_Ct()
        
        # 2) Risk-driver
        self.compute_risk_drivers()
        
        # 3) Solar transform
        self.fit_transform()
        
        # 4) Seasonal mean
        self.fit_seasonal_model_Yt()
        
        # Center the mean to be exactly equal to zero
        self.fit_monthly_mean()
        
        # 5) ARMA model
        self.fit_ARMA()
        
        # 6) Seasonal variance
        self.fit_seasonal_variance()
        
        # 7) GARCH variance
        self.fit_GARCH()
        
        # 8) Mixture model
        self.fit_NM_model()
        self.update_NM_classification()
        
        # Update the moments
        self.update_moments()
        
        # Update the log-likelihoods
        self.update_logLik()
        
        print("Model fitting complete!")
    
    def fit_seasonal_model_Ct(self):
        """
        Initialize and fit a seasonal clear sky model.
        """
        control = self._spec.clearsky
        data = self._data[self._data['isTrain'] & (self._data['weights'] != 0)].copy()
        
        # Initialize a seasonal model for clear sky radiation
        seasonal_model_Ct = SeasonalClearsky(control=control)
        
        # Fit the seasonal model
        seasonal_model_Ct.fit(
            GHI=data[self._spec.target].values,
            date=data['date'].values,
            lat=self._spec.coords['lat'],
            clearsky=data['clearsky'].values
        )
        
        # Add seasonal clear sky to data
        self._data['Ct'] = seasonal_model_Ct.predict(self._data['n'].values)
        
        # Store the model
        self._seasonal_model_Ct = copy.deepcopy(seasonal_model_Ct)
    
    def compute_risk_drivers(self):
        """
        Compute the risk drivers and impute observations that exceed clear sky level.
        """
        target = self._spec.target
        control = self._spec
        data = self._data.copy()
        transform = self._transform
        
        if control.stochastic_clearsky:
            # Risk driver
            data['Xt'] = transform.X(data[target].values, data['clearsky'].values)
            
            # Detect and impute outliers
            outliers = clearsky_outliers(
                data['GHI'].values,
                data['clearsky'].values,
                date=data['date'].values,
                quiet=control.quiet
            )
            
            # Update GHI
            data['GHI'] = outliers['x']
            
            # Risk driver
            data['Xt'] = transform.X(data['GHI'].values, data['clearsky'].values)
        else:
            # Detect and impute outliers
            outliers = clearsky_outliers(
                data['GHI'].values,
                data['Ct'].values,
                date=data['date'].values,
                quiet=control.quiet
            )
            
            # Update GHI
            data['GHI'] = outliers['x']
            
            # Add computed risk driver
            data['Xt'] = transform.X(data['GHI'].values, data['Ct'].values)
        
        # Update data
        self._data['Xt'] = data['Xt']
        
        # Store outliers data
        self._outliers = outliers
    
    def fit_transform(self):
        """
        Fit the parameters of the SolarTransform object.
        """
        control = self._spec.transform
        data = self._data[self._data['isTrain'] & (self._data['weights'] != 0)].copy()
        outliers = self._outliers
        
        # Fit transformation parameters
        params = self._transform.fit(
            data['Xt'].values,
            threshold=control['threshold'],
            min_pos=control['min_pos'],
            max_pos=control['max_pos']
        )
        
        # Update transform parameters
        self._transform.update(params['alpha'], params['beta'])
        
        data = self._data.copy()
        
        # Rescale minimum to avoid extreme values
        idx_Xt_min = data['Xt'] <= params['Xt_min']
        data.loc[idx_Xt_min, 'Xt'] = params['Xt_min'] * (1 + control['delta'])
        
        # Rescale maximum
        idx_Xt_max = data['Xt'] >= params['Xt_max']
        data.loc[idx_Xt_max, 'Xt'] = params['Xt_max'] * (1 - control['delta'])
        
        # Store the index of imputed values
        if outliers is None:
            outliers = {'index_type': {}}
        if 'index_type' not in outliers:
            outliers['index_type'] = {}
        
        outliers['index_type']['transform'] = list(np.where(idx_Xt_min | idx_Xt_max)[0])
        
        # Update outliers index and dates
        all_indices = set(outliers.get('index', []))
        all_indices.update(outliers['index_type']['transform'])
        outliers['index'] = sorted(list(all_indices))
        outliers['date'] = data.iloc[outliers['index']]['date'].values
        
        # Compute the transformed variable
        data['Yt'] = self._transform.Y(self._transform.X_prime(data['Xt'].values))
        
        # Update private data
        self._data['Yt'] = data['Yt']
        self._outliers = outliers
    
    def fit_seasonal_model_Yt(self):
        """
        Fit a seasonal model on the transformed variable (Yt) and compute deseasonalized series.
        """
        target = self._spec.target
        control = self._spec.seasonal_mean
        data = self._data.copy()
        
        # Train data
        data_train = data[data['isTrain'] & (data['weights'] != 0)].copy()
        
        # Initialize the seasonal model for Yt
        seasonal_model_Yt = SeasonalModel(
            order=control['order'],
            period=control['period']
        )
        
        # Fit with appropriate formula
        seasonal_model_Yt.fit(
            y=data_train['Yt'].values,
            t=data_train['n'].values,
            include_intercept=control['include_intercept'],
            include_trend=control['include_trend']
        )
        
        # Compute Yt_bar
        data['Yt_bar'] = seasonal_model_Yt.predict(data['n'].values)
        
        # Compute Yt_tilde
        data['Yt_tilde'] = data['Yt'] - data['Yt_bar']
        
        # Store seasonal model for Yt
        self._seasonal_model_Yt = copy.deepcopy(seasonal_model_Yt)
        
        # Update data
        self._data['Yt_tilde'] = data['Yt_tilde']
        self._data['Yt_bar'] = data['Yt_bar']
        self._data[f'{target}_bar'] = self._transform.iRY(
            data['Yt_bar'].values,
            data['Ct'].values
        )
    
    def fit_monthly_mean(self):
        """
        Correct the deseasonalized series by subtracting its monthly mean.
        """
        control = self._spec
        
        if control.seasonal_mean['monthly_mean']:
            data = self._data.copy()
            
            # Train data
            train_data = data[data['isTrain'] & (data['weights'] != 0)]
            
            # Compute monthly unconditional mean
            monthly_mean = train_data.groupby('Month')['Yt_tilde'].mean().reset_index()
            monthly_mean.columns = ['Month', 'Yt_tilde_uncond']
            
            # Add unconditional mean to the dataset
            data = data.merge(monthly_mean, on='Month', how='left')
            
            # Update Yt_tilde
            data['Yt_tilde'] = data['Yt_tilde'] - data['Yt_tilde_uncond']
            
            # Update monthly data
            self._monthly_data['Yt_tilde_uncond'] = monthly_mean['Yt_tilde_uncond'].values
            
            # Update private data
            self._data['Yt_tilde'] = data['Yt_tilde']
    
    def fit_ARMA(self):
        """
        Fit an ARMA model on Yt_tilde and compute ARMA residuals (eps).
        """
        control = self._spec.mean_model
        data = self._data.copy()
        
        # Train data
        data_train = data[data['isTrain']].copy()
        
        # Initialize an ARMA model
        ARMA = ARMAModel(
            ar_order=control['arOrder'],
            ma_order=control['maOrder'],
            include_intercept=control['include_intercept']
        )
        
        # Fit ARMA model
        ARMA.fit(data_train['Yt_tilde'].values)
        
        # Fitted Yt_tilde
        data['Yt_tilde_hat'] = ARMA.filter(data['Yt_tilde'].values)
        
        # Fitted residuals
        data['eps'] = data['Yt_tilde'] - data['Yt_tilde_hat']
        
        # Store ARMA model
        self._ARMA = copy.deepcopy(ARMA)
        
        # Update data
        self._data['Yt_tilde_hat'] = data['Yt_tilde_hat']
        self._data['eps'] = data['eps']
    
    def fit_seasonal_variance(self):
        """
        Fit a seasonal model on ARMA squared residuals and compute deseasonalized residuals.
        """
        control = self._spec.seasonal_variance
        data = self._data.copy()
        
        # Train data
        data_train = data[data['isTrain'] & (data['weights'] != 0)].copy()
        
        # Squared residuals
        data_train['eps2'] = data_train['eps'] ** 2
        
        # Initialize the seasonal model for eps2
        seasonal_variance = SeasonalModel(
            order=control['order'],
            period=control['period']
        )
        
        # Fit seasonal variance
        seasonal_variance.fit(
            y=data_train['eps2'].values,
            t=data_train['n'].values,
            include_intercept=True,
            include_trend=control['include_trend']
        )
        
        # Fitted seasonal standard deviation
        self._data['sigma_bar'] = np.sqrt(seasonal_variance.predict(data['n'].values))
        
        # Compute standardized residuals
        self._data['eps_tilde'] = self._data['eps'] / self._data['sigma_bar']
        
        # Store seasonal variance model
        self._seasonal_variance = copy.deepcopy(seasonal_variance)
        
        # Compute monthly corrective variance
        self.fit_monthly_variance()
        
        # Correction of parameters to ensure unitary variance
        self.correct_seasonal_variance()
    
    def fit_monthly_variance(self):
        """
        Correct the standardized series by dividing by its monthly standard deviation.
        """
        control = self._spec.seasonal_variance
        data = self._data.copy()
        
        if control['monthly_mean']:
            # Train data
            train_data = data[data['isTrain'] & (data['weights'] != 0)]
            
            # Compute monthly unconditional std deviation
            monthly_data = train_data.groupby('Month')['eps_tilde'].std().reset_index()
            monthly_data.columns = ['Month', 'sigma_uncond']
            
            # Add unconditional variance to monthly data
            self._monthly_data['sigma_uncond'] = monthly_data['sigma_uncond'].values
            
            # Merge with data
            data = data.merge(monthly_data, on='Month', how='left')
            
            # Updated standardized residuals
            self._data['eps_tilde'] = self._data['eps'] / (
                self._data['sigma_bar'] * data['sigma_uncond']
            )
    
    def correct_seasonal_variance(self):
        """
        Correct the parameters of the seasonal variance to ensure unitary variance.
        """
        control = self._spec.seasonal_variance
        data = self.data.copy()
        
        # Initialize correction factor
        if not hasattr(self._seasonal_variance, 'extra_params'):
            self._seasonal_variance.extra_params = {}
        self._seasonal_variance.extra_params['correction'] = 1.0
        
        if control['correction']:
            seasonal_variance = copy.deepcopy(self._seasonal_variance)
            
            # Update train data
            data_train = data[data['isTrain'] & (data['weights'] != 0)]
            
            # Calculate variance correction factor
            correction_factor = np.var(
                data_train['eps'] / (data_train['sigma_bar'] * data_train['sigma_uncond'])
            )
            seasonal_variance.extra_params['correction'] = correction_factor
            
            # Correct parameters
            corrected_coefs = seasonal_variance.coefficients * correction_factor
            seasonal_variance.update(corrected_coefs)
            
            # Update sigma_bar
            self._data['sigma_bar'] = np.sqrt(
                seasonal_variance.predict(self._data['n'].values)
            )
            
            # Updated standardized residuals
            self._data['eps_tilde'] = self._data['eps'] / (
                self._data['sigma_bar'] * self.data['sigma_uncond']
            )
            
            # Update seasonal variance model
            self._seasonal_variance = seasonal_variance
    
    def fit_GARCH(self):
        """
        Fit a GARCH model on the deseasonalized residuals (eps_tilde).
        """
        control = self._spec
        data = self._data.copy()
        
        # Train data
        data_train = data[data['isTrain'] & (data['weights'] != 0)]
        
        # GARCH specification
        GARCH_spec = control.variance_model
        
        # Initialize a GARCH model
        GARCH_model = sGARCH(
            arch_order=GARCH_spec['archOrder'],
            garch_order=GARCH_spec['garchOrder']
        )
        
        # Control for garch variance
        if control.garch_variance:
            # Fit the model
            GARCH_model.fit(
                data_train['eps_tilde'].values,
                weights=data_train['weights'].values
            )
            
            # Fitted std. deviation
            data['sigma'] = np.sqrt(GARCH_model.filter(data['eps_tilde'].values))
            
            # Fitted standardized residuals
            data['u_tilde'] = data['eps_tilde'] / data['sigma']
        else:
            # No GARCH
            data['sigma'] = 1.0
            data['u_tilde'] = data['eps_tilde']
        
        # Store GARCH parameters
        self._GARCH = copy.deepcopy(GARCH_model)
        
        # Update data
        self._data['sigma'] = data['sigma']
        self._data['u_tilde'] = data['u_tilde']
    
    def fit_NM_model(self):
        """
        Initialize and fit a Gaussian Mixture model.
        """
        control = self._spec.mixture_model
        outliers = self._outliers
        
        # Train data
        data_train = self._data[
            self._data['isTrain'] & (self._data['weights'] != 0)
        ].copy()
        
        # Gaussian Mixture fitted only on train data
        NM_model = SolarMixture(
            components=2,
            abstol=control['abstol'],
            maxit=control['maxit'],
            maxrestarts=control['maxrestarts']
        )
        
        # Match moments if specified
        if control['match_expectation'] or control['match_variance']:
            # Target empirical moments
            target = data_train.groupby('Month')['u_tilde'].agg([
                ('mu_target', 'mean'),
                ('var_target', 'var')
            ]).reset_index()
            
            # Match or not empirical moments
            if not control['match_empiric']:
                target['mu_target'] = 0.0
                target['var_target'] = 1.0
            
            # Match or not expectation
            if not control['match_expectation']:
                target['mu_target'] = None
            
            # Match or not variance
            if not control['match_variance']:
                target['var_target'] = None
            
            NM_model.fit(
                x=data_train['u_tilde'].values,
                date=data_train['date'].values,
                weights=data_train['weights'].values,
                mu_target=target['mu_target'].values,
                var_target=target['var_target'].values,
                method=control['method']
            )
        else:
            NM_model.fit(
                x=data_train['u_tilde'].values,
                date=data_train['date'].values,
                weights=data_train['weights'].values,
                method=control['method']
            )
        
        # Store Gaussian Mixture parameters
        self._NM_model = copy.deepcopy(NM_model)
    
    def update_NM_classification(self, filter=False):
        """
        Update the classification of the Bernoulli random variable.
        
        Parameters:
        -----------
        filter : bool
            When True, update mixture classification before classifying
        """
        if filter:
            self._NM_model.filter()
        
        # Complete data
        data = self._data.copy()
        
        # Remove existing classification columns
        data = data.drop(['B', 'z1', 'z2'], axis=1, errors='ignore')
        
        # Classify the series by month
        df_list = []
        for nmonth in range(1, 13):
            df_m = data[data['Month'] == nmonth].copy()
            
            # Classify
            classification = self._NM_model.model[nmonth].classify(
                df_m['u_tilde'].values
            )
            
            df_m['B'] = classification['B1']
            df_m['z1'] = classification['z1']
            df_m['z2'] = classification['z2']
            
            df_list.append(df_m[['date', 'B', 'z1', 'z2']])
        
        # Merge classifications back
        df_classifications = pd.concat(df_list, ignore_index=True)
        data = data.merge(df_classifications, on='date', how='left')
        
        # Update private data
        self._data['B'] = data['B']
        self._data['z1'] = data['z1']
        self._data['z2'] = data['z2']
    
    def update(self, params):
        """
        Update the parameters inside the model.
        
        Parameters:
        -----------
        params : dict or array-like
            Model parameters to update
        """
        if params is None:
            print("`params` is missing, nothing to update!")
            return
        
        if not isinstance(params, dict):
            params = solarModel_match_params(params, self.coefficients)
        
        # Update transform parameters
        self._transform.update(params['params']['alpha'], params['params']['beta'])
        
        # Update clear sky model
        self._seasonal_model_Ct.update(params['seasonal_model_Ct'])
        
        # Update seasonal mean model
        self._seasonal_model_Yt.update(params['seasonal_model_Yt'])
        
        # Update ARMA mean model
        self._ARMA.update(params['ARMA'])
        
        # Update seasonal variance model
        self._seasonal_variance.update(params['seasonal_variance'])
        
        # Update GARCH variance model
        self._GARCH.update(params['GARCH'])
        
        # Update Gaussian Mixture model
        means = np.column_stack([
            params['NM_mu_up'],
            params['NM_mu_dw']
        ])
        sd = np.column_stack([
            params['NM_sd_up'],
            params['NM_sd_dw']
        ])
        p = np.column_stack([
            params['NM_p_up'],
            1 - params['NM_p_up']
        ])
        
        self._NM_model.update(means=means, sd=sd, p=p)
        
        # Set log-likelihood to NA
        self._loglik = None
    
    def update_moments(self):
        """Update the moments inside the model."""
        # Update conditional moments
        self._moments['conditional'] = solarMoments_conditional(
            self.data,
            control_model=self._spec
        )
        
        # Update unconditional moments
        self._moments['unconditional'] = solarMoments_unconditional(
            self.data,
            ARMA=self._ARMA,
            GARCH=self._GARCH
        )
    
    def update_logLik(self):
        """Update the log-likelihood inside the model."""
        # Update the log-likelihoods
        self._data['loglik'] = self.logLik()
        
        # Update total log-likelihood
        self._loglik = self._data['loglik'].sum()
    
    def update_risk_drivers(self):
        """Update the clear sky and risk drivers."""
        # Update clear sky
        self._data['Ct'] = self._seasonal_model_Ct.predict(self._data['n'].values)
        
        # Update risk driver
        self.compute_risk_drivers()
        
        # Update solar transform and compute Yt
        self.fit_transform()
    
    def filter(self, fit=True):
        """
        Filter the time series when new parameters are supplied via update().
        
        Parameters:
        -----------
        fit : bool
            When True, monthly mean and variances will be re-estimated
        """
        control = self._spec
        target = self._spec.target
        
        # Update seasonal mean of Yt
        self._data['Yt_bar'] = self._seasonal_model_Yt.predict(self._data['n'].values)
        
        # Update seasonal mean of target variable
        self._data[f'{target}_bar'] = self._transform.iRY(
            self._data['Yt_bar'].values,
            self._data['Ct'].values
        )
        
        # Update Yt_tilde
        self._data['Yt_tilde'] = self._data['Yt'] - self._data['Yt_bar']
        
        # Fit the corrective mean
        if fit:
            self.fit_monthly_mean()
        
        # Update Yt_tilde_hat
        self._data['Yt_tilde_hat'] = self._ARMA.filter(self._data['Yt_tilde'].values)
        
        # Update ARMA residuals
        self._data['eps'] = self._data['Yt_tilde'] - self._data['Yt_tilde_hat']
        
        # Update seasonal std. deviation
        self._data['sigma_bar'] = np.sqrt(
            self._seasonal_variance.predict(self._data['n'].values)
        )
        
        # Update standardized residuals
        self._data['eps_tilde'] = self._data['eps'] / self._data['sigma_bar']
        
        if fit:
            # Compute corrective monthly variance
            self.fit_monthly_variance()
            
            # Correct the seasonal variance
            self.correct_seasonal_variance()
        
        # Update GARCH standard deviation
        if self._spec.garch_variance:
            self._data['sigma'] = np.sqrt(
                self._GARCH.filter(self._data['eps_tilde'].values)
            )
        else:
            self._data['sigma'] = 1.0
        
        # Fitted standardized residuals
        self._data['u_tilde'] = self._data['eps_tilde'] / self._data['sigma']
    
    def Moments(self, t_now, t_hor, theta=0, quiet=False):
        """
        Compute the conditional moments.
        
        Parameters:
        -----------
        t_now : str
            Today's date
        t_hor : str or list
            Horizon date(s)
        theta : float
            Shift parameter for the mixture
        quiet : bool
            Suppress verbose messages
        
        Returns:
        --------
        pd.DataFrame
            DataFrame with computed moments
        """
        if isinstance(t_hor, str):
            t_hor = [t_hor]
        
        results = []
        for t_h in t_hor:
            moment = solarMoments(
                t_now, t_h,
                self.data, self._ARMA,
                self._GARCH, self._NM_model,
                self._transform,
                theta=theta,
                quiet=quiet
            )
            results.append(moment)
        
        return pd.concat(results, ignore_index=True)
    
    def VaR(self, moments=None, t_now=None, t_hor=None, theta=0, ci=0.05):
        """
        Compute Value at Risk for the solar model.
        
        Parameters:
        -----------
        moments : pd.DataFrame, optional
            Pre-computed moments
        t_now : str
            Today's date
        t_hor : str
            Horizon date
        theta : float
            Shift parameter
        ci : float
            Confidence interval (one tail)
        
        Returns:
        --------
        pd.DataFrame
            DataFrame with VaR computations
        """
        # Compute moments if not provided
        if moments is None:
            t_seq = pd.date_range(start=t_now, end=t_hor, freq='D')[1:]
            moments = pd.concat([
                self.Moments(t_now, str(t.date()), theta, quiet=False)
                for t in t_seq
            ])
        
        # Add realized GHI
        moments = moments.merge(
            self.data[['date', 'GHI']],
            on='date',
            how='left'
        )
        
        # Initialize VaR
        moments['VaR'] = np.nan
        
        link = self._spec.transform['link']
        
        for i in range(len(moments)):
            df_n = moments.iloc[i]
            
            # CDF of Yt
            def cdf_Y(x):
                return pmixnorm(
                    x,
                    mean=[df_n['M_Y1'], df_n['M_Y0']],
                    sd=[df_n['S_Y1'], df_n['S_Y0']],
                    alpha=[df_n['p1'], 1 - df_n['p1']]
                )
            
            # Value at Risk
            moments.loc[i, 'VaR'] = qsolarGHI(
                ci, df_n['Ct'], df_n['alpha'],
                df_n['beta'], cdf_Y, link=link
            )
        
        # Violations of VaR
        moments['et'] = (moments['GHI'] < moments['VaR']).astype(int)
        
        # Select only relevant variables
        return moments[['date', 'Year', 'Month', 'Day', 'GHI', 'VaR', 'et']]
    
    def logLik(self, moments=None, target="Yt", quasi=False):
        """
        Compute the log-likelihood of the model.
        
        Parameters:
        -----------
        moments : pd.DataFrame, optional
            Dataset containing the moments
        target : str
            Target variable ("Yt" or "GHI")
        quasi : bool
            When True, compute pseudo-likelihood with Gaussian link
        
        Returns:
        --------
        np.ndarray
            Log-likelihood values
        """
        # Default argument
        if moments is None:
            moments = self._moments['conditional'].copy()
        
        # Add target and weights
        moments = moments.merge(
            self._data[['date', target, 'weights', 'isTrain']],
            on='date',
            how='left'
        )
        
        # Compute log-likelihood
        moments['loglik'] = 0.0
        
        link = self._spec.transform['link']
        
        if target == "GHI":
            for i in range(len(moments)):
                df_n = moments.iloc[i]
                
                if df_n['weights'] == 0 and df_n['isTrain']:
                    moments.loc[i, 'loglik'] = 0.0
                    continue
                
                if quasi:
                    # Normal density
                    def pdf_quasi(x):
                        return dnorm(x, df_n['e_Yt'], df_n['sd_Yt'])
                    
                    # Quasi log-likelihood
                    moments.loc[i, 'loglik'] = np.log(
                        dsolarGHI(
                            df_n['GHI'], df_n['Ct'],
                            df_n['alpha'], df_n['beta'],
                            pdf_quasi, link=link
                        )
                    )
                else:
                    # True mixture density
                    def pdf_Yt(x):
                        return dmixnorm(
                            x,
                            mean=[df_n['M_Y1'], df_n['M_Y0']],
                            sd=[df_n['S_Y1'], df_n['S_Y0']],
                            alpha=[df_n['p1'], 1 - df_n['p1']]
                        )
                    
                    # True log-likelihood
                    moments.loc[i, 'loglik'] = np.log(
                        dsolarGHI(
                            df_n['GHI'], df_n['Ct'],
                            df_n['alpha'], df_n['beta'],
                            pdf_Yt, link=link
                        )
                    )
        else:  # target == "Yt"
            if quasi:
                # Standardize the time series
                moments['z'] = (moments['Yt'] - moments['e_Yt']) / moments['sd_Yt']
                
                # Pseudo log-likelihood
                moments['loglik'] = np.log(
                    dnorm(moments['z']) / moments['sd_Yt']
                )
            else:
                # Standardize into components
                moments['z1'] = (moments['Yt'] - moments['M_Y1']) / moments['S_Y1']
                moments['z0'] = (moments['Yt'] - moments['M_Y0']) / moments['S_Y0']
                
                # True mixture log-likelihood
                moments['loglik'] = np.log(
                    dnorm(moments['z1']) / moments['S_Y1'] * moments['p1'] +
                    dnorm(moments['z0']) / moments['S_Y0'] * (1 - moments['p1'])
                )
        
        return moments['loglik'].values
    
    def __repr__(self):
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
    
    # ==================== Properties ====================
    
    @property
    def place(self):
        """Location name."""
        return self._spec.place
    
    @property
    def model_name(self):
        """Model name string."""
        return (f"{self._spec.transform['link']}-"
                f"ARMA({self._ARMA.order[0]}, {self._ARMA.order[1]})"
                f"GARCH({self._GARCH.order[0]}, {self._GARCH.order[1]})")
    
    @property
    def data(self):
        """Get complete data with seasonal and monthly parameters."""
        seasonal_data = self._seasonal_data.drop('n', axis=1, errors='ignore')
        result = self._data.merge(seasonal_data, on=['Month', 'Day'], how='left')
        return result
    
    @property
    def seasonal_data(self):
        """Get seasonal data."""
        return self._seasonal_data.merge(self.monthly_data, on='Month', how='left')
    
    @property
    def monthly_data(self):
        """Get monthly data with mixture parameters."""
        if self._NM_model is not None:
            return self._monthly_data.merge(
                self._NM_model.coefficients,
                on='Month',
                how='left'
            )
        return self._monthly_data
    
    @property
    def loglik(self):
        """Get log-likelihood."""
        return self._loglik
    
    @property
    def spec(self):
        """Get model specification."""
        return self._spec
    
    @property
    def location(self):
        """Get location information."""
        return pd.DataFrame({
            'place': [self.place],
            'target': [self._spec.target],
            **self._spec.coords
        })
    
    @property
    def transform(self):
        """Get solar transform object."""
        return self._transform
    
    @property
    def seasonal_model_Ct(self):
        """Get seasonal clear sky model."""
        return self._seasonal_model_Ct
    
    @property
    def seasonal_model_Yt(self):
        """Get seasonal mean model."""
        return self._seasonal_model_Yt
    
    @property
    def ARMA(self):
        """Get ARMA model."""
        return self._ARMA
    
    @property
    def seasonal_variance(self):
        """Get seasonal variance model."""
        return self._seasonal_variance
    
    @property
    def GARCH(self):
        """Get GARCH model."""
        return self._GARCH
    
    @property
    def NM_model(self):
        """Get Gaussian mixture model."""
        return self._NM_model
    
    @property
    def moments(self):
        """Get conditional and unconditional moments."""
        moments = self._moments.copy()
        
        # Extra data to add
        data_extra = self.data[['date', 'GHI_bar', 'Ct', 'p1']].copy()
        data_extra['alpha'] = self._transform.alpha
        data_extra['beta'] = self._transform.beta
        
        # Conditional moments
        if moments['conditional'] is not None:
            moments['conditional'] = moments['conditional'].merge(
                data_extra, on='date', how='left'
            )
        
        # Unconditional moments
        if moments['unconditional'] is not None:
            moments['unconditional'] = moments['unconditional'].merge(
                data_extra, on='date', how='left'
            )
        
        return moments
    
    @property
    def coefficients(self):
        """Get model parameters as a dictionary."""
        # Extract coefficients from all models
        coefs = {
            'location': self.location,
            'params': pd.DataFrame({
                'alpha': [self._transform.alpha],
                'beta': [self._transform.beta]
            }),
            'seasonal_model_Ct': pd.DataFrame(
                self._seasonal_model_Ct.coefficients
            ),
            'seasonal_model_Yt': pd.DataFrame(
                self._seasonal_model_Yt.coefficients
            ),
            'ARMA': pd.DataFrame(self._ARMA.coefficients),
            'seasonal_variance': pd.DataFrame(
                self._seasonal_variance.coefficients
            ),
            'GARCH': pd.DataFrame(self._GARCH.coefficients),
        }
        
        # Add mixture model parameters
        nm_coefs = self._NM_model.coefficients
        coefs['NM_mu_up'] = pd.DataFrame({
            f'mu_up_{i+1}': [nm_coefs.loc[i, 'mu1']]
            for i in range(12)
        })
        coefs['NM_mu_dw'] = pd.DataFrame({
            f'mu_dw_{i+1}': [nm_coefs.loc[i, 'mu2']]
            for i in range(12)
        })
        coefs['NM_sd_up'] = pd.DataFrame({
            f'sd_up_{i+1}': [nm_coefs.loc[i, 'sd1']]
            for i in range(12)
        })
        coefs['NM_sd_dw'] = pd.DataFrame({
            f'sd_dw_{i+1}': [nm_coefs.loc[i, 'sd2']]
            for i in range(12)
        })
        coefs['NM_p_up'] = pd.DataFrame({
            f'p_{i+1}': [nm_coefs.loc[i, 'p1']]
            for i in range(12)
        })
        
        return coefs
    
    @property
    def var_theta(self):
        """Variance-covariance matrix with robust standard errors."""
        H = self._hessian
        J = self._jacobian
        return np.linalg.inv(H) @ J.T @ J @ np.linalg.inv(H)


# Note: This translation assumes the existence of helper classes:
# - SolarTransform, SeasonalClearsky, SeasonalModel, ARMAModel, 
# - sGARCH, SolarMixture, and various helper functions
# These would need to be translated separately
