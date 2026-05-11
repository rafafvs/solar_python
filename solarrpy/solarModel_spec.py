import pandas as pd
import numpy as np
from datetime import datetime
from typing import Optional, Dict, Any, Union

from .seasonalClearsky import clearsky_outliers, control_seasonalClearsky

# Assuming clearsky_outliers and control_seasonal_clearsky are defined in 
# your other modules (as per previous files)
# from .clearsky_utils import clearsky_outliers, control_seasonal_clearsky

class SolarModelSpec:
    def __init__(self):
        """Initialize the specification with default parameters."""
        # Initialize internal state (Private)
        self._place = None
        self._coords = {}
        self._data = None
        self._dates = {}
        self._target = None
        
        self._transform = {}
        self._clearsky = {}
        self._seasonal_mean = {}
        self._mean_model = {}
        self._seasonal_variance = {}
        self._variance_model = {}
        self._mixture_model = {}
        
        self._garch_variance = False
        self._stochastic_clearsky = False
        self._clearsky_threshold = 1.01
        self._quiet = False

        # Set default values via setters
        self.set_params()
        self.set_transform()
        self.set_clearsky()
        self.set_seasonal_mean()
        self.set_mean_model()
        self.set_seasonal_variance()
        self.set_variance_model()
        self.set_mixture_model()

    def specification(self,
                      place: str,
                      data: pd.DataFrame,
                      target: str = "GHI",
                      min_date: Optional[str] = None,
                      max_date: Optional[str] = None,
                      from_date: Optional[str] = None,
                      to_date: Optional[str] = None,
                      coords: Optional[Dict[str, float]] = None):
        """
        Main entry point for configuring the data and modeling window.
        """
        if target not in ["GHI", "clearsky"]:
            raise ValueError("target must be 'GHI' or 'clearsky'")

        # Handle Date Logic (equivalent to R's missing check)
        df = data.copy()
        df['date'] = pd.to_datetime(df['date'])

        actual_min = df['date'].min()
        actual_max = df['date'].max()

        limit_min = pd.to_datetime(min_date) if min_date else actual_min
        limit_max = pd.to_datetime(max_date) if max_date else actual_max
        train_from = pd.to_datetime(from_date) if from_date else actual_min
        train_to = pd.to_datetime(to_date) if to_date else actual_max

        # Filtering
        df = df[(df['date'] >= limit_min) & (df['date'] <= limit_max)].copy()

        # Pre-processing
        df['clearsky'] = df['clearsky'] * self._clearsky_threshold
        
        # Impute outliers
        outliers = clearsky_outliers(df['GHI'], df['clearsky'], date=df['date'], threshold=0, quiet=self._quiet)
        df['GHI'] = outliers['x']

        # Labeling train/test sets
        df['isTrain'] = (df['date'] >= train_from) & (df['date'] <= train_to)
        df['weights'] = df['isTrain'].astype(float)
        df['weights'] = df['weights'] / df['weights'].sum()

        # Stats calculation
        nobs_train = df['isTrain'].sum()
        perc_train = nobs_train / len(df)
        nobs_test = (~df['isTrain']).sum()
        perc_test = nobs_test / len(df)

        self._dates = {
            "data": {"from": limit_min, "to": limit_max, "nobs": len(df), "perc": 1.0},
            "train": {"from": train_from, "to": train_to, "nobs": nobs_train, "perc": perc_train},
            "test": {"from": train_to, "to": limit_max, "nobs": nobs_test, "perc": perc_test}
        }

        # Store results (Mimicking R attributes)
        self._place = place
        # Coords resolution order: explicit kwarg > df.attrs['coords'] >
        # legacy ad-hoc attribute > unknown sentinel. Storing "NA" strings
        # broke downstream float() conversions; use NaN floats instead so
        # callers can detect missing values in a uniform way.
        if coords is not None:
            self._coords = dict(coords)
        elif "coords" in getattr(df, "attrs", {}):
            self._coords = dict(df.attrs["coords"])
        elif hasattr(df, "coords"):
            self._coords = dict(getattr(df, "coords"))
        else:
            self._coords = {"lat": float("nan"), "lon": float("nan"), "alt": float("nan")}
        self._data = df
        self._target = target

    # =================================================================
    # Control/Setter Methods
    # =================================================================

    def set_params(self, stochastic_clearsky=False, clearsky_threshold=1.01, quiet=False):
        self._stochastic_clearsky = stochastic_clearsky
        self._clearsky_threshold = clearsky_threshold
        self._quiet = quiet

    def set_transform(self, min_pos=1, max_pos=1, link="invgumbel", delta=0.05, threshold=0.01):
        valid_links = ["invgumbel", "gumbel", "logis", "norm"]
        if link not in valid_links:
            raise ValueError(f"link must be one of {valid_links}")
        self._transform = {
            "min_pos": min_pos, "max_pos": max_pos, "delta": delta,
            "threshold": threshold, "link": link
        }

    def set_clearsky(self, control=None, **kwargs):
        """Configure the clear-sky control dict.

        Parameters
        ----------
        control : dict or None
            Pre-built control dict (typically the return value of
            ``control_seasonalClearsky``). If None, a default is built using
            ``control_seasonalClearsky()`` with any keyword overrides.
        **kwargs : optional
            Forwarded to ``control_seasonalClearsky`` when ``control`` is None.
        """
        if control is None:
            self._clearsky = control_seasonalClearsky(**kwargs)
        else:
            self._clearsky = control

    def set_seasonal_mean(self, order=1, period=365, include_trend=False, include_intercept=True, monthly_mean=False):
        self._seasonal_mean = {
            "order": order, "period": period, "include_trend": include_trend,
            "include_intercept": include_intercept, "monthly_mean": monthly_mean
        }

    def set_mean_model(self, arOrder=1, maOrder=0, include_intercept=False):
        self._mean_model = {"arOrder": arOrder, "maOrder": maOrder, "include_intercept": include_intercept}

    def set_seasonal_variance(self, order=1, period=365, include_trend=False, correction=False, monthly_mean=False):
        self._seasonal_variance = {
            "order": order, "period": period, "correction": correction,
            "include_trend": include_trend, "monthly_mean": monthly_mean
        }

    def set_variance_model(self, archOrder=1, garchOrder=1, garch_variance=True):
        self._garch_variance = garch_variance and not (archOrder == 0 and garchOrder == 0)
        self._variance_model = {"archOrder": archOrder, "garchOrder": garchOrder}

    def set_mixture_model(self, abstol=1e-20, match_expectation=False, match_variance=False,
                          match_empiric=False, maxit=5000, maxrestarts=500):

        self._mixture_model = {
            "abstol": abstol, "match_expectation": match_expectation,
            "match_variance": match_variance, "match_empiric": match_empiric,
            "maxit": maxit, "maxrestarts": maxrestarts
        }

    # =================================================================
    # Active Bindings (Properties)
    # =================================================================

    @property
    def place(self): return self._place

    @property
    def target(self): return self._target

    @property
    def coords(self): return self._coords

    @property
    def dates(self): return self._dates

    @property
    def data(self): return self._data

    @property
    def transform(self): return self._transform

    @property
    def clearsky(self): return self._clearsky

    @property
    def seasonal_mean(self): return self._seasonal_mean

    @property
    def mean_model(self): return self._mean_model

    @property
    def seasonal_variance(self): return self._seasonal_variance

    @property
    def variance_model(self): return self._variance_model

    @property
    def mixture_model(self): return self._mixture_model

    @property
    def garch_variance(self): return self._garch_variance

    @property
    def clearsky_threshold(self): return self._clearsky_threshold

    @property
    def stochastic_clearsky(self): return self._stochastic_clearsky

    @property
    def quiet(self): return self._quiet

    # =================================================================
    # Print Representation
    # =================================================================

    def __repr__(self):
        """Replicates the colorful R print output."""
        # ANSI Colors
        PURPLE = '\033[95m'
        BLUE = '\033[94m'
        GREEN = '\033[92m'
        RED = '\033[91m'
        END = '\033[0m'
        
        def msg_col(val): return f"{GREEN}{val}{END}" if val else f"{RED}{val}{END}"

        output = []
        if self._place:
            d = self._dates['data']
            tr = self._dates['train']
            te = self._dates['test']
            
            output.append(f"{'-'*21} solarModel ({PURPLE}{self._place}{END}) {'-'*21}")
            output.append(f"{PURPLE}Target{END}: {self._target}")
            output.append(f"{PURPLE}Coordinates{END}: ({PURPLE}Lat{END}: {self._coords.get('lat')}, "
                          f"{PURPLE}Lon{END}: {self._coords.get('lon')}, {PURPLE}Alt{END}: {self._coords.get('alt')})")
            output.append(f"{PURPLE}Observations{END}: {d['nobs']}")
            output.append("-" * 63)
            output.append(f"{BLUE}Dates{END}: {d['from'].date()} - {d['to'].date()}")
            output.append(f"  - {BLUE}Train{END}: {tr['from'].date()} - {tr['to'].date()} ({tr['nobs']} points ~ {tr['perc']:.2%})")
            output.append(f"  - {BLUE}Test{END}: {te['from'].date()} - {te['to'].date()} ({te['nobs']} points ~ {te['perc']:.2%})")

        output.append(f"{'-'*24}{PURPLE} Specification {END}{'-'*24}")
        model_name = f"S({self._seasonal_mean['order']}, {self._seasonal_variance['order']})-ARMA({self._mean_model['arOrder']}, {self._mean_model['maOrder']})"
        if self._garch_variance:
            model_name += f"-GARCH({self._variance_model['archOrder']}, {self._variance_model['archOrder']})"
        output.append(model_name)
        output.append(f" - Stochastic clearsky: {msg_col(self._stochastic_clearsky)}")
        
        output.append(f"{'-'*25}{PURPLE} Mean Models {END}{'-'*25}")
        output.append(f" - Trend: {msg_col(self._seasonal_mean['include_trend'])}")
        output.append(f" - Intercept: {msg_col(self._seasonal_mean['include_intercept'])}")
        
        return "\n".join(output)