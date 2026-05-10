"""
solarrpy — Python translation of the solarr R package.

This package re-exports the public API used by the SoRad / SoRadIDX
reproduction notebooks. Only modules that already exist in the partial
translation are exposed here. Phases that have not been implemented yet
(sGARCH, solarMixture, solarMoments, dsolarGHI, dmixnorm, sorad, soradidx,
buyer) will be added to this re-export list as they land.
"""

# --- Calendar / numerical helpers
from .zzz import number_of_day

# --- Bounded transforms
from .boundTransform import BoundTransform
from .solarTransform import SolarTransform

# --- Solar geometry / clear sky
from .seasonalSolarFunctions import SeasonalSolarFunctions
from .seasonalClearsky import (
    SeasonalClearsky,
    control_seasonalClearsky,
    clearsky_outliers,
    clearsky_optimizer,
    clearsky_delta_optimizer,
)

# --- Seasonal regression
from .seasonalModel import (
    SeasonalModel,
    add_seasonal_features,
    add_seasonal_features_dt,
    seasonalModel_params_to_zeta,
    seasonalModel_params_to_phi,
    seasonalModel_params_to_zeta_jacobian,
)

# --- ARMA mean dynamics
from .ARMA_model import ARMAModel

# --- Gaussian mixture (EM)
from .gaussianMixture import GaussianMixtureModel
from .gaussianMixture_internals import (
    gm_moments,
    gm_moments_match,
    gm_loglik,
    gm_fit,
    gm_fit_moments_match,
)

# --- Continuous-time radiation model + Appendix A.1 internals
from .radiationModel import RadiationModel
from .radiationModel_internals import (
    create_monthly_sequence,
    martingale_method_seasonal,
    reparam_seasonal_function,
    integral_sigma_numeric,
    integral_sigma2_formula,
)

# --- SolarModel orchestrator + spec
from .solarModel_spec import SolarModelSpec
from .SolarModel import SolarModel

__all__ = [
    # Helpers
    "number_of_day",
    # Transforms
    "BoundTransform", "SolarTransform",
    # Solar geometry / clear sky
    "SeasonalSolarFunctions",
    "SeasonalClearsky", "control_seasonalClearsky", "clearsky_outliers",
    "clearsky_optimizer", "clearsky_delta_optimizer",
    # Seasonal regression
    "SeasonalModel",
    "add_seasonal_features", "add_seasonal_features_dt",
    "seasonalModel_params_to_zeta", "seasonalModel_params_to_phi",
    "seasonalModel_params_to_zeta_jacobian",
    # ARMA
    "ARMAModel",
    # Gaussian mixture
    "GaussianMixtureModel",
    "gm_moments", "gm_moments_match", "gm_loglik", "gm_fit", "gm_fit_moments_match",
    # Radiation model
    "RadiationModel",
    "create_monthly_sequence", "martingale_method_seasonal",
    "reparam_seasonal_function", "integral_sigma_numeric", "integral_sigma2_formula",
    # SolarModel
    "SolarModelSpec", "SolarModel",
]
