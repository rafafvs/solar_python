"""Radiation model class."""

import warnings
import numpy as np
import pandas as pd
import time

from scipy.stats import norm
from scipy.integrate import quad


from .zzz import number_of_day
from .radiationModel_internals import shift, create_monthly_sequence, martingale_method_seasonal, reparam_seasonal_function, integral_sigma_numeric, integral_sigma2_formula 

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# RadiationModel
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------


class RadiationModel:
    """
    Radiation model for solar GHI under a mean-reverting SDE framework.
    Core engine managing the transition from discrete fits to continuous-time stochastic integrals.
    """

    # ------------------------------------------------------------------
    # Public fields
    # ------------------------------------------------------------------
    version = "1.0.0"
    theta = None
    lambda_S = 0.0

    # ------------------------------------------------------------------
    # Constructor
    # ------------------------------------------------------------------

    def __init__(self, model, correction = False):
        self._model = model.clone(deep=True)
        self._measure = "P"
        self._lambda = 0.0
        self._k1 = 1.0
        self._k2 = 1.0
        self._integral_variance = None
        self._integral_expectation = None
        self._seasonal_variance = None

        # 1) Estimate mean-reversion parameter θ
        df_fit = self._model.data
        df_fit = df_fit[df_fit["isTrain"] & (df_fit["weights"] != 0)]
        self.theta = martingale_method_seasonal(df_fit["Yt"], df_fit["Yt_bar"])

        # Update discrete AR(1) coefficient from θ, then re-filter / re-fit
        self._model.update({"phi_1": np.exp(-self.theta)})
        self._model.filter()
        self._model.fit_seasonal_variance()
        self._model.filter()
        self._model.fit_NM_model()
        self._model.update_NM_classification()
        self._model.update_moments()
        self._model.update_logLik()

        # 2) Reparametrise seasonal variance in continuous time
        self.parametrize_seasonal_variance()

        # 3) Optionally correct mixture means for integral discrepancy
        if correction:
            self.correct_NM_coefficients()

    # ------------------------------------------------------------------
    # Public methods
    # ------------------------------------------------------------------

    def parametrize_seasonal_variance(self) -> None:
        """Compute the parameters of the seasonal variance given OLS estimates."""
        omega = 2 * np.pi / 365
        reparam = reparam_seasonal_function(
            self._model.seasonal_variance.coefficients,
            self.theta,
            omega=omega,
        )
        # Clone and update the seasonal variance object
        self._seasonal_variance = self.model.seasonal_variance.clone(deep=True)
        self._seasonal_variance.extra_params["reparam"] = reparam

        # Rename c_ to match the original coefficient names, then update
        c_ = reparam["c_"]
        if hasattr(c_, "index"):                    # pd.Series
            c_.index = self.model.seasonal_variance.coefficients.index
        self._seasonal_variance.update(c_)

        # Store integral callables
        self._integral_variance = integral_sigma2_formula(
            self.theta, reparam["gamma"], omega=omega
        )
        self._integral_expectation = integral_sigma_numeric(
            self.theta, c_, omega=omega
        )

    def correct_NM_coefficients(self) -> None:
        """
        Correct mixture means and std deviations for the discrepancy between
        the integral seasonal std. deviation and variance.
        """
        t = np.arange(1, 366)          # days 1 … 365

        sigma_bar_J = np.sqrt(self._integral_variance(t - 1, t, t))
        sigma_bar_I = self._integral_expectation(t - 1, t, t)
        self._k1 = float(np.mean(sigma_bar_I / sigma_bar_J))

        # Original variance and parameters
        v      = self.model.NM_model.moments["variance"]
        sigma  = self.model.NM_model.sd.copy()
        probs  = self.model.NM_model.p

        # Adjusted means (multiply by k1)
        means = self.model.NM_model.means * self._k1

        # New expectation and k2
        m1 = means[:, 0] * probs[:, 0] + means[:, 1] * probs[:, 1]
        k2 = (
            v + m1 - (means[:, 0] - means[:, 1]) ** 2 * probs[:, 0] * probs[:, 1]
        ) / (sigma[:, 0] ** 2 * probs[:, 0] + sigma[:, 1] ** 2 * probs[:, 1])
        self._k2 = np.sqrt(k2)

        # Update mixture parameters
        self._model.NM_model.update(means=means)

    def change_measure(self, measure: str) -> None:
        """Change the reference probability measure."""
        if measure not in ("P", "Q"):
            raise ValueError("`measure` must be 'P' or 'Q'.")
        self._measure = measure
        self._lambda = self.lambda_S if measure == "Q" else 0.0

    # ------------------------------------------------------------------
    # Point-in-time seasonal functions
    # ------------------------------------------------------------------

    def Ct(self, t_now) -> float:
        """Clear-sky radiation for day t_now."""
        return self._model.seasonal_model_Ct.predict(number_of_day(t_now))

    def Yt_bar(self, t_now) -> float:
        """Seasonal mean of the transformed variable Y_t for day t_now."""
        return self._model.seasonal_model_Yt.predict(number_of_day(t_now))

    def Rt_bar(self, t_now) -> float:
        """Seasonal mean of solar radiation R_t for day t_now."""
        return self._model.Y_to_R(self.Yt_bar(t_now), t_now)

    def sigma_bar(self, t_now) -> float:
        """Instantaneous seasonal std. deviation σ̄_t for day t_now."""
        return float(np.sqrt(self._seasonal_variance.predict(number_of_day(t_now))))

    def mu_B(self, t_now, B: int = 1) -> float:
        """Mixture drift μ_B for component B (1 or 2) at day t_now."""
        if B == 1:
            return self._model.NM_model.mu1.predict(t_now)
        return self._model.NM_model.mu2.predict(t_now)

    def sigma_B(self, t_now, B: int = 1) -> float:
        """
        Mixture diffusion σ_B for component B (1 or 2) at day t_now.
        """
        if B == 1:
            return self._model.NM_model.sd1.predict(t_now)
        return self._model.NM_model.sd2.predict(t_now)

    def mu_Y(self, Yt, t_now, B = 1, dt = 1.0) -> float:
        """Drift μ_Y of the transformed process Y_t."""
        n = number_of_day(t_now)
        dYt_bar = (
            self._model.seasonal_model_Yt.predict(n + dt)
            - self._model.seasonal_model_Yt.predict(n)
        )
        sigma_bar = self.sigma_bar(n)
        drift_Y = dYt_bar - self.theta * (Yt - self.Yt_bar(n)) + sigma_bar * self.mu_B(n, B)
        # !!! Change not already implemented for QT measure
        q_adj = sigma_bar * self.sigma_B(n, B) * self._lambda * (1 if self._measure == "Q" else 0)
        return drift_Y + q_adj

    def sigma_Y(self, t_now, B = 1) -> float:
        """Diffusion σ_Y of the transformed process Y_t."""
        return self.sigma_bar(t_now) * self.sigma_B(t_now, B)

    def mu_R(self, Rt, t_now, B = 1, dt = 1.0) -> float:
        """Drift μ_R of the solar radiation process R_t (Itô drift via chain rule)."""
        Yt  = self._model.R_to_Y(Rt, t_now)
        n   = number_of_day(t_now)
        Ct  = self.Ct(n)
        dCt = self.Ct(n + dt) - Ct

        alpha = self._model.transform.alpha
        beta  = self._model.transform.beta

        Kt  = 1 - alpha - beta * np.exp(-np.exp(Yt))
        mu_Y_val  = self.mu_Y(Yt, t_now, B, dt)
        sig_Y_val = self.sigma_Y(t_now, B)

        return (
            Kt * dCt / dt
            + Ct * beta * np.exp(Yt - np.exp(Yt))
            * (mu_Y_val + 0.5 * (1 - np.exp(Yt)) * sig_Y_val ** 2)
        )

    def sigma_R(self, Rt, t_now, B = 1) -> float:
        """Diffusion σ_R of the solar radiation process R_t."""
        Yt = self._model.R_to_Y(Rt, t_now)
        return (
            self.Ct(t_now)
            * self.model.transform.beta
            * np.exp(Yt - np.exp(Yt))
            * self.sigma_Y(t_now, B)
        )

    # ------------------------------------------------------------------
    # Integral formulas
    # ------------------------------------------------------------------

    def integral_expectation(self, t_now, t_hor, df_date = None, last_day = True) -> pd.DataFrame:
        """Compute the integral ∫ σ̄(s) ds for the expectation."""
        if df_date is None:
            df_date = create_monthly_sequence(t_now, t_hor, last_day=last_day)
        df_date = df_date.copy()
        df_date["int_sigma"] = self._integral_expectation(
            df_date["n"], df_date["N"], df_date["tau"]
        )
        return df_date

    def integral_variance(self, t_now, t_hor, df_date = None, last_day = True) -> pd.DataFrame:
        """Compute the integral ∫ σ̄²(s) ds for the variance."""
        if df_date is None:
            df_date = create_monthly_sequence(t_now, t_hor, last_day=last_day)
        df_date = df_date.copy()
        df_date["int_sigma2"] = self._integral_variance(
            df_date["n"], df_date["N"], df_date["tau"]
        )
        return df_date

    def e_mix_drift(self, t_now = None, t_hor = None, df_date = None) -> dict:
        """Integral mixture drift for both components of Y_t."""
        if df_date is None:
            df_date = create_monthly_sequence(t_now, t_hor, last_day=True)

        df_params = self._model.NM_model.coefficients.copy()
        df_params["e_mu_B"] = df_params["mu1"] * df_params["p1"] + df_params["mu2"] * df_params["p2"]
        df_params["e_sigma_B"] = df_params["sd1"] * df_params["p1"] + df_params["sd2"] * df_params["p2"]

        df = pd.merge(df_date, df_params, on="Month", how="left").sort_values("n")
        df["int"] = self._integral_expectation(df["n"], df["N"], df["tau"])

        body  = df.iloc[:-1]    # rows t … T-1
        last  = df.iloc[-1]     # row T

        # Mean-drift integrals
        sum_mu   = float((body["int"] * body["e_mu_B"]).sum())
        e_mu1    = sum_mu + float(last["int"] * last["mu1"])
        e_mu2    = sum_mu + float(last["int"] * last["mu2"])

        # Diffusion-drift integrals
        sum_sd   = float((body["int"] * body["e_sigma_B"]).sum())
        e_sd1    = sum_sd + float(last["int"] * last["sd1"])
        e_sd2    = sum_sd + float(last["int"] * last["sd2"])

        lam = self._lambda
        return {
            "e_drift_1": e_mu1 + e_sd1 * lam, 
            "e_drift_2": e_mu2 + e_sd2 * lam,
            "e_mu1": e_mu1, "e_mu2": e_mu2, 
            "e_sd1": e_sd1, "e_sd2": e_sd2,
        }

    def e_mix_diffusion(self, t_now = None, t_hor = None, df_date = None) -> dict:
        """Integral mixture diffusion (variance) for both components of Y_t."""
        if df_date is None:
            df_date = create_monthly_sequence(t_now, t_hor, last_day=True)

        df_params = self._model.NM_model.coefficients.copy()
        df_params["mu_diff"]    = df_params["mu1"]  - df_params["mu2"]
        df_params["sigma_diff"] = df_params["sd1"]  - df_params["sd2"]
        df_params["sigma2_1"]   = df_params["sd1"] ** 2
        df_params["sigma2_2"]   = df_params["sd2"] ** 2

        df = pd.merge(df_date, df_params, on="Month", how="left").sort_values("n")
        df["int_sigma2"] = self._integral_variance(df["n"], df["N"], df["tau"])

        df_t = df.iloc[:-1]   # t … T-1
        df_T = df.iloc[-1]    # T (last row)

        variance_drift_P  = float((df_t["int_sigma2"] * df_t["mu_diff"]    ** 2 * df_t["p1"] * df_t["p2"]).sum())
        variance_drift_Q  = float((df_t["int_sigma2"] * df_t["sigma_diff"] ** 2 * df_t["p1"] * df_t["p2"]).sum())
        variance_diffusion = float((df_t["int_sigma2"] * (df_t["sigma2_1"] * df_t["p1"] + df_t["sigma2_2"] * df_t["p2"])).sum())

        common_variance = variance_drift_P + variance_diffusion + variance_drift_Q * self._lambda ** 2

        last_1 = float(df_T["int_sigma2"] * df_T["sigma2_1"])
        last_2 = float(df_T["int_sigma2"] * df_T["sigma2_2"])

        return {
            "variance_1":        common_variance + last_1,
            "variance_2":        common_variance + last_2,
            "variance_drift_P":  variance_drift_P,
            "variance_drift_Q":  variance_drift_Q,
            "common_variance":   common_variance + float(df_T["int_sigma2"]),
            "variance_diffusion": variance_diffusion,
            "last_1":            last_1,
            "last_2":            last_2,
        }

    # ------------------------------------------------------------------
    # Conditional moments
    # ------------------------------------------------------------------

    def M_Y(self, Rt, t_now, t_hor, df_date = None) -> dict:
        """Conditional expectation of Y_T given Y_t = R_to_Y(Rt)."""
        if df_date is None:
            df_date = create_monthly_sequence(t_now, t_hor, last_day=True)

        Yt  = self.model.R_to_Y(Rt, t_now)
        tau = float(df_date["tau"].max() - df_date["n"].min())
        mix_drift  = self.e_mix_drift(df_date=df_date)

        Y_forecast = (
            self.Yt_bar(t_hor)
            + (Yt - self.Yt_bar(t_now)) * np.exp(-self.theta * tau)
        )
        return {
            "up":   Y_forecast + mix_drift["e_drift_1"],
            "dw":   Y_forecast + mix_drift["e_drift_2"],
            "e_Yt": Y_forecast,
        }

    def S_Y(self, t_now, t_hor, df_date = None) -> dict:
        """Conditional variance of Y_T."""
        if df_date is None:
            df_date = create_monthly_sequence(t_now, t_hor, last_day=True)

        mix_diff = self.e_mix_diffusion(df_date=df_date)
        return {"up": mix_diff["variance_1"], "dw": mix_diff["variance_2"]}

    def Moments(self, t_now, t_hor, quiet = False) -> pd.DataFrame:
        """Compute the conditional moments table over historical dates."""
        t_now = pd.to_datetime(t_now)
        t_hor = pd.to_datetime(t_hor)
        tau   = int((t_hor.iloc[0] - t_now.iloc[0]).days) if hasattr(t_hor, "iloc") else int((t_hor - t_now).days)

        data = self._model.data.copy()
        data["K_Y"]      = self._model.R_to_Y(data["GHI_bar"], data["date"])
        data["payoff"]   = (data["GHI_bar"] - data["GHI"]) * ((data["GHI_bar"] > data["GHI"]).astype(float))
        data["L1_Yt_bar"] = shift(data["Yt_bar"], tau)
        data["L1_Yt"]    = shift(data["Yt"], tau)
        data = data[(data["date"] >= t_hor.min()) & (data["date"] <= t_hor.max())]

        # 1) Y forecast
        nstep = 1
        t0 = time.time()
        if not quiet:
            print(f"Fitting Y forecast ({nstep}/3) ...", end="", flush=True)
        data["Y_forecast"] = (
            data["Yt_bar"]
            + (data["L1_Yt"] - data["L1_Yt_bar"]) * np.exp(-self.theta * tau)
        )
        if not quiet:
            print(f"Done in {time.time()-t0:.2f} secs", end="\r", flush=True)
        nstep += 1

        # 2) Mixture drifts
        t0 = time.time()
        if not quiet:
            print(f"Fitting Mixture drift ({nstep}/3) ...", end="", flush=True)
        drifts = [self.e_mix_drift(tn, th) for tn, th in zip(t_now, t_hor)]
        data["e_mu1"] = [d["e_mu1"] for d in drifts]
        data["e_mu2"] = [d["e_mu2"] for d in drifts]
        data["e_sd1"] = [d["e_sd1"] for d in drifts]
        data["e_sd2"] = [d["e_sd2"] for d in drifts]
        if not quiet:
            print(f"Done in {time.time()-t0:.2f} secs", end="\r", flush=True)
        nstep += 1

        # 3) Mixture diffusions
        t0 = time.time()
        if not quiet:
            print(f"Fitting Mixture diffusions ({nstep}/3) ...", end="", flush=True)
        diffusions = [self.e_mix_diffusion(tn, th)for tn, th in zip(t_now, t_hor)]
        data["variance_drift_P"]  = [d["variance_drift_P"]  for d in diffusions]
        data["variance_drift_Q"]  = [d["variance_drift_Q"]  for d in diffusions]
        data["variance_diffusion"] = [d["variance_diffusion"] for d in diffusions]
        data["last_1"] = [d["last_1"] for d in diffusions]
        data["last_2"] = [d["last_2"] for d in diffusions]
        if not quiet:
            print(f"Done in {time.time()-t0:.2f} secs", end="\r", flush=True)

        data["alpha"]     = self.model.transform.alpha
        data["beta"]      = self.model.transform.beta
        data["date_from"] = data["date"] - pd.Timedelta(days=tau)

        data.rename(columns={"date_from": "t_now", "date": "t_hor"}, inplace=True)

        return data[
            "t_now", "t_hor", "Y_forecast", "e_mu1", "e_mu2", "e_sd1", "e_sd2", 
            "variance_drift_P", "variance_drift_Q", "variance_diffusion", "last_1", 
            "last_2", "p1", "alpha", "beta", "payoff", "Ct", "Yt_bar", "GHI", "GHI_bar", "K_Y"
        ]

    # ------------------------------------------------------------------
    # Conditional distributions
    # ------------------------------------------------------------------

    def pdf_Y(self, Rt, t_now, t_hor, B = None):
        """Conditional density of Y_T | Y_t."""
        # Create once the sequence of dates
        df_date = create_monthly_sequence(t_now, t_hor, last_day=True)
        M       = self.M_Y(Rt, t_now, t_hor, df_date=df_date)
        S       = self.S_Y(t_now, t_hor, df_date=df_date)
        p       = self._model.NM_model.prob.predict(t_hor)

        mu1, mu2 = M["up"], M["dw"]
        sd1, sd2 = np.sqrt(S["up"]), np.sqrt(S["dw"])

        if B is None:
            return lambda x: p * norm.pdf(x, mu1, sd1) + (1 - p) * norm.pdf(x, mu2, sd2)
        elif B == 1:
            return lambda x: norm.pdf(x, mu1, sd1)
        else:
            return lambda x: norm.pdf(x, mu2, sd2)

    def cdf_Y(self, Rt, t_now, t_hor, B = None):
        """Conditional CDF of Y_T | Y_t."""
        df_date = create_monthly_sequence(t_now, t_hor, last_day=True)
        M       = self.M_Y(Rt, t_now, t_hor, df_date=df_date)
        S       = self.S_Y(t_now, t_hor, df_date=df_date)
        p       = self._model.NM_model.prob.predict(t_hor)

        mu1, mu2 = M["up"], M["dw"]
        sd1, sd2 = np.sqrt(S["up"]), np.sqrt(S["dw"])

        if B is None:
            return lambda x: p * norm.cdf(x, mu1, sd1) + (1 - p) * norm.cdf(x, mu2, sd2)
        elif B == 1:
            return lambda x: norm.cdf(x, mu1, sd1)
        else:
            return lambda x: norm.cdf(x, mu2, sd2)

    def pdf_R(self, Rt, t_now, t_hor, B = None):
        """Conditional density of R_T | R_t."""
        C_T   = self.Ct(t_hor)
        pdf_y = self.pdf_Y(Rt, t_now, t_hor, B)
        alpha = self._model.transform.alpha
        beta  = self._model.transform.beta

        def density(x):
            x   = np.asarray(x, dtype=float)
            z_x = (1 - x / C_T - alpha) / beta
            with np.errstate(invalid="ignore", divide="ignore"):
                u_x = np.log(-np.log(z_x))
            num = pdf_y(u_x)
            den = C_T * beta * np.log(z_x ** z_x)
            probs = -num / den
            probs = np.where(np.isnan(probs), 0.0, probs)
            return probs

        return density

    def cdf_R(self, Rt, t_now, t_hor, B = None):
        """Conditional CDF of R_T | R_t."""
        C_T   = self.Ct(t_hor)
        cdf_y = self.cdf_Y(Rt, t_now, t_hor, B)
        alpha = self._model.transform.alpha
        beta  = self._model.transform.beta

        def distribution(x):
            x   = np.asarray(x, dtype=float)
            z_x = (1 - x / C_T - alpha) / beta
            with np.errstate(invalid="ignore", divide="ignore"):
                u_x = np.log(-np.log(z_x))
            probs = cdf_y(u_x)
            probs = np.where(np.isnan(probs), 0.0, probs)
            return probs

        return distribution

    def e_GHI(self, Rt, t_now, t_hor, B = None, moment = 1) -> float:
        """Conditional moment E[R_T^moment | R_t]."""
        pdf_y = self.pdf_Y(Rt, t_now, t_hor, B)
        C_T   = self.Ct(t_hor)
        alpha = self._model.transform.alpha
        beta  = self._model.transform.beta

        def integrand(y):
            return (C_T * (1 - alpha - beta * np.exp(-np.exp(y)))) ** moment * pdf_y(y)

        result, _ = quad(integrand, -np.inf, np.inf)
        return result

    def v_GHI(self, Rt, t_now, t_hor, B = None) -> float:
        """Conditional variance Var[R_T | R_t] = E[R_T²] - E[R_T]²."""
        return self.e_GHI(Rt, t_now, t_hor, B, 2) - self.e_GHI(Rt, t_now, t_hor, B, 1) ** 2

    def __str__(self):
        """Print method for the class `radiationModel`."""
        msg = f"MRP: {self._lambda}\n"
        msg += f"Measure: {self._measure}\n"
        msg += f"Version: {self.version}\n"
        return msg

    def __repr__(self):
        return self.__str__()

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def model(self):
        """The underlying SolarModel object."""
        return self._model

    @property
    def measure(self) -> str:
        """Current probability measure ('P' or 'Q')."""
        return self._measure

    @property
    def lambda_(self) -> float:
        """Market risk premium in use (0 under P, lambda_S under Q)."""
        return self._lambda