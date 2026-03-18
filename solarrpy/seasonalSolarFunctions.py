import numpy as np
import pandas as pd
import datetime

# Helper function to replicate R's implied `number_of_day` function
# which is used in the R code to convert dates/inputs to day of the year (1-365)
def number_of_day(x):
    """
    Simulates the external R function `number_of_day`.
    Improved to handle Series that are already numeric.
    """
    # 1. Handle Pandas objects
    if isinstance(x, (pd.DatetimeIndex, pd.Series)):
        # Check if the Series/Index is already numeric (int or float)
        if pd.api.types.is_numeric_dtype(x):
            return x
        # If it's datetimelike, extract the day of the year
        return x.dt.dayofyear
    
    # 2. Handle single Python datetime objects
    if isinstance(x, (pd.Timestamp, datetime.datetime, datetime.date)):
        return x.timetuple().tm_yday
    
    # 3. Fallback for scalars (ints, floats)
    return x


class SeasonalSolarFunctions:
    """
    Solar seasonal functions
    Version 1.0.0
    References: Duffie, Solar Engineering of Thermal Processes Fourth Edition.
    """

    def __init__(self, method="spencer", legal_hour=True):
        """
        Initialize a `seasonalSolarFunctions` object
        :param method: str, method type for computations. Can be `cooper` or `spencer`.
        :param legal_hour: bool, when `True` the clock time will be corrected for the legal hour.
        """
        valid_methods = ["spencer", "cooper"]
        if method not in valid_methods:
            raise ValueError(f"method must be one of {valid_methods}")
        
        self.__method_ = method
        self.legal_hour = legal_hour
        
        # Private fields equivalent
        self.__version = "1.0.0"
        self.__Gsc = 1367

    @property
    def Gsc(self):
        """solar constant, i.e., `1367`."""
        return self.__Gsc

    def update_method(self, x=None):
        """
        Extract or update the method used for computations.
        """
        valid_methods = ["spencer", "cooper"]
        if x is None:
            return self.__method_
        else:
            if x not in valid_methods:
                raise ValueError(f"method must be one of {valid_methods}")
            self.__method_ = x
            return self.__method_

    def B(self, n):
        """
        Seasonal adjustment parameter.
        Eq. 1.4.2 from Duffie (4th edition)
        """
        return (2 * np.pi * n) / 365

    def degree(self, x):
        """
        Convert angles in radiant into an angles in degrees.
        """
        return x * 180 / np.pi

    def radiant(self, x):
        """
        Convert angles in degrees into an angles in radiant.
        """
        return x * np.pi / 180

    def E(self, n):
        """
        Compute the time adjustment in minutes.
        Eq. 1.5.3 from Duffie (4th edition)
        Returns: The time adjustment as timedelta (minutes)
        """
        n_day = number_of_day(n)
        B = self.B(n_day - 1)
        # Time adjustment in minutes
        E_val = 229.2 * (0.000075 + 0.001868 * np.cos(B) - 0.032077 * np.sin(B) - 
                         0.014615 * np.cos(2 * B) - 0.04089 * np.sin(2 * B))
        # Time adjustment in seconds equivalent (pandas Timedelta)
        return pd.to_timedelta(E_val, unit='m')

    def elevation(self, alt):
        """
        Compute the angle in the degree given a certain altitude in meters.
        """
        Re = 6371 * 1000
        alpha = np.arccos(Re / (Re + alt))
        return self.degree(alpha)

    def solar_time(self, x, lon, lon_st=15, tz="Europe/Rome"):
        """
        Compute the solar time from a clock time.
        Eq. 1.5.2 from Duffie (4th edition)
        """
        # Convert to datetime (pandas handles vectors and scalars)
        if not isinstance(x, (pd.DatetimeIndex, pd.Series)):
            date_h = pd.to_datetime(x).tz_localize(tz) if pd.to_datetime(x).tz is None else pd.to_datetime(x).tz_convert(tz)
        else:
            date_h = x
            
        # The R code has an `if (FALSE)` block here which is dead code, but we keep the structure
        if False:
            # Replicating the dead code logic exactly as in R
            year_val = date_h[0].year if isinstance(date_h, pd.DatetimeIndex) else date_h.year
            start_legal_hour = pd.to_datetime(f"{year_val}-03-10")
            end_legal_hour = pd.to_datetime(f"{year_val}-10-30")
            is_legal_hour = (pd.to_datetime(date_h.date) >= start_legal_hour) & (pd.to_datetime(date_h.date) <= end_legal_hour)
            date_h = np.where(is_legal_hour, date_h - pd.Timedelta(hours=1), date_h)
            date_h = pd.DatetimeIndex(date_h)

        # Solar hour for the selected day-time
        return date_h + pd.to_timedelta(4 * (lon - lon_st), unit='m') + self.E(date_h)

    def solar_hour(self, LST):
        """
        Compute the solar hour for a specific clock time.
        """
        LST = pd.to_datetime(LST)
        if isinstance(LST, (pd.DatetimeIndex, pd.Series)):
            return LST.dt.hour + LST.dt.minute / 60 + LST.dt.second / 3600
        else:
            return LST.hour + LST.minute / 60 + LST.second / 3600

    def hour_angle(self, LST):
        """
        Compute the solar angle for a specific hour of the day.
        Eq. 1.42 from Comini (2013)
        """
        solar_hour = self.solar_hour(LST)
        return 15 * (solar_hour - 12)

    def incidence_angle(self, LST, lat, alt=0, beta=0, gamma=0):
        """
        Compute the incidence angle
        """
        beta = self.radiant(beta)
        gamma = self.radiant(gamma)
        phi = self.radiant(lat)
        omega = self.radiant(self.hour_angle(LST))
        delta = self.radiant(self.declination(LST))

        T_ = np.sin(delta) * (np.sin(phi) * np.cos(beta) - np.cos(phi) * np.sin(beta) * np.cos(gamma))
        U_ = np.cos(delta) * (np.cos(phi) * np.cos(beta) + np.sin(phi) * np.sin(beta) * np.cos(gamma))
        V_ = np.cos(delta) * np.sin(beta) * np.sin(gamma)

        cos_theta_z = T_ + U_ * np.cos(omega) + V_ * np.sin(omega)
        
        # Equivalent to R's cos_theta_z[cos_theta_z < 0] <- 0
        cos_theta_z = np.maximum(cos_theta_z, 0)
        
        return cos_theta_z

    def azimut_angle(self, LST, lat, alt, beta=0, gamma=0):
        """
        Compute the solar azimuth angle for a specific time of the day.
        Eq. 1.6.6 from Duffie (4th edition)
        """
        phi = self.radiant(lat)
        delta = self.radiant(self.declination(LST))
        omega = self.radiant(self.hour_angle(LST))
        theta_z = self.radiant(self.incidence_angle(LST, lat, alt, beta=0, gamma=0))

        numer = (np.cos(theta_z) * np.sin(phi) - np.sin(delta))
        denom = np.sin(theta_z) * np.cos(phi)
        
        # Using numpy clip instead of pmin/pmax
        ratio = np.clip(numer / denom, -1, 1)
        
        # NOTE: R code references `angle` which is undefined, meaning it likely meant `ratio`. 
        # To strictly replicate while allowing execution, we set angle = ratio.
        angle = ratio 
        gamma_s = np.sign(omega) * np.abs(np.arccos(angle))
        
        return self.degree(gamma_s)

    def Gon(self, n, deriv=False):
        """
        Compute the solar constant adjusted for the day of the year.
        Eq. 1.4.1a or 1.4.1b from Duffie (4th edition)
        """
        n_day = number_of_day(n)
        Gsc = self.__Gsc

        if self.__method_ == "cooper":
            B = self.B(n_day)
            if not deriv:
                Gn = Gsc * (1 + 0.033 * np.cos(B))
            else:
                Gn = -Gsc * (2 * np.pi / 365) * 0.033 * np.sin(B)
        elif self.__method_ == "spencer":
            B = self.B(n_day - 1)
            if not deriv:
                Gn = Gsc * (1.000110 + 0.034221 * np.cos(B) + 0.001280 * np.sin(B) + 
                            0.000719 * np.cos(2 * B) + 0.000077 * np.sin(2 * B))
            else:
                Gn = Gsc * (2 * np.pi / 365) * (-0.034221 * np.sin(B) + 0.001280 * np.cos(B) - 
                                                0.001438 * np.sin(2 * B) + 0.000154 * np.cos(2 * B))
        return Gn

    def declination(self, n, deriv=False):
        """
        Compute solar declination in degrees.
        Eq. 1.6.1a or 1.6.1b from Duffie (4th edition)
        """
        n_day = number_of_day(n)
        if self.__method_ == "cooper":
            if not deriv:
                declination = 23.45 * np.sin(2 * np.pi * (284 + n_day) / 365)
            else:
                declination = 23.45 * np.cos(2 * np.pi * (284 + n_day) / 365) * (2 * np.pi / 365)
        elif self.__method_ == "spencer":
            B = self.B(n_day - 1)
            if not deriv:
                declination = (180 / np.pi) * (0.006918 - 0.399912 * np.cos(B) + 0.070257 * np.sin(B) - 
                                               0.006758 * np.cos(2 * B) + 0.000907 * np.sin(2 * B) - 
                                               0.002697 * np.cos(3 * B) + 0.00148 * np.sin(3 * B))
            else:
                declination = (360 / 365) * (0.399912 * np.sin(B) + 0.070257 * np.cos(B) + 
                                             0.013516 * np.sin(2 * B) + 0.001814 * np.cos(2 * B) + 
                                             0.008091 * np.sin(3 * B) + 0.00444 * np.cos(3 * B))
        return declination

    def Hon(self, n, lat, alt, deriv=False):
        """
        Compute the solar extraterrestrial radiation
        Eq. 1.10.3 from Duffie (4th edition)
        """
        phi = self.radiant(lat)
        delta = self.radiant(self.declination(n))
        omega_s = self.radiant(self.sunset_hour_angle(n, lat, alt))
        A = (24 * 3600) / (np.pi * 3600000)
        Gon = self.Gon(n)
        B_n = (np.cos(phi) * np.cos(delta) * np.sin(omega_s) + omega_s * np.sin(phi) * np.sin(delta))
        
        if not deriv:
            Hon = A * Gon * B_n
        else:
            d_delta = self.radiant(self.declination(n, deriv=True))
            d_omega_s = self.radiant(self.sunset_hour_angle(n, lat, alt, deriv=True))
            d_B_n_d_delta = -np.sin(delta) * np.cos(phi) * np.sin(omega_s) + omega_s * np.cos(delta) * np.sin(phi)
            d_B_n_d_omega = np.cos(delta) * np.cos(phi) * np.cos(omega_s) + np.sin(delta) * np.sin(phi)
            d_B_n_d_n = d_B_n_d_delta * d_delta + d_B_n_d_omega * d_omega_s
            d_Gon_dn = self.Gon(n, deriv=True)
            Hon = A * (B_n * d_Gon_dn + Gon * d_B_n_d_n)
            
        return Hon

    def sunset_hour_angle(self, n, lat, alt=None, deriv=False):
        """
        Compute solar angle at sunset in degrees
        Eq. 1.6.10 from Duffie (4th edition)
        """
        phi = self.radiant(lat)
        declination = self.radiant(self.declination(n))

        if alt is None:
            arg = -np.tan(declination) * np.tan(phi)
            if not deriv:
                arg = np.clip(arg, -1, 1)
                omega_s = self.degree(np.arccos(arg))
            else:
                d_declination = self.declination(n, deriv=True)
                d_arg = -np.tan(phi) / np.cos(declination)**2
                omega_s = -(1 / np.sqrt(1 - arg**2)) * d_arg * d_declination
            return omega_s

        alt_rad = self.radiant(self.elevation(alt))
        arg = (np.sin(alt_rad) - np.sin(declination) * np.sin(phi)) / (np.cos(phi) * np.cos(declination))
        
        if not deriv:
            arg = np.clip(arg, -1, 1)
            omega_s = self.degree(np.arccos(arg))
        else:
            d_declination = self.declination(n, deriv=True)
            d_arg = (np.sin(alt_rad) * np.sin(declination) - np.sin(phi)) / (np.cos(phi) * np.cos(declination)**2)
            omega_s = -(1 / np.sqrt(1 - arg**2)) * d_arg * d_declination
            
        return omega_s

    def sun_hours(self, n, lat, alt):
        """
        Compute number of sun hours for a day n.
        Eq. 1.6.11 from Duffie (4th edition)
        """
        sun_hours_val = self.sunset_hour_angle(n, lat, alt) * (2 / 15)
        return pd.to_timedelta(sun_hours_val, unit='h')

    def solar_altitude(self, n, lat):
        """
        Compute solar altitude in degrees
        """
        phi = self.radiant(lat)
        declination = self.radiant(self.declination(n))
        arg = np.sin(declination) * np.sin(phi) + np.cos(declination) * np.cos(phi)
        arg = np.clip(arg, -1, 1)
        return self.degree(np.arcsin(arg))

    def solar_angles(self, x, lat, lon, alt, lon_st=15, beta=0, gamma=0, by="1min", tz="Europe/Rome"):
        """
        Compute the solar angle for a latitude in different dates.
        """
        if isinstance(x, (pd.DatetimeIndex, pd.Series, list)):
            x = x[0]
            
        day_date = pd.to_datetime(x).date()
        start_date = pd.to_datetime(f"{day_date} 00:00:00").tz_localize(tz)
        end_date = pd.to_datetime(f"{day_date + datetime.timedelta(days=1)} 00:00:00").tz_localize(tz)
        
        # Format pandas freq from R's `1 min` to `1min` string equivalent
        freq = str(by).replace(" ", "")
        day_date_seq = pd.date_range(start=start_date, end=end_date, freq=freq)
        
        phi = self.radiant(lat)
        declination = self.declination(pd.to_datetime(day_date))
        omega_max = self.sunset_hour_angle(pd.to_datetime(day_date), lat, alt)
        omega_min = -omega_max
        sun_hours_val = self.sun_hours(pd.to_datetime(day_date), lat, alt)
        E_val = self.E(pd.to_datetime(day_date))
        
        LST = self.solar_time(day_date_seq, lon, lon_st, tz=tz)
        omega = self.hour_angle(LST)
        theta = self.incidence_angle(LST, lat, alt, beta, gamma)
        
        # NOTE: Original R code binds "solartime = solartime" which throws an error since LST holds the calculation.
        # Adjusted to map to LST here to prevent breaking the code execution.
        output = pd.DataFrame({
            'date': day_date,
            'clocktime': day_date_seq,
            'solartime': LST, 
            'lat': lat,
            'lon': lon,
            'omega': omega,
            'declination': declination,
            'omega_min': omega_min,
            'omega_max': omega_max,
            'sun_hours': sun_hours_val,
            'theta': theta,
            'E': E_val
        })
        return output

    def clearsky(self, cosZ=None, G0=None, alt=None, clime="No Correction"):
        """
        Hottel clearsky
        """
        valid_climes = ["No Correction", "Summer", "Winter", "Subartic Summer", "Tropical"]
        
        # Handle match.arg logic from R (uses first element if vector/list is passed)
        if isinstance(clime, (list, tuple, np.ndarray, pd.Series)):
            clime = clime[0]
        if clime not in valid_climes:
            raise ValueError(f"clime must be one of {valid_climes}")

        altitude = 0 if alt is None else alt
        
        if altitude > 2.5:
            a0_star = 0.6 * (1 - np.exp(-0.214 * (altitude - 1.12)))
        else:
            a0_star = 0.4237 - 0.00821 * (6.0 - altitude)**2
            
        a1_star = 0.5055 - 0.00595 * (6.5 - altitude)**2
        a2_star = 0.2711 - 0.01858 * (2.5 - altitude)**2

        a = np.array([a0_star, a1_star, a2_star])
        
        if clime == "Summer":
            correction_factor = np.array([0.97, 0.99, 1.02])
            a = a * correction_factor
        elif clime == "Winter":
            correction_factor = np.array([1.03, 1.01, 1.00])
            a = a * correction_factor
        elif clime == "Subartic Summer":
            correction_factor = np.array([0.99, 0.99, 1.01])
            a = a * correction_factor
        elif clime == "Tropical":
            correction_factor = np.array([0.95, 0.98, 1.02])
            a = a * correction_factor

        # Note: 0-based indexing in Python vs 1-based in R
        a0 = a[0]
        a1 = a[1]
        a2 = a[2]

        tau_beam = a0 + a1 * np.exp(-a2 / cosZ)
        tau_diffuse = 0.271 - 0.294 * tau_beam
        
        output = pd.DataFrame({'tau_beam': tau_beam, 'tau_diffuse': tau_diffuse})
        
        output['tau_beam'] = np.where((output['tau_beam'] > 1) | (output['tau_beam'] < 0), 0, output['tau_beam'])
        output['tau_diffuse'] = np.where((output['tau_diffuse'] > 1) | (output['tau_diffuse'] < 0), 0, output['tau_diffuse'])
        
        skymax = G0 * output['tau_beam'] + G0 * output['tau_diffuse']
        return skymax