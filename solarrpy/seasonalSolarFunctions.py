import numpy as np
import pandas as pd
from .zzz import number_of_day

class SeasonalSolarFunctions:

    def __init__(self, method='spencer', legal_hour=True):  
        # Method of computation
        if method not in ["spencer", "cooper"]:
            raise ValueError("Method must be 'spencer' or 'cooper'")
        
        self.version = '1.0.0'
        self._Gsc = 1367
        self.method_ = method
        self.legal_hour = legal_hour

    
    @property
    def Gsc(self):
        """Solar constant, i.e., 1367."""
        return self._Gsc
        
    @property
    def method(self):
        return self._method
        
    @method.setter
    def method(self, val):
        if val not in ["spencer", "cooper"]:
            raise ValueError("Method must be 'spencer' or 'cooper'")
        self._method = val

    # =================================================================
    # Base Math & Angles
    # =================================================================

    def B(self, n):
        """Seasonal adjustment parameter (Eq. 1.4.2)."""
        return (2 * np.pi * n) / 365

    def degree(self, x):
        """Convert radians to degrees."""
        return np.degrees(x)

    def radiant(self, x):
        """Convert degrees to radians."""
        return np.radians(x)

    def E(self, n):
        """Compute the time adjustment in minutes."""
        n = number_of_day(n)
            
        B = self.B(n - 1)
        # Time adjustment in minutes (Eq 1.5.3)
        E_min = 229.2 * (0.000075 + 0.001868 * np.cos(B) - 0.032077 * np.sin(B) 
                         - 0.014615 * np.cos(2 * B) - 0.04089 * np.sin(2 * B))
        
        # Convert to Pandas Timedelta
        return pd.to_timedelta(E_min, unit='m')

    def elevation(self, alt):
        """Compute the angle in degrees given altitude in meters."""
        Re = 6371 * 1000
        alpha = np.arccos(Re / (Re + alt))
        return self.degree(alpha)

    # =================================================================
    # Solar Time & Angles
    # =================================================================

    def solar_time(self, x, lon, lon_st=15, tz="Europe/Rome"):
        """Compute the solar time from a clock time."""
        if not isinstance(x, (pd.DatetimeIndex, pd.Series)):
            date_h = pd.to_datetime(x).tz_localize(tz) if pd.to_datetime(x).tz is None else pd.to_datetime(x).tz_convert(tz)
        else:
            date_h = x

        # Solar hour for the selected day-time
        return date_h + pd.to_timedelta(4 * (lon - lon_st), unit='m') + self.E(date_h)

    def solar_hour(self, LST):
        """Compute the decimal solar hour for a specific clock time."""
        LST = pd.to_datetime(LST)
        return LST.hour + (LST.minute / 60) + (LST.second / 3600)

    def hour_angle(self, LST):
        """Compute the solar angle in degrees for a specific hour of the day."""
        solar_hr = self.solar_hour(LST)
        return 15 * (solar_hr - 12)

    def incidence_angle(self, LST, lat, alt=0, beta=0, gamma=0):
        """Compute the incidence angle for a specific time of the day."""
        # Altitude of the surface in radiant
        beta = self.radiant(beta)
        # Orientation of the surface in radiant
        gamma = self.radiant(gamma)
        # Latitude in radiant
        phi = self.radiant(lat)
        # Solar hour angle
        omega = self.radiant(self.hour_angle(LST))
        # Declination in radiant
        delta = self.radiant(self.declination(LST))
        # Components
        T_ = np.sin(delta) * (np.sin(phi) * np.cos(beta) - np.cos(phi) * np.sin(beta) * np.cos(gamma))
        U_ = np.cos(delta) * (np.cos(phi) * np.cos(beta) + np.sin(phi) * np.sin(beta) * np.cos(gamma))
        V_ = np.cos(delta) * np.sin(beta) * np.sin(gamma)
        # Cosine of the angle of incidence
        cos_theta_z = T_ + U_ * np.cos(omega) + V_ * np.sin(omega)
        
        return np.maximum(cos_theta_z, 0)

    def azimut_angle(self, LST, lat, alt=0, beta=0, gamma=0):
        """Compute the azimut angle for a specific time of the day."""
        # Latitude in radiant
        phi = self.radiant(lat)
        # Declination in radiant
        delta = self.radiant(self.declination(LST))
        # Solar hour angle
        omega = self.radiant(self.hour_angle(LST))
        # The angle of incidence is the zenith angle of the sun
        theta_z = self.radiant(self.incidence_angle(LST, lat, alt, beta=0, gamma=0))
        # Azimut angle
        numer = (np.cos(theta_z) * np.sin(phi) - np.sin(delta))
        denom = np.sin(theta_z) * np.cos(phi)
        ratio = np.clip(numer / denom, -1, 1)
        gamma_s = np.sign(omega) * np.abs(np.arccos(ratio))
        
        return self.degree(gamma_s)

    # =================================================================
    # Solar Constant & Declination (with Calculus)
    # =================================================================

    def Gon(self, n, deriv=False):
        """Compute the adjusted solar constant (with optional derivative)."""
        n = number_of_day(n)
        
        # Convert Gsc to float if it is a list, tuple, series, or array
        Gsc = self._Gsc[0] if isinstance(self._Gsc, (list, tuple, pd.Series, np.ndarray)) else self._Gsc
        Gsc = float(Gsc)
        
        if self.method_ == "cooper":
            B = self.B(n)
            if not deriv:
                return Gsc * (1 + 0.033 * np.cos(B))
            else:
                return -Gsc * (2 * np.pi / 365) * 0.033 * np.sin(B)
                
        elif self.method_ == "spencer":
            B = self.B(n - 1)
            if not deriv:
                return Gsc * (1.000110 + 0.034221 * np.cos(B) + 0.001280 * np.sin(B) 
                              + 0.000719 * np.cos(2 * B) + 0.000077 * np.sin(2 * B))
            else:
                return Gsc * (2 * np.pi / 365) * (-0.034221 * np.sin(B) + 0.001280 * np.cos(B) 
                                                  - 0.001438 * np.sin(2 * B) + 0.000154 * np.cos(2 * B))


    def declination(self, n, deriv=False):
        """Compute solar declination in degrees (with optional derivative)."""
        n = number_of_day(n)
        
        if self.method_ == "cooper":
            if not deriv:
                return 23.45 * np.sin(2 * np.pi * (284 + n) / 365)
            else:
                return 23.45 * np.cos(2 * np.pi * (284 + n) / 365) * (2 * np.pi / 365)
                
        elif self.method_ == "spencer":
            B = self.B(n - 1)
            if not deriv:
                return (180 / np.pi) * (0.006918 - 0.399912 * np.cos(B) + 0.070257 * np.sin(B) 
                                        - 0.006758 * np.cos(2 * B) + 0.000907 * np.sin(2 * B) 
                                        - 0.002697 * np.cos(3 * B) + 0.00148 * np.sin(3 * B))
            else:
                return (360 / 365) * (0.399912 * np.sin(B) + 0.070257 * np.cos(B) 
                                      + 0.013516 * np.sin(2 * B) + 0.001814 * np.cos(2 * B) 
                                      + 0.008091 * np.sin(3 * B) + 0.00444 * np.cos(3 * B))

    # =================================================================
    # Advanced Geometry & Extraterrestrial Radiation
    # =================================================================

    def sunset_hour_angle(self, n, lat, alt, deriv=False):
        """Compute solar angle at sunset in degrees."""
        phi = self.radiant(lat)
        declination = self.radiant(self.declination(n))
        
        if alt is None:
            arg = -np.tan(declination) * np.tan(phi)
            if not deriv:
                arg = np.clip(arg, -1, 1)
                return self.degree(np.arccos(arg))
            else:
                d_declination = self.declination(n, deriv=True)
                d_arg = -np.tan(phi) / (np.cos(declination)**2)
                return -(1 / np.sqrt(1 - arg**2)) * d_arg * d_declination
                
        # Generalized version with altitude
        alt_rad = self.radiant(self.elevation(alt))
        arg = (np.sin(alt_rad) - np.sin(declination) * np.sin(phi)) / (np.cos(phi) * np.cos(declination))
        
        if not deriv:
            arg = np.clip(arg, -1, 1)
            return self.degree(np.arccos(arg))
        else:
            d_declination = self.declination(n, deriv=True)
            d_arg = (np.sin(alt_rad) * np.sin(declination) - np.sin(phi)) / (np.cos(phi) * (np.cos(declination)**2))
            return -(1 / np.sqrt(1 - arg**2)) * d_arg * d_declination

    def Hon(self, n, lat, alt=None, deriv=False):
        """Compute the extraterrestrial radiation."""
        # Latitude in radiant
        phi = self.radiant(lat)
        # Declination in radiant
        delta = self.radiant(self.declination(n))
        # Sunset hour angle in radiant
        omega_s = self.radiant(self.sunset_hour_angle(n, lat, alt))
        # Scale factor -> (kilowatt hour)/m^2 per day
        A = (24 * 3600) / (np.pi * 3600000)
        # Extraterrestrial radiation in daily joules/m^2 per day
        Gon = self.Gon(n)
        B_n = (np.cos(phi) * np.cos(delta) * np.sin(omega_s) + omega_s * np.sin(phi) * np.sin(delta))
        
        if not deriv:
            return A * Gon * B_n
        else:
            # Derivative of the declination in radiant
            d_delta = self.radiant(self.declination(n, deriv=True))
            # Derivative of the sunset hour angle in radiant
            d_omega_s = self.radiant(self.sunset_hour_angle(n, lat, alt, deriv=True))
            # Derivative of the intern angle with respect to declination (delta)
            d_B_n_d_delta = -np.sin(delta) * np.cos(phi) * np.sin(omega_s) + omega_s * np.cos(delta) * np.sin(phi)
            # Derivative of the intern angle with respect to sunset hour angle (omega_s)
            d_B_n_d_omega = np.cos(delta) * np.cos(phi) * np.cos(omega_s) + np.sin(delta) * np.sin(phi)
            # Derivative of the intern angle with respect to time (n)
            d_B_n_d_n = d_B_n_d_delta * d_delta + d_B_n_d_omega * d_omega_s
            # Derivative of the extraterrestrial radiation with respect to time
            d_Gon_dn = self.Gon(n, deriv=True)
            # Total derivative (scaled)
            return A * (B_n * d_Gon_dn + Gon * d_B_n_d_n)

    def sun_hours(self, n, lat, alt):
        """Compute number of sun hours for day n."""
        sun_hrs = self.sunset_hour_angle(n, lat, alt) * (2 / 15)
        return pd.to_timedelta(sun_hrs, unit='h')

    def solar_altitude(self, n, lat):
        """Compute solar altitude in degrees."""
        phi = self.radiant(lat)
        delta = self.radiant(self.declination(n))
        arg = np.sin(delta) * np.sin(phi) + np.cos(delta) * np.cos(phi)
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
        end_date = start_date + pd.Timedelta(days=1) # Simplified and corrected
        
        # Format pandas freq from R's `1 min` to `1min` string equivalent
        freq = str(by).replace(" ", "")
        day_date_seq = pd.date_range(start=start_date, end=end_date, freq=freq)
    
        # Solar time
        LST = self.solar_time(day_date_seq, lon, lon_st, tz=tz)
        
        output = pd.DataFrame({
            'date': day_date,
            'clocktime': day_date_seq,
            'solartime': LST,
            'lat': lat,
            'lon': lon,
            'omega': self.hour_angle(LST), # Solar angle
            'declination': self.declination(day_date), # Declination in radiant
            'omega_min': -self.sunset_hour_angle(day_date, lat, alt), # Sunrise hour angle in radiant
            'omega_max': self.sunset_hour_angle(day_date, lat, alt), # Sunset hour angle in radiant
            'sun_hours': self.sun_hours(day_date, lat, alt), # Number of sun hours
            'theta': self.incidence_angle(LST, lat, alt, beta, gamma), # Incidence angle
            'E': self.E(day_date) # Time adjustment in seconds
        })
        return output

    # =================================================================
    # Hottel Clear Sky Model
    # =================================================================

    def clearsky(self, cosZ=None, G0=None, alt=None, clime="No Correction"):
        """Hottel clearsky model."""
        valid_climes = ["No Correction", "Summer", "Winter", "Subartic Summer", "Tropical"]
        
        # Handle match.arg logic from R (uses first element if vector/list is passed)
        if isinstance(clime, (list, tuple, np.ndarray, pd.Series)):
            clime = clime[0]
        if clime not in valid_climes:
            raise ValueError(f"clime must be one of {valid_climes}")

        # Altitude must be converted from metre to km
        altitude = 0 if alt is None else alt

        if altitude > 2.5:
            a0_star = 0.6 * (1 - np.exp(-0.214 * (altitude - 1.12)))
        else:
            a0_star = 0.4237 - 0.00821 * (6.0 - altitude)**2

        a1_star = 0.5055 - 0.00595 * (6.5 - altitude)**2
        a2_star = 0.2711 - 0.01858 * (2.5 - altitude)**2
        
        # Correction factors
        factors = {
            "Summer": [0.97, 0.99, 1.02],
            "Winter": [1.03, 1.01, 1.00],
            "Subartic Summer": [0.99, 0.99, 1.01],
            "Tropical": [0.95, 0.98, 1.02],
            "No Correction": [1.0, 1.0, 1.0]
        }
        
        a = np.array([a0_star, a1_star, a2_star]) * np.array(factors[clime])
        
        tau_beam = a[0] + a[1] * np.exp(-a[2] / cosZ)
        tau_diffuse = 0.271 - 0.294 * tau_beam
        
        output = pd.DataFrame({'tau_beam': tau_beam, 'tau_diffuse': tau_diffuse})
        
        output['tau_beam'] = np.where((output['tau_beam'] > 1) | (output['tau_beam'] < 0), 0, output['tau_beam'])
        output['tau_diffuse'] = np.where((output['tau_diffuse'] > 1) | (output['tau_diffuse'] < 0), 0, output['tau_diffuse'])
        
        skymax = G0 * output['tau_beam'] + G0 * output['tau_diffuse']
        return skymax