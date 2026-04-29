import pandas as pd
import numpy as np

def number_of_day(x):
    """
    Compute the number of the day of the year given a vector of dates or numbers.
    """
    # 1. Check if the input is already numeric
    # Handles single numbers
    if isinstance(x, (int, float, np.number)):
        return x
    
    # Handles Pandas Series or NumPy arrays of numbers
    if isinstance(x, (pd.Series, np.ndarray)) and pd.api.types.is_numeric_dtype(x):
        return x
        
    # 2. Convert to datetime and extract the day of the year
    # pd.to_datetime safely handles single strings, lists of strings, and arrays
    dt_x = pd.to_datetime(x)
    
    # If it's a collection (like a Series or DatetimeIndex), use the .dt accessor
    if hasattr(dt_x, 'dt'):
        return dt_x.dt.dayofyear
    else:
        # If it's a single Timestamp object
        return dt_x.dayofyear