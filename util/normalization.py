import pandas as pd

def minmax_normalize(dataframe):
    """
    Min-Max Normalize a Pandas DataFrame.

    Parameters:
    - dataframe: Pandas DataFrame to be normalized.

    Returns:
    - normalized_dataframe: Min-Max normalized DataFrame.
    - min_values: Dictionary containing minimum values for each column.
    - max_values: Dictionary containing maximum values for each column.
    """
    min_values = dataframe.min()
    max_values = dataframe.max()
    normalized_dataframe = (dataframe - min_values) / (max_values - min_values)
    return normalized_dataframe

def minmax_denormalize(normalized_dataframe, min_values, max_values):
    """
    Convert a Min-Max normalized Pandas DataFrame back to the original scale.

    Parameters:
    - normalized_dataframe: Min-Max normalized DataFrame.
    - min_values: Dictionary containing minimum values for each column.
    - max_values: Dictionary containing maximum values for each column.

    Returns:
    - denormalized_dataframe: DataFrame in the original scale.
    """
    denormalized_dataframe = normalized_dataframe * (max_values - min_values) + min_values
    return denormalized_dataframe

