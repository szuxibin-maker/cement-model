import time
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

def measure_and_plot_time(model, data_dict, model_name):
    sizes = []
    times = []

    for label, data in data_dict.items():
        # Check if data is a DataFrame, a dictionary, or a list of dictionaries
        if isinstance(data, pd.DataFrame):
            data_size = len(data)
        elif isinstance(data, dict):
            data_size = sum(len(df) if isinstance(df, pd.DataFrame) else 1 for df in data.values())
        elif isinstance(data, list) and all(isinstance(df, dict) for df in data):
            data_size = sum(len(df) if isinstance(df, pd.DataFrame) else 1 for df in data)
        else:
            raise ValueError("Invalid data format. Supported formats: DataFrame, dictionary, or list of dictionaries.")

        sizes.append(data_size)

        # Measure the time taken by the model
        start_time = time.time()
        model(data)
        elapsed_time = time.time() - start_time
        times.append(elapsed_time)

    # Plotting
    plt.plot(sizes, times, marker='o')
    plt.title("Model Time Measurement")
    plt.xlabel("Dataset Size")
    plt.ylabel("Time (seconds)")

    # Create the time_plot folder if it doesn't exist
    if not os.path.exists("time_plot"):
        os.makedirs("time_plot")

    # Save the plot with the model name
    model_name = model.__name__
    plot_filename = f"time_plot/{model_name}_time_plot.png"
    plt.savefig(plot_filename)
    plt.show()


    return times
