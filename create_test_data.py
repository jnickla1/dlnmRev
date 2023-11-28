import pandas as pd
import numpy as np


# Assuming your data is loaded from a CSV file named 'your_data.csv'
data = pd.read_csv('~/Documents/dlnmRevData/Mar8_envs_tem.csv')

# Remove rows with any NaN values
data.dropna(inplace=True)

# Create a new column 'Outcome' with all 0s
data['Outcome'] = 0



def sample_flat_data(data):
    sampled_data = data.copy()  # Create a copy of the original data

    for index, row in sampled_data.iterrows():
        for col in np.arange(4,30):  # Assuming columns are named t0, t1, ... t29
            if np.random.rand() < 0.005:  # 0.5% chance of sampling
                sampled_data.at[index, 'Outcome'] = 1  # Put 1 in the 'Outcome' column

                # Replace cells after the sampled time with NaN
                for later_col in range(col + 1, 30):
                    sampled_data.at[index, f't{later_col}'] = np.nan
                break  # Move to the next row after sampling one row

    return sampled_data

# Assuming 'modified_data.csv' contains the modified data from the previous step
sampled_data = sample_flat_data(data)


# Save the sampled data to a new CSV file
sampled_data.to_csv('~/Documents/dlnmRevData/sampled_flat_data.csv', index=False)
