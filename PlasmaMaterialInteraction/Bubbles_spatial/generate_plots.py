import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Define timestamps and distances
timestamps = [2e-6, 4e-6]
distances = np.linspace(0, 4e-7, num=20) # Equidistant distance array

# Load the Excel file
file_path = "./results_spatial.xlsx"  # Update this with your file path
output_dir = "distance_plots"
os.makedirs(output_dir, exist_ok=True)  # Create output directory if it doesn't exist

# Read available sheet names
xls = pd.ExcelFile(file_path)
sheet_names = sorted([name for name in xls.sheet_names])  # Convert to float and sort

# Function to find the closest sheet name
def find_closest_sheet(timestamp, sheet_names):
    return min(sheet_names, key=lambda x: abs(float(x) - timestamp))

# Process each timestamp
for timestamp in timestamps:
    closest_sheet = find_closest_sheet(timestamp, sheet_names)
    sheet_name = str(closest_sheet)  # Convert back to string for lookup

    try:
        df = pd.read_excel(file_path, sheet_name=sheet_name, header=None)  # Read closest sheet
        first_row = df.iloc[0, :len(distances)]  # Get first row data

        # Plot data
        plt.figure(figsize=(8, 5))
        plt.plot(distances, first_row, marker='o', linestyle='-')
        plt.xlabel("Distance")
        plt.ylabel("Value")
        plt.title(f"t = {closest_sheet}")
        plt.grid(True)

        # Save plot using the exact sheet name
        plot_filename = os.path.join(output_dir, f"{sheet_name}.png")
        plt.savefig(plot_filename)
        plt.close()
        print(f"Saved plot: {plot_filename}")

    except Exception as e:
        print(f"Error processing sheet '{sheet_name}': {e}")

print("All plots generated successfully.")
