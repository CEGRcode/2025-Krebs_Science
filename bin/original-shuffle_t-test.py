import pandas as pd
from scipy import stats
import os
import argparse

# Set up argument parser to accept folder path and output file as arguments
parser = argparse.ArgumentParser(description="Process CSV files and perform T-tests")
parser.add_argument('folder_path', type=str, help="Path to the folder containing CSV files")
parser.add_argument('output_file', type=str, help="Output file to save the results")
args = parser.parse_args()

# Initialize a list to store the results
results = []

# Loop through all CSV files in the folder
for filename in os.listdir(args.folder_path):
    if filename.endswith('_score_peak.csv'):
        file_path = os.path.join(args.folder_path, filename)
        
        try:
            # Read the CSV file into a DataFrame
            data = pd.read_csv(file_path)

            # Clean column names (strip whitespace)
            data.columns = data.columns.str.strip()

            # Check for required columns
            if 'enrichment' not in data.columns or 'Nucleosomephase' not in data.columns:
                print(f"Warning: Missing required columns in {filename}. Skipping.")
                continue  # Skip this file

            # Clean the data (ensure 'enrichment' is numeric)
            data['enrichment'] = pd.to_numeric(data['enrichment'], errors='coerce')

            # Drop rows where 'enrichment' is NaN
            data = data.dropna(subset=['enrichment'])

            # Loop through each pair of locations and perform t-test
            locations = data['Nucleosomephase'].unique()

            if len(locations) == 2:  # Ensure there are exactly two groups for comparison
                group1 = data[data['Nucleosomephase'] == locations[0]]['enrichment']
                group2 = data[data['Nucleosomephase'] == locations[1]]['enrichment']

                # Ensure both groups have at least two data points
                if len(group1) < 2 or len(group2) < 2:
                    print(f"Warning: Not enough data in groups for t-test in {filename}. Skipping.")
                    continue

                # Perform the t-test (Welch's t-test for unequal variances)
                t_stat, p_value = stats.ttest_ind(group1, group2, equal_var=False)

                # Format the p-value in scientific notation (e.g., 1.23e-10)
                formatted_p_value = f"{p_value:.3e}"

                # Check for significance and report result
                alpha_05 = 0.05  # significance threshold for *
                alpha_01 = 0.01  # significance threshold for **
                alpha_001 = 0.001  # significance threshold for ***

                if p_value < alpha_001:
                    significance = '***'  # p-value < 0.001
                elif p_value < alpha_01:
                    significance = '**'   # p-value < 0.01
                elif p_value < alpha_05:
                    significance = '*'    # p-value < 0.05
                else:
                    significance = 'ns'   # p-value >= 0.05 (not significant)

                # Append the results (with formatted p-value in scientific notation)
                results.append([filename, locations[0], locations[1], t_stat, formatted_p_value, significance])

        except Exception as e:
            print(f"Error processing file {filename}: {e}")
            continue  # Skip this file if there's an error

# Convert the results list to a DataFrame (without Q and Phase)
results_df = pd.DataFrame(results, columns=['Filename', 'Nucleosomephase 1', 'Nucleosomephase 2', 'T-statistic', 'P-value', 'Significance'])

# Save the results to the specified output file
if not os.path.exists(args.output_file):
    results_df.to_csv(args.output_file, index=False)  # Write header if file does not exist
else:
    results_df.to_csv(args.output_file, index=False, mode='a', header=False)  # Append to file if exists
