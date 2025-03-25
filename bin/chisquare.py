import pandas as pd
from scipy import stats
import os
import argparse

# Set up argument parser to accept folder path and output file as arguments
parser = argparse.ArgumentParser(description="Process CSV files and perform Chi-square tests")
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

            # Check for required columns
            if 'enrichment' not in data.columns or 'Nucleosomephase' not in data.columns:
                print(f"Warning: Missing required columns in {filename}. Skipping.")
                continue  # Skip this file

            # Clean the data (ensure 'enrichment' is numeric)
            data['enrichment'] = pd.to_numeric(data['enrichment'], errors='coerce')

            # Drop rows where 'enrichment' is NaN
            data = data.dropna(subset=['enrichment'])

            # Bin the enrichment data into categories (low, medium, high)
            bins = [0, 0.1, 0.3, 1]  # Example: low (0-0.1), medium (0.1-0.3), high (0.3-1)
            labels = ['Low', 'Medium', 'High']
            data['enrichment_binned'] = pd.cut(data['enrichment'], bins=bins, labels=labels, right=False)

            # Create a contingency table for the Chi-square test
            contingency_table = pd.crosstab(data['Nucleosomephase'], data['enrichment_binned'])

            # Perform the Chi-square test
            chi2_stat, p_value, dof, expected = stats.chi2_contingency(contingency_table)

            # Format the p-value in scientific notation (e.g., 1.23e-10)
            formatted_p_value = f"{p_value:.3e}"  # Format the p-value to scientific notation

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
            results.append([filename, chi2_stat, formatted_p_value, dof, significance])

        except Exception as e:
            print(f"Error processing file {filename}: {e}")
            continue  # Skip this file if there's an error

# Convert the results list to a DataFrame (without Q and Phase)
results_df = pd.DataFrame(results, columns=['Filename', 'Chi-squared Statistic', 'P-value', 'Degrees of Freedom', 'Significance'])

# Save the results to the specified output file
with open(args.output_file, 'a') as f:  # 'a' to append to the file
    results_df.to_csv(f, index=False, header=f.tell() == 0)  # Write header only if file is empty
