import pandas as pd
from scipy import stats
import os
import argparse

# Set up argument parser to accept the folder path and output file name as arguments
parser = argparse.ArgumentParser(description="Process CSV files and perform Chi-square tests")
parser.add_argument('folder_path', type=str, help="Path to the folder containing CSV files")
parser.add_argument('output_file', type=str, help="Output file to save the results")
args = parser.parse_args()

# Initialize an empty list to store the results
results = []

# Function to extract Q and phase from the filename
def extract_q_phase(filename):
    # Extract the Q and phase (sense or anti) from the filename
    parts = filename.split('_')
    q = parts[1]  # Extract Q1, Q2, Q3, or Q4
    phase = parts[2]  # sense or anti
    return q, phase

# Loop through all CSV files in the folder
for filename in os.listdir(args.folder_path):
    if filename.endswith('.csv'):
        file_path = os.path.join(args.folder_path, filename)
        
        # Read the CSV file into a DataFrame
        data = pd.read_csv(file_path)
        
        # Check if 'enrichment' column exists
        if 'enrichment' not in data.columns:
            print(f"Warning: 'enrichment' column not found in {filename}. Skipping this file.")
            continue  # Skip this file and move to the next one

        # Clean the data (ensure 'enrichment' is numeric)
        data['enrichment'] = pd.to_numeric(data['enrichment'], errors='coerce')

        # Bin the enrichment data into categories (low, medium, high)
        bins = [0, 0.1, 0.3, 1]  # Example: low (0-0.1), medium (0.1-0.3), high (0.3-1)
        labels = ['Low', 'Medium', 'High']
        data['enrichment_binned'] = pd.cut(data['enrichment'], bins=bins, labels=labels, right=False)

        # Create a contingency table for the Chi-square test
        contingency_table = pd.crosstab(data['Nucleosomephase'], data['enrichment_binned'])

        # Perform the Chi-square test
        chi2_stat, p_value, dof, expected = stats.chi2_contingency(contingency_table)

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

        # Append the results to the list
        results.append([filename, chi2_stat, p_value, dof, significance])

# Convert the results list to a DataFrame for easy viewing
results_df = pd.DataFrame(results, columns=['Filename', 'Chi-squared Statistic', 'P-value', 'Degrees of Freedom', 'Significance'])

# Sort the results by the custom order (Q1, Q2, Q3, Q4 for sense/anti)
results_df['Q'], results_df['Phase'] = zip(*results_df['Filename'].apply(extract_q_phase))

# Define a custom sorting order
order = ['Q1', 'Q2', 'Q3', 'Q4']
results_df['Q_order'] = results_df['Q'].apply(lambda x: order.index(x))

# Sort by Q_order, Phase, and Filename
sorted_results_df = results_df.sort_values(by=['Q_order', 'Phase', 'Filename'])

# Print the final sorted report
print(sorted_results_df)

# Save the results to the specified output file
output_file = args.output_file
sorted_results_df.to_csv(output_file, index=False)
