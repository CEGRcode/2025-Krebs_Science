import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Read the CSV file into a DataFrame
data = pd.read_csv('NFIA_Nuc_Engagement.csv')

# Filter out transcripts outside the desired range
filtered_data = data[(data['NFIA-Nuc_engagement'] >= 0) & (data['NFIA-Nuc_engagement'] <= 3)]

# Bin the 'NFIA-Nuc_engagement' data into categories (low, medium, high)
bins = [0, 1, 2, 3]  # You can change the bin ranges based on your needs
labels = ['Low', 'Medium', 'High']
filtered_data['engagement_binned'] = pd.cut(filtered_data['NFIA-Nuc_engagement'], bins=bins, labels=labels, right=False)

# Define a custom color palette with alternating colors for categories
custom_palette = ["yellow", "pink", "lightblue"]

# Create the violin plot with the custom color palette
sns.set(style="whitegrid")
plt.figure(figsize=(12, 8))
ax = sns.violinplot(x='Location', y='NFIA-Nuc_engagement', data=filtered_data, palette=custom_palette)

# Set y-axis limits to 0-3 (based on the engagement bins)
ax.set_ylim(0, 3)

# Perform Chi-square test of independence for 'Location' vs 'engagement_binned'
# Create a contingency table for the Chi-square test
contingency_table = pd.crosstab(filtered_data['Location'], filtered_data['engagement_binned'])

# Perform the Chi-square test
chi2_stat, p_value, dof, expected = stats.chi2_contingency(contingency_table)

# Format the p-value in scientific notation
formatted_p_value = f"{p_value:.3e}"

# Add significance indicator based on the p-value
alpha = 0.05
if p_value < 0.001:
    significance = '***'
elif p_value < 0.01:
    significance = '**'
elif p_value < 0.05:
    significance = '*'
else:
    significance = 'ns'  # Not significant

# Plot significance result on the plot
plt.title(f'Nucleosome Engagement by Location (Chi-squared p-value: {formatted_p_value})')
plt.xlabel('Location')
plt.ylabel('NFIA-Nuc Engagement')
plt.text(0.5, 2.5, f'Significance: {significance}', ha='center', va='top', color='black', fontsize=14)

# Show the plot
plt.show()

# Print the results of the Chi-square test
print(f"Chi-squared Statistic: {chi2_stat:.3f}")
print(f"P-value: {formatted_p_value}")
print(f"Degrees of Freedom: {dof}")
print(f"Expected Frequencies: \n{expected}")
