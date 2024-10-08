import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Read the CSV file into a DataFrame
data = pd.read_csv('Random_NFIA_Nuc_Engagement.csv')

# Filter out transcripts outside the desired range
filtered_data = data[(data['NFIA-Nuc_engagement'] >= 0) & (data['NFIA-Nuc_engagement'] <= 3)]

# Define a custom color palette with alternating colors for categories
custom_palette = ["yellow", "pink"]

# Create the violin plot with the custom color palette
sns.set(style="whitegrid")
plt.figure(figsize=(12, 8))
ax = sns.violinplot(x='Location', y='NFIA-Nuc_engagement', data=filtered_data, palette=custom_palette)

# Plot the median points
#median_data = filtered_data.groupby('Location')['NFIA-Nuc_engagement'].median().reset_index()
#sns.scatterplot(x='Location', y='NFIA-Nuc_engagement', data=median_data, color='black', s=100, label='Median', edgecolor='w')

# Add median annotations
#for index, row in median_data.iterrows():
    #ax.text(row['Location'], row['NFIA-Nuc_engagement'] + 2, f"{row['NFIA-Nuc_engagement']:.1f}", ha='center', va='bottom', color='black')

# Set y-axis limits to 0-100
y_min, y_max = ax.get_ylim()
ax.set_ylim(0, 3)

# Perform pairwise t-tests for 'NFIA-Nuc_engagement' between 'Location' levels
classes = filtered_data['Location'].unique()
x_labels = [label for label in range(len(classes))]  # Get x-axis positions of categories
num_comparisons = len(classes) * (len(classes) - 1) / 2  # Number of pairwise comparisons
alpha = 0.05
adjusted_alpha = alpha / num_comparisons  # Bonferroni correction

# Define a function to calculate a different y-position for each annotation
def get_y_position(index, start_y, spacing=10):
    """Calculate different y-positions for significance annotations starting from a specific y-position."""
    return start_y + index * spacing

# Perform t-tests for all pairs of classes and display significant indicators
annotation_index = 0
start_y = 1  # Starting y-position for annotations
spacing = 1   # Vertical spacing between annotations

for i, class1 in enumerate(classes):
    for j, class2 in enumerate(classes):
        if i < j:
            subset1 = filtered_data[filtered_data['Location'] == class1]['NFIA-Nuc_engagement']
            subset2 = filtered_data[filtered_data['Location'] == class2]['NFIA-Nuc_engagement']
            t_stat, p_value = stats.ttest_ind(subset1, subset2)
            
            # Add significance indicators based on adjusted p-value
            if p_value < adjusted_alpha:
                significance = '***'
            elif p_value < alpha:
                significance = '**'
            else:
                significance = ''
            
            # Add significance indicators to the plot
            x1, x2 = x_labels[i], x_labels[j]
            y_position = get_y_position(annotation_index, start_y, spacing)  # Compute y position
            
            # Ensure the y-position is within the y-axis limits
            if y_position > y_max:
                break  # Exit if the y-position exceeds y-axis limit
            
            y_height = 1.5  # Height of the significance line
            
            # Plot significance line
            plt.plot([x1, x1, x2, x2], [y_position, y_position + y_height, y_position + y_height, y_position], lw=1, c='k')
            plt.text((x1 + x2) * 0.5, y_position + y_height, significance, ha='center', va='bottom', color='k')
            
            annotation_index += 1  # Increment index for each annotation

plt.title('Nucleosome Engagement by Location')
plt.xlabel('Location')
plt.ylabel('NFIA-Nuc Engagement')
plt.legend()  # Show legend for median points
plt.show()
