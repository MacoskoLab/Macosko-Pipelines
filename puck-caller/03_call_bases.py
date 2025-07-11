"""
This script does the following:
1) Read in master mask and puck base images (harmonized orientation)
2) Call bases
3) Output: spatial barcode whitelist

@author: kimgaram
"""
# read in files that were previously generated
import pickle
with open('~/dir/expanded_labeled_mask.pkl', 'rb') as file:
    expanded_labeled_mask = pickle.load(file)
with open('~/dir/transformed_images.pkl', 'rb') as file:
    transformed_images = pickle.load(file)  
    
%%  
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# make binary mask of segmentations
import cv2
import numpy as np                   #lib for numerial computations
import matplotlib.pyplot as plt      #lib for plotting / displaying images
from skimage import io
from skimage.filters import threshold_otsu, threshold_multiotsu

binary_masks = {}

for base_var, image in transformed_images.items():
    
    # Get the segmentation mask and the original image for the current key
    segmentation_mask = expanded_labeled_mask
    original_image = transformed_images[base_var]
    
    # Apply filter out background, i.e. where signal is 0
    int_filter = np.ma.masked_where(segmentation_mask == 0, original_image[:, :, 0])

    # Store binary mask in 'binary_masks' dictionary
    binary_masks[base_var] = segmentation_mask>0


# Display a test binary 
plt.imshow(binary_masks['base01'][3000:3500, 3000:3500], cmap='gray')
plt.title('Binary Mask Image')
plt.axis('off')
plt.show()


# Save the dictionary to a file
import pickle
with open('~/dir/binary_masks_ALL_transformed.pkl', 'wb') as file:
    pickle.dump(binary_masks, file)

#%%
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Run regionprops / save highest percentile
 
import pickle
import pandas as pd
import numpy as np
from skimage.measure import regionprops_table    
from scipy.stats import percentileofscore

# Functions for base calling and percentile rank
def basecall_percentile_vectorized(intensity_means, flattened_data_sorted, precomputed_percentiles):
    # Use searchsorted to find indices for intensity means in the sorted flattened data
    indices = np.searchsorted(flattened_data_sorted, intensity_means, side='right')
    
    # Clip the indices to ensure they are within the valid range
    indices = np.clip(indices, 0, len(precomputed_percentiles) - 1)
    
    # Retrieve percentiles using clipped indices
    percentiles = precomputed_percentiles[indices]
    return percentiles

# Dictionary to store DataFrames for each base (image)
base_call_dfs = {}

# Iterate through each base image and corresponding binary mask
for base_var, binary_mask in binary_masks.items():
    if base_var in transformed_images:
        # Get the corrected image and segmentation
        corrected_image = transformed_images[base_var]
        segmentation = expanded_labeled_mask
    
    
        # Initialize individual DataFrames for each color channel
        df_Cy5_T, df_TexasRed_G, df_Cy3_C, df_FAM_A = None, None, None, None
    
        # Loop through each channel (Cy5, TexasRed, Cy3, FAM)
        for i, channel_name in enumerate(['Cy5 (T)', 'TexasRed (G)', 'Cy3 (C)', 'FAM (A)']):
            channel = corrected_image[:, :, i]  # Extract the channel
            print(f'{base_var} {channel_name} image read in')
            
            # print(f"{base_var} corrected_image shape: {corrected_image.shape}")
            # print(f"{base_var} segmentation shape: {segmentation.shape}")

            # Use regionprops_table to get centroids and mean intensities for all regions
            props = regionprops_table(segmentation, intensity_image=channel, 
                                      properties=['label', 'intensity_mean', 'centroid'])
            print(f'{base_var} {channel_name} regionprops complete')
            
            # Convert to DataFrame
            df_channel = pd.DataFrame.from_dict(props)
    
        
            # Flatten the channel image to get all pixel intensities within the mask
            flattened_data = channel[binary_mask].flatten()
            flattened_data_sorted = np.sort(flattened_data)
            precomputed_percentiles = np.linspace(0, 100, len(flattened_data_sorted))

            print(f'{base_var} {channel_name} data flattened and percentiles precomputed')
    
            # Apply vectorized base calling
            intensity_means = df_channel['intensity_mean'].values
            df_channel[f'{channel_name}_pct_threshold'] = basecall_percentile_vectorized(
                intensity_means, flattened_data_sorted, precomputed_percentiles
            )
            
            print(f'{base_var} {channel_name} base calling and percentile calculation complete')
            
        
            # Store the DataFrame for the corresponding channel
            if channel_name == 'Cy5 (T)':
                df_Cy5_T = df_channel
            elif channel_name == 'TexasRed (G)':
                df_TexasRed_G = df_channel
            elif channel_name == 'Cy3 (C)':
                df_Cy3_C = df_channel
            elif channel_name == 'FAM (A)':
                df_FAM_A = df_channel
    
        # Set index by 'label' and rename columns appropriately for merging
        df_Cy5_T = df_Cy5_T.set_index('label')
        #print(df_Cy5_T.columns)
        df_Cy5_T.columns = ['intensity_mean_Cy5(T)', 'centroid_x', 'centroid_y', 'pct_Cy5(T)']
    
        df_TexasRed_G = df_TexasRed_G.set_index('label')
        df_TexasRed_G.columns = ['intensity_mean_TexasRed(G)', 'centroid_x_TexasRed(G)', 'centroid_y_TexasRed(G)', 'pct_TexasRed(G)']
    
        df_Cy3_C = df_Cy3_C.set_index('label')
        df_Cy3_C.columns = ['intensity_mean_Cy3(C)', 'centroid_x_Cy3(C)', 'centroid_y_Cy3(C)', 'pct_Cy3(C)']
    
        df_FAM_A = df_FAM_A.set_index('label')
        df_FAM_A.columns = ['intensity_mean_FAM(A)', 'centroid_x_FAM(A)', 'centroid_y_FAM(A)', 'pct_FAM(A)']
    
        # Merge all the DataFrames (for each color channel)
        merged = df_Cy5_T.merge(df_FAM_A, right_index=True, left_index=True)
        merged = merged.merge(df_TexasRed_G, right_index=True, left_index=True)
        merged = merged.merge(df_Cy3_C, right_index=True, left_index=True)
    
        # Calculate the sum of bases called for each segmentation
        merged['sum'] = (merged['pct_Cy5(T)'] > 0) + (merged['pct_TexasRed(G)'] > 0) + (merged['pct_Cy3(C)'] > 0) + (merged['pct_FAM(A)'] > 0)
    
        # Calculate base calling percentages based on percentile scores
        merged['pct_call'] = merged[['pct_Cy5(T)', 'pct_TexasRed(G)', 'pct_Cy3(C)', 'pct_FAM(A)']].idxmax(axis=1)
    
        # Drop unnecessary 'centroid...' columns except for centroid_x and centroid_y
        centroid_columns_to_keep = ['centroid_x', 'centroid_y']
        merged = merged.loc[:, ~merged.columns.str.startswith('centroid_') | merged.columns.isin(centroid_columns_to_keep)]
        
        # Calculate largest and second-largest percentage values
        merged['largest_pct_value'] = merged[['pct_Cy5(T)', 'pct_TexasRed(G)', 'pct_Cy3(C)', 'pct_FAM(A)']].apply(lambda row: row.nlargest(1).iloc[0], axis=1)
        merged['second_pct_value'] = merged[['pct_Cy5(T)', 'pct_TexasRed(G)', 'pct_Cy3(C)', 'pct_FAM(A)']].apply(lambda row: row.nlargest(2).iloc[1], axis=1)
        
        # Calculate the difference between the largest and second-largest percentages
        merged['pct_diff_top2'] = merged['largest_pct_value'] - merged['second_pct_value']
        
        # Create a histogram for the difference in percentages
        plt.hist(merged['pct_diff_top2'], bins=1000, range=(0, 100), color='gray')
        plt.xlabel('Difference between largest and second-largest pct')
        plt.ylabel('Frequency')
        plt.title(f'{base_var}: Histogram of Percentage Differences between Top 2 Values')
        plt.show()
    
        # Store the final merged DataFrame for this image in the dictionary
        base_call_dfs[base_var] = merged
    
    else:
        print(f"Warning: {base_var} not found in corrected_base_images or segmentations")



# Save dictionary to a file
import pickle
with open('~/dir/base_call_dfs.pkl', 'wb') as file:
    pickle.dump(base_call_dfs, file)


# to load dictionary from file:
import pickle
with open('~/dir/base_call_dfs', 'rb') as file:
    base_call_dfs = pickle.load(file)
    

#%%
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# Call an actual base for each centroid

#Assign as follows: 
    #if pct_call = 'pct_Cy3(C)', then base_call = 'C'
    #if pct_call = 'pct_Cy5(T)', then base_call = 'T'
    #if pct_call = 'pct_TexasRed(G) ', then base_call = 'G'
    #if pct_call = 'pct_FAM(A)', then base_call = 'A'
    #if pct_diff_top2 < 15, then base_call = 'N'
    
base_call_dfs_FINAL_ALL = base_call_dfs
# Iterate over each DataFrame in the dictionary
for base_var, df in base_call_dfs_FINAL_ALL.items():
    # Apply conditions to assign base_call values
    df['base_call'] = df['pct_call'].apply(
        lambda x: 'C' if x == 'pct_Cy3(C)' else
                  'T' if x == 'pct_Cy5(T)' else
                  'G' if x == 'pct_TexasRed(G)' else
                  'A' if x == 'pct_FAM(A)' else 'N'
    )
    
    # Apply additional condition for pct_diff_top2 #CAN CHANGE THIS VALUE AS DESIRED
    df.loc[df['pct_diff_top2'] < 1, 'base_call'] = 'N'

    # Update the dictionary with the modified DataFrame
    base_call_dfs_FINAL_ALL[base_var] = df


# Look at base distribution at each position
# output TABLES and a BARPLOT
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# List of consistent base call order
base_order = ['A', 'C', 'G', 'T', 'N']

# Colors for each base
base_colors = {
    'A': '#FF9999',  # Light red
    'C': '#99CCFF',  # Light blue
    'G': '#99FF99',  # Light green
    'T': '#CC99FF',  # Light purple
    'N': '#FFCC99'   # Light orange
}

# Dictionary to store the summary tables for each base_var
summary_tables = {}

# List to store the counts and base_vars for plotting
plot_data = {base: [] for base in base_order}  # Store base percentages across positions

# Iterate over each DataFrame in the dictionary
for base_var, df in base_call_dfs_FINAL_ALL.items():
    # Count occurrences of each base_call and reindex to ensure consistent order
    base_call_counts = df['base_call'].value_counts().reindex(base_order, fill_value=0)
    
    # Calculate percentage
    base_call_percentage = (base_call_counts / len(df)) * 100
    
    # Create a summary DataFrame
    summary_df = pd.DataFrame({
        'base_call': base_call_counts.index,
        'count': base_call_counts.values,
        'percentage': base_call_percentage.values
    })
    
    # Store the summary table in the dictionary
    summary_tables[base_var] = summary_df
    
    # Append data for plotting, but group by base, not position
    for base in base_order:
        plot_data[base].append(summary_df.loc[summary_df['base_call'] == base, 'percentage'].values[0])

# Convert plot_data into a format suitable for grouped bar plot
num_positions = len(base_call_dfs_FINAL_ALL)
labels = list(base_call_dfs_FINAL_ALL.keys())  # The label locations for the positions
x = np.arange(num_positions)  # The label locations for the positions
width = 0.15  # Width of the bars

# Create a grouped bar chart
fig, ax = plt.subplots(figsize=(16, 8))  # Increased figure size for clarity

# Plot each base's counts across positions with a consistent offset for grouped bars
for i, base in enumerate(base_order):
    ax.bar(x + i * width, plot_data[base], width, color=base_colors[base], label=base)

# Annotate each bar with its percentage value
for i, base in enumerate(base_order):
    for j, percent in enumerate(plot_data[base]):
        if percent > 0:  # Only annotate bars with a percentage over 5 to reduce clutter
            ax.text(x[j] + i * width, percent + 0.5, f'{percent:.1f}%', ha='center', va='bottom')

# Add labels, title, and legend
ax.set_xlabel('Position')
ax.set_ylabel('Percentage')
ax.set_title('Base Percentage Distribution at Each Position')
ax.set_xticks(x + width * (len(base_order) - 1) / 2)  # Adjust xticks to align
ax.set_xticklabels(labels)
ax.legend(title='Base', bbox_to_anchor=(1.05, 1), loc='upper left')

# Rotate x-tick labels for clarity
plt.xticks(rotation=0)


# Show plot
plt.tight_layout()
plt.show()


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

#%%
def create_final_output(base_calls_with_coords):
    merged_df = None

    # Iterate through each base DataFrame and align them on 'centroid_x' and 'centroid_y'
    for base_var, df in base_calls_with_coords.items():
        # Rename the 'base_call' column to include the base name for distinction
        df_renamed = df.rename(columns={'base_call': f'base_call_{base_var}'})
        
        if merged_df is None:
            # For the first base, initialize the merged_df with its data
            merged_df = df_renamed[['centroid_x', 'centroid_y', f'base_call_{base_var}']]
        else:
            # Merge subsequent DataFrames on 'centroid_x' and 'centroid_y'
            merged_df = pd.merge(merged_df, df_renamed[['centroid_x', 'centroid_y', f'base_call_{base_var}']],
                                 on=['centroid_x', 'centroid_y'], how='outer')

    return merged_df


# run the function
final_output = create_final_output(base_call_dfs_FINAL_ALL)

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#%%
# FINALLY, create and save the sb_whitelist
def generate_sb_whitelist(final_output, total_bases=14):
    # Step 1: Initialize the list to store spatial barcodes
    sb_whitelist_data = []

    # Step 2: Iterate through each row of final_output
    for _, row in final_output.iterrows():
        # Initialize an array for the 14-base string, defaulting to 'N'
        sb_bases = ['N'] * total_bases

        # Iterate through the base columns to construct the spatial barcode
        for i in range(1, total_bases + 1):  # base01 to base14
            base_column = f'base_call_base{i:02d}'
            if base_column in row:
                sb_bases[i - 1] = row[base_column]  # Place the base call at the correct position

        # Concatenate the 14 bases into the spatial barcode (sb)
        sb = ''.join(sb_bases)

        # Append the data to the list
        sb_whitelist_data.append({
            'sb': sb,
            'x': row['centroid_x'],
            'y': row['centroid_y']
        })

    # Convert the list into a DataFrame
    sb_whitelist_df = pd.DataFrame(sb_whitelist_data)

    return sb_whitelist_df


# Assuming 'final_output' is already created
sb_whitelist = generate_sb_whitelist(final_output)
sb_whitelist




# Save to CSV if needed
# sb_whitelist.to_csv('~/dir/Puck_name.csv', index=False)


""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

# how many duplicated 'sb' values are there???
def find_duplicate_sb_counts(sb_whitelist):
    # Find duplicated 'sb' values
    duplicated_sb = sb_whitelist[sb_whitelist.duplicated('sb', keep=False)]

    # Count the number of occurrences for each duplicated 'sb'
    sb_counts = duplicated_sb['sb'].value_counts()

    # Filter out 'sb' values that appear only once (unique values)
    duplicated_sb_counts = sb_counts[sb_counts > 1]

    return duplicated_sb_counts

# Example usage:
duplicate_counts = find_duplicate_sb_counts(sb_whitelist)
print(f"Number of duplicated 'sb' instances: {len(duplicate_counts)}")


#Investigate: distance between duplicates (i.e. is it a bead that got counted twice?)
def compute_duplicate_distances(sb_whitelist):
    """
    For each duplicated 'sb', compute the Euclidean distance between consecutive duplicate rows.
    Returns a DataFrame with each pair and its distance.
    """
    # Filter rows with duplicated 'sb' (keeping all occurrences)
    duplicated_df = sb_whitelist[sb_whitelist.duplicated('sb', keep=False)]
    
    # Group by 'sb' so we process each set of duplicates together
    duplicate_groups = duplicated_df.groupby('sb')
    
    # List to store computed distances
    rows = []
    for sb, group in duplicate_groups:
        # Sort by index (or any other criterion relevant to your data)
        group = group.sort_index()
        # Loop over consecutive rows within the group
        for i in range(len(group) - 1):
            # Get coordinates for the pair of rows
            x1, y1 = group.iloc[i]['x'], group.iloc[i]['y']
            x2, y2 = group.iloc[i+1]['x'], group.iloc[i+1]['y']
            # Compute Euclidean distance
            distance = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
            rows.append({
                'sb': sb,
                'index1': group.index[i],
                'index2': group.index[i+1],
                'x1': x1,
                'y1': y1,
                'x2': x2,
                'y2': y2,
                'distance': distance
            })
    # Create a DataFrame from the list
    duplicates_table = pd.DataFrame(rows)
    # Sort by distance (from smallest to largest)
    duplicates_table_sorted = duplicates_table.sort_values(by='distance').reset_index(drop=True)
    return duplicates_table_sorted

# Compute and preview the duplicates table with distances
duplicates_table = compute_duplicate_distances(sb_whitelist)
print("Table of duplicates with distances (sorted by distance):")
print(duplicates_table.head(250)) #if within distance 17, merge

###############################################################################
# Optional: Merge duplicates if they are immediately next to each other
###############################################################################
def merge_adjacent_duplicates(sb_whitelist, distance_threshold):
    """
    Merge duplicate rows for each 'sb' if the distance between consecutive rows 
    is less than or equal to the specified threshold.
    For merged groups, the x and y values are replaced by their mean.
    """
    merged_rows = []
    # Process each barcode group
    for sb, group in sb_whitelist.groupby('sb'):
        # Sort the group (using index; change if another order is more meaningful)
        group = group.sort_index().reset_index(drop=True)
        if len(group) > 1:
            # Compute distances between consecutive rows
            distances = np.sqrt(np.diff(group['x'])**2 + np.diff(group['y'])**2)
            # If all consecutive distances are within the threshold, merge them into one entry
            if (distances <= distance_threshold).all():
                merged_x = group['x'].mean()
                merged_y = group['y'].mean()
                merged_rows.append({'sb': sb, 'x': merged_x, 'y': merged_y})
            else:
                # Otherwise, keep rows individually (or you could refine this logic to merge only close pairs)
                for _, row in group.iterrows():
                    merged_rows.append({'sb': row['sb'], 'x': row['x'], 'y': row['y']})
        else:
            # For unique barcodes, just keep the row as is
            row = group.iloc[0]
            merged_rows.append({'sb': row['sb'], 'x': row['x'], 'y': row['y']})
    
    return pd.DataFrame(merged_rows)

# Example: merge duplicates if their consecutive distance is <= a chosen threshold (e.g., 5)
distance_threshold = 17  # You can adjust this threshold based on your criteria
merged_sb_whitelist = merge_adjacent_duplicates(sb_whitelist, distance_threshold)
print("\nMerged whitelist (if duplicates are within the distance threshold):")
print(merged_sb_whitelist.head())

#
def count_sbs_with_more_than_two_Ns(sb_whitelist):
    # Create a boolean mask where the number of 'N's in 'sb' is greater than 2
    more_than_two_Ns = sb_whitelist['sb'].apply(lambda x: x.count('N') > 2)
    
    # Count the number of rows that match the condition
    num_rows_with_more_than_two_Ns = more_than_two_Ns.sum()
    
    return num_rows_with_more_than_two_Ns

# Example usage:
num_rows_with_more_than_two_Ns = count_sbs_with_more_than_two_Ns(sb_whitelist)
print(f"Number of rows with more than 2 'N's in 'sb': {num_rows_with_more_than_two_Ns}")

## Drop rows with duplicate 'sb' values
def filter_unique_sb_whitelist(sb_whitelist):
    # Drop duplicates based on the 'sb' column, keeping only the first occurrence
    sb_whitelist_unique = sb_whitelist.drop_duplicates(subset='sb', keep='first').reset_index(drop=True)
    
    return sb_whitelist_unique

# Apply:
sb_whitelist_unique = filter_unique_sb_whitelist(sb_whitelist)
print(sb_whitelist_unique) 

# Optionally, save result back to a CSV
sb_whitelist_unique.to_csv('~/dir/Puck_name.csv', header=False, index=False, sep=',')


# To import back in:
import pandas as pd
file_path = '~/dir/Puck_name.csv'
sb_whitelist_unique = pd.read_csv(file_path, header=None)
sb_whitelist_unique.columns = ['sb', 'x', 'y']
sb_whitelist_unique