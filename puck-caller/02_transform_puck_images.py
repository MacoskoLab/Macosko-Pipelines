"""
This script does the following:
1) Read in ICP transformation matrices
2) Transform all original puck base images to be in same orientation/coord system
    
@author: kimgaram
"""

#%%
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#read in corrected_base_images
import pickle
with open('~/dir/corrected_base_images.pkl', 'rb') as file:
    corrected_base_images = pickle.load(file)    


#read in ICP_transformation_matrices
import pickle
with open('~/dir/ICP_transformation_matrices.pkl', 'rb') as file:
    ICP_transformation_matrices = pickle.load(file)
    
    
#look at current ICP transformation matrix structure
ICP_transformation_matrices['base02'].shape #(4,4)


import numpy as np 
def condense_matrix_4x4_to_3x3(matrix_4x4):
    # Extract the top-left 2x2 block (rotation + scaling)
    rotation_scaling = matrix_4x4[:2, :2]
    
    # Extract the first two elements from the last column (translation)
    translation = matrix_4x4[:2, 3]
    
    # Create the condensed 3x3 matrix
    condensed_matrix = np.eye(3)  # Start with an identity 3x3 matrix
    condensed_matrix[:2, :2] = rotation_scaling
    condensed_matrix[:2, 2] = translation
    
    return condensed_matrix

# Apply the function to condense all transformation matrices
condensed_matrices = {base: condense_matrix_4x4_to_3x3(matrix) for base, matrix in ICP_transformation_matrices.items()}

# Example of one condensed matrix
print(condensed_matrices['base02'])

condensed_matrices['base02'].shape #(3,3)

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#%%
"""Apply ICP_transformation_matrices to each image and channel"""
# corrected_base_images # original images, 4 channels each
# ICP_transformation_matrices # 4x4 transformation matrix
# condensed_matrices #3x3 transformation matrix
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
import numpy as np
import matplotlib.pyplot as plt

transformed_images = {}

for base_var, corrected_image in corrected_base_images.items():
    if base_var == 'base02':
        transformed_images[base_var] = corrected_image
    else:
        img = corrected_image
        transformation_matrix = condensed_matrices[base_var]  # 3x3
        
        # Use only the first 2 rows for affine transformation
        affine_matrix = transformation_matrix[:2, :]
        
        # Initialize a placeholder for the transformed image
        transformed_img = np.zeros_like(img)  # This will be a 3D array (height, width, channels)
        
        # Generate a grid of x and y coordinates
        x, y = np.mgrid[0:img.shape[0], 0:img.shape[1]]
        coords = np.stack([x.ravel(), y.ravel(), np.ones_like(x.ravel())], axis=0)  # Shape: (3, num_pixels)
        
        # Apply the affine transformation to all coordinates at once
        transformed_coords = np.dot(affine_matrix, coords).astype(int)  # Shape: (2, num_pixels)
        
        # Keep only valid coordinates (inside image boundaries)
        valid_mask = (
            (transformed_coords[0] >= 0) & (transformed_coords[0] < img.shape[0]) &
            (transformed_coords[1] >= 0) & (transformed_coords[1] < img.shape[1])
        )
        
        # Reshape to 2D
        transformed_coords = transformed_coords[:, valid_mask]
        orig_coords = coords[:2, valid_mask]  # Corresponding original coordinates
    
        # Apply the transformation to each channel
        for i, channel_name in enumerate(['Cy5 (T)', 'TexasRed (G)', 'Cy3 (C)', 'FAM (A)']):
            print(f"Applying transformation to {channel_name} channel of {base_var}")
            transformed_img[transformed_coords[0], transformed_coords[1], i] = img[orig_coords[0], orig_coords[1], i]
        
        # Store the transformed image
        transformed_images[base_var] = transformed_img
    
        # Display the transformed image using matplotlib
        plt.imshow(transformed_img[:, :, 0], cmap='gray')  # Display the first channel (Cy5 (T))
        plt.title(f'Transformed Image - {base_var} (Cy5 (T))')
        plt.show()
        
        # Visualize original vs transformed images, zooming into a 500x500 region for the first channel (Cy5)
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Original (Cy5 channel)
        axes[0].imshow(img[3000:3500, 3000:3500, 0], cmap='gray')  # Display the first channel (Cy5 (T))
        axes[0].set_title(f'{base_var} - Original Cy5 (500x500)')
        
        # Transformed (Cy5 channel)
        axes[1].imshow(transformed_img[3000:3500, 3000:3500, 0], cmap='gray')  # Display the first channel (Cy5 (T))
        axes[1].set_title(f'{base_var} - Transformed Cy5 (500x500)')
        
        plt.show()


# to save dictionary to a file:
import pickle
with open('~/dir/transformed_images.pkl', 'wb') as file:
    pickle.dump(transformed_images, file)


# to load dictionary from file:
import pickle
with open('~/dir/transformed_images.pkl', 'rb') as file:
    transformed_images = pickle.load(file)
   