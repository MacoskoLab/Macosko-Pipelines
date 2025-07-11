"""
This script does the following:
1) Read in all base images
2) Create segmentations from all base images
3) Create bleedthrough-corrected images called `corrected_base_images`
        - this requires creating a `binary mask` from `each segmentation`, eroding the segmentations, then looking at color  
4) Call centroids for each base image
5) [To use in next step]: Generate [ICP] transformation matrices to get all images into same coord system
6) [To use in next step]: Create 'master' mask to use for base calling
    
@author: kimgaram
"""

#order of spatial barcode sequence, from base 1 to 14 will be:
    #T-1, T, T+1, T+2, T+3, T+4, 3UP, 3UP-1, UP-1, UP, UP+1, UP+2, UP+3, UP+4
#%%
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#STEP 1: call in every image (i.e. for every base) & create masks for every base image
import os
import glob
from skimage.io import imread  
import tifffile as tif    

# Define the directory
directory = '~/dir/01_raw_images'

# Define the character mapping for base assignment
base_mapping = {'T-1': 'base01', 'T': 'base02', "T+1": 'base03', 'T+2': 'base04',
                'TruR2-1': 'base01', 'TruR2': 'base02', "TruR2+1": 'base03', 'TruR2+2': 'base04',
                'T+3': 'base05', 'T+4': 'base06', '3UP': 'base07', '3UP-1': 'base08',
                'UP-1': 'base09', 'UP': 'base10', 'UP+1': 'base11', 'UP+2': 'base12',
                'UP+3': 'base13', 'UP+4': 'base14'}

# Get a list of all .tif files in the directory
tif_files = glob.glob(os.path.join(directory, '*.tif'))

# Create a dictionary to store base images in the correct order
base_tif_images = {f'base{i:02d}': None for i in range(1, 15)}

# Loop through each file and assign it to the appropriate key in the dictionary
for file in tif_files:
    file_basename = os.path.basename(file)  # Get only the filename, not the full path
    for key, base_var in base_mapping.items():
        # Check for an exact match of the key as a substring in the filename
        if f'_{key}_' in file_basename or file_basename.startswith(f'{key}_') or file_basename.endswith(f'_{key}.tif') or file_basename == f'{key}.tif':
            # Read the tif image and store it in the dictionary in the correct order
            base_tif_images[base_var] = imread(file)
            print(f"Assigned {file_basename} to base_tif_images['{base_var}']")
            break

# Verify the ordering of the images in the dictionary
for base, image in base_tif_images.items():
    if image is not None:
        print(f"{base} is loaded.")
    else:
        print(f"{base} is missing or not matched to any file.")
        
# Now, `base_tif_images` is a dictionary; can access each image by its key (e.g. base_tif_images['base01'])

#%%
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#STEP 2A: convert all images to grayscale by averaging across the channels
import numpy as np                   #lib for numerial computations
import matplotlib.pyplot as plt      #lib for plotting / displaying images

# Convert the images to grayscale (b+w) by averaging across the channels
for base_var, image in base_tif_images.items():
    image = np.asarray(image)  # Ensure it's a NumPy array
    if image.ndim == 3:  # If it's (H, W, C), convert to grayscale
        base_bw_images[base_var] = np.mean(image, axis=-1)
    elif image.ndim == 2:  # Already grayscale
        base_bw_images[base_var] = image
    else:
        print(f"Warning: Image {base_var} has unexpected shape {image.shape}, skipping.")

#base_bw_images = {base_var: np.mean(image, axis=-1) for base_var, image in base_tif_images.items()}

#STEP 2B: display a test grayscale image to verify that it is functioning properly
plt.imshow(base_bw_images['base13'], cmap='gray')
plt.title('Grayscale Image')
plt.axis('off')
plt.show()

#%%
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#STEP 3: segment all images using CellPose

# First, check that GPU is activated and seen
!nvcc --version
!nvidia-smi

import os, shutil
import numpy as np
import matplotlib.pyplot as plt
from cellpose import core, utils, io, models, metrics
from glob import glob

use_GPU = core.use_gpu()
yn = ['NO', 'YES']
print(f'>>> GPU activated? {yn[use_GPU]}') #YES

# LOAD PACKAGES
from cellpose import models, denoise, io                   #import models module from cellpose package
from scipy.ndimage import center_of_mass, label            #import center_of_mass and label functions
from skimage.color import label2rgb
import matplotlib.pyplot as plt
from PIL import Image
import tifffile
import numpy as np

# DEFINE CELLPOSE MODEL
        # model_type="cyto3" or "nuclei", or other model
        # restore_type: "denoise_cyto3", "deblur_cyto3", "upsample_cyto3", "denoise_nuclei", "deblur_nuclei", "upsample_nuclei"
model = denoise.CellposeDenoiseModel(gpu=True, model_type="cyto3", restore_type="deblur_cyto3")

# Run CellPose segmentation (for all bases)
segmentations = {}

# Specify where masks should be saved
base_directory = '/dir'
mask_directory = f"{base_directory}/02_masks"
mask_directory = os.path.expanduser(mask_directory) # Add '~' to base_directory
os.makedirs(mask_directory, exist_ok=True) # Create mask_directory if it doesn't exist

# Specify the path for the CSV file
csv_file = f"{base_directory}/segmentation_results.csv"

# Create the CSV file and write the header (if it doesn't already exist)
import csv
if not os.path.exists(csv_file):
    with open(csv_file, mode='w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)
        writer.writerow(['filename', 'number_of_segmentations'])

# Run Cellpose w/ diameter=15pixels and 'cyto3' from model zoo
for base_var, gray_image in base_bw_images.items():
    
    # Run CellPose segmentation
    result = model.eval(gray_image, 
                        diameter=15,
                        flow_threshold=2.0,  # Increase flow threshold to separate adjacent circles
                        channels=[0, 0])     # Assuming grayscale image, [nucleus_channel, cytoplasm_channel]
    # Unpack the results from the model
    if len(result) == 4:
        masks, flows, styles, diams = result
    else:
        masks, flows, styles = result  # No diameters returned
    
    # Store the masks in the segmentations dictionary
    segmentations[base_var] = masks

    # Save the mask as a numpy array file
    np.save(f'{mask_directory}/{base_var}_masks.npy', masks)

    # Display mask in color to see distinct segmentations
    rgb_image = label2rgb(masks, bg_label=0)  # Convert to RGB with a black background
    
    # Print the number of segmentations
    num_segmentations = masks.max()
    print(f"Number of segmentations for {base_var}: {num_segmentations}")

    # Add number of segmentations to csv file
    with open(csv_file, mode='a', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([base_var, num_segmentations])

    # Rescale the image to 8-bit format for saving as PNG or TIFF
    rgb_image_8bit = (rgb_image * 255).astype(np.uint8)

    # Convert the RGB image to PIL Image format
    image_pil = Image.fromarray(rgb_image_8bit)
    
    # Save mask as tiff
    tifffile.imwrite(f'{mask_directory}/{base_var}_RGBmasks.tiff', rgb_image_8bit, photometric='rgb')

#%%
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# if need to load back in the segmentations:
import os
import numpy as np
import matplotlib.pyplot as plt

# Define the directory where the masks are saved
mask_directory = '~/dir/02_masks'

# Initialize the dictionary to store the loaded segmentations
segmentations = {}

# Loop through files in mask directory and load .npy files into dictionary
for file in sorted(os.listdir(mask_directory)):
    if file.endswith('_masks.npy'):
        base_var = file.split('_masks.npy')[0]  # Extract the base_var from the filename
        segmentations[base_var] = np.load(os.path.join(mask_directory, file))        

# Now `segmentations` dictionary contains the loaded masks
plt.imshow(segmentations['base03'], cmap='gray')
plt.show()

#%%
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#STEP 4: create binary masks for every image

#4A) first, isolate signal from only within segmented regions and create a binary mask

import cv2
import numpy as np                   #lib for numerial computations
import matplotlib.pyplot as plt      #lib for plotting / displaying images
from skimage import io
from skimage.filters import threshold_otsu, threshold_multiotsu

binary_masks = {}

for base_var in segmentations:
    if base_var in base_tif_images:
        # Get the segmentation mask and the original image for the current key
        segmentation_mask = segmentations[base_var]
        original_image = base_tif_images[base_var]
        
        # Apply filter out background, i.e. where signal is 0
        if original_image.ndim == 3:  # RGB image (H, W, C)
            int_filter = np.ma.masked_where(segmentations[base_var] == 0, original_image[:, :, 0])
        elif original_image.ndim == 2:  # Grayscale image (H, W)
            int_filter = np.ma.masked_where(segmentations[base_var] == 0, original_image)  # No need for `[:, :, 0]`
        else:
            raise ValueError(f"Unexpected image dimensions: {original_image.shape}")

        # Store binary mask in 'binary_masks' dictionary
        binary_masks[base_var] = segmentations[base_var]>0
    
    else:
        print(f"Warning: {base_var} not found in base_tif_images")

#4B) display a test binary 
plt.imshow(binary_masks['base14'], cmap='gray')
plt.title('Binary Mask Image')
plt.axis('off')
plt.show()

To save dictionary to a file
import pickle
with open('~/dir/binary_masks.pkl', 'wb') as file:
    pickle.dump(binary_masks, file)

# to load dictionary from file:
import pickle
with open('~/dir/binary_masks.pkl', 'rb') as file:
    binary_masks = pickle.load(file)

#%%
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#4C) QC: overlay pairs of mask images on each other (to check for images not having identical positions)
# plt.imshow(binary_masks['base04'], cmap='Reds', alpha=0.5)
# plt.imshow(binary_masks['base13'], cmap='Blues', alpha=0.3)
# plt.show()

# plt.imshow(binary_masks['base01'], cmap='Reds', alpha=0.5)
# plt.imshow(binary_masks['base12'], cmap='Blues', alpha=0.25)
# plt.show()

#%%
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# STEP 5A: Erode segmentations to create new binary masks (this will be used to make sure only 'pure' pixels are used for base calling)
import pickle
import numpy as np
import matplotlib.pyplot as plt
from skimage.segmentation import find_boundaries
from scipy.ndimage import distance_transform_edt

eroded_binary_masks = {}

for base_var in segmentations:
    if base_var in base_tif_images:
        # Get the segmentation mask and the original image for the current key
        segmentation = segmentations[base_var]
        original_image = base_tif_images[base_var]
        binary_mask = binary_masks[base_var]
        
        # Parameters: adjust the shrinkage value as desired (in pixels)
        shrink_px = 1.2

        # Step 1: Compute the boundary image of the segmentation
        # This creates a boolean image where True values indicate cell boundaries.
        boundaries = find_boundaries(segmentation, mode='inner')

        # Step 2: Compute the distance transform from the boundaries.
        # Here, we compute the EDT on the inverse of the boundaries so that 
        # each pixel value represents its distance to the nearest boundary.
        dt = distance_transform_edt(~boundaries)

        # Step 3: Threshold the distance transform to shrink the segmentation.
        # Only keep pixels that are at least `shrink_px` away from a boundary.
        mask_shrink = dt >= shrink_px

        # Step 4: Intersect the thresholded mask with the binary background mask.
        # Ensure the binary mask is boolean.
        binary_mask_bool = binary_mask.astype(bool)
        combined_mask = mask_shrink & binary_mask_bool

        # Step 5: Apply the combined mask to the segmentation.
        # Multiply the segmentation by the mask so that only the eroded regions retain their label.
        eroded_segmentation = segmentation * combined_mask
        
        # Create & store binary mask in 'eroded_masks' dictionary
        eroded_binary_masks[base_var] = eroded_segmentation > 0
        print(f'eroded mask created for {base_var}')

        # Visualize the result
        plt.figure(figsize=(8, 6))
        plt.imshow(binary_mask_bool[5500:6000, 5500:6000], cmap='grey', interpolation='nearest')
        plt.title(f'{base_var}: Bool Segmentation Mask')
        plt.colorbar(label='Cell Label')
        plt.show()
        
        plt.figure(figsize=(8, 6))
        plt.imshow(eroded_segmentation[5500:6000, 5500:6000], cmap='grey', interpolation='nearest')
        plt.title(f'{base_var}: Eroded Segmentation Mask')
        plt.colorbar(label='Cell Label')
        plt.show()

# Save dictionary to file    
# import pickle
# with open('~/dir/eroded_binary_masks.pkl', 'wb') as file:
#     pickle.dump(eroded_binary_masks, file)

# to load dictionary from file:
import pickle
with open('~/dir/eroded_binary_masks.pkl', 'rb') as file:
    eroded_binary_masks = pickle.load(file)

#%%
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# STEP 5B: Correct for ALL cross-channel bleedthroughs!!
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

def correct_bleedthroughs(
    original_image,
    binary_mask,
    fallback_slope_txr=0.5,  # fallback for Cy3→TxR (i.e. bleedthrough of Cy3 into TxR)
    fallback_slope_cy3=0.5,  # fallback for TxR→Cy3 (i.e. bleedthrough of TxR into Cy3)
    show_plots=True,
    base_name='Base'
):
    """
    Corrects bleedthrough as follows:
      - For TxR: corrected_txr = original TxR – (slope_cy3_to_txr × Cy3) – (slope_FAM_to_TxR × FAM)
      - For Cy3: corrected_cy3 = original Cy3 – (slope_txr_to_cy3 × TxR) – (slope_FAM_to_Cy3 × FAM)
      - For Cy5: corrected_cy5 = original Cy5 – (slope_txr_to_Cy5 × TxR)
      - For FAM: corrected_FAM = original FAM – (slope_Cy3_to_FAM × Cy3)
    
    (Density plots, histograms, and zoomed‐in views are produced for all channels.)
    
    Parameters
    ----------
    original_image : np.ndarray
        4-channel image with channels [Cy5, TxR, Cy3, FAM].
    binary_mask : np.ndarray
        Eroded mask (shape (H,W)) used for selecting pixels.
    fallback_slope_txr : float
        Fallback slope for Cy3→TxR if regression fails.
    fallback_slope_cy3 : float
        Fallback slope for TxR→Cy3 if regression fails.
    show_plots : bool
        If True, displays plots.
    base_name : str
        Name of the base image (used in figure titles).
    
    Returns
    -------
    corrected_image_uint16 : np.ndarray
        The final corrected 4-channel image scaled to uint16.
    slope_txr_to_cy3 : float
        Slope from Cy3→TxR (used in TxR correction).
    slope_cy3_to_txr : float
        Slope from TxR→Cy3 (used in Cy3 correction).
    """
    # --- Ensure the mask is Boolean ---
    binary_mask = binary_mask.astype(bool)
    
    # ------------------------------
    # 1) Extract channels
    # ------------------------------
    ch_cy5 = original_image[:, :, 0]  # Cy5
    ch_txr = original_image[:, :, 1]  # Texas Red (TR)
    ch_cy3 = original_image[:, :, 2]  # Cy3
    ch_FAM = original_image[:, :, 3]  # FAM

    # ------------------------------
    # 2) Apply the mask to each channel
    # ------------------------------
    masked_cy3 = ch_cy3[binary_mask]
    masked_txr = ch_txr[binary_mask]
    masked_cy5 = ch_cy5[binary_mask]
    masked_FAM = ch_FAM[binary_mask]

    # ------------------------------
    # 3) Calculate Slopes for cross-channel talk
    # ------------------------------
    valid_idx = (masked_cy3 > 0) & (masked_txr > 0)
    valid_cy3 = masked_cy3[valid_idx]
    valid_txr = masked_txr[valid_idx]
    
    if len(valid_cy3) < 2 or len(valid_txr) < 2:
        print(f"{base_name}: Not enough valid pixels. Using fallback slopes.")
        slope_cy3_to_txr = fallback_slope_txr
        slope_txr_to_cy3 = fallback_slope_cy3
    else:
        # --- Correction for bleedthrough from Cy3 → TxR ---
        cy3_75th = np.percentile(valid_cy3, 75)
        cy3_80th = np.percentile(valid_cy3, 80)
        cy3_85th = np.percentile(valid_cy3, 85)
        cy3_100th = np.percentile(valid_cy3, 100)
        txr_50th = np.percentile(valid_txr, 50)
        
        cy3_75th = np.percentile(valid_cy3, 75)
        subset_mask = valid_cy3 >= cy3_75th
        subset_cy3 = valid_cy3[subset_mask]
        subset_txr = valid_txr[subset_mask]
        
        if len(subset_cy3) < 2:
            print(f"{base_name}: Not enough points in Cy3→TxR subset. Using fallback slope.")
            slope_cy3_to_txr = fallback_slope_txr
        else:
            subset_ratio = subset_cy3 / subset_txr
            ratio_threshold = np.percentile(subset_ratio, 60)
            highlight_mask = subset_ratio >= ratio_threshold
            highlight_cy3 = subset_cy3[highlight_mask]
            highlight_txr = subset_txr[highlight_mask]
            if len(highlight_cy3) < 2:
                print(f"{base_name}: Not enough highlighted points for Cy3→TxR regression. Using fallback slope.")
                slope_cy3_to_txr = fallback_slope_txr
            else:
                coeffs = np.polyfit(highlight_cy3, highlight_txr, 1)
                slope_cy3_to_txr = coeffs[0]
        
        # --- Correction for bleedthrough from TxR → Cy3 ---
        txr_75th = np.percentile(valid_txr, 75)
        txr_80th = np.percentile(valid_txr, 80)
        txr_85th = np.percentile(valid_txr, 85)
        txr_100th = np.percentile(valid_txr, 100)
        cy3_50th = np.percentile(valid_cy3, 50)
        
        txr_75th = np.percentile(valid_txr, 75)
        subset_mask2 = valid_txr >= txr_75th
        subset_txr_2 = valid_txr[subset_mask2]
        subset_cy3_2 = valid_cy3[subset_mask2]
        
        if len(subset_txr_2) < 2:
            print(f"{base_name}: Not enough points in TxR→Cy3 subset. Using fallback slope.")
            slope_txr_to_cy3 = fallback_slope_cy3
        else:
            subset_ratio2 = subset_txr_2 / subset_cy3_2
            ratio_threshold2 = np.percentile(subset_ratio2, 60)
            highlight_mask2 = subset_ratio2 >= ratio_threshold2
            highlight_txr_2 = subset_txr_2[highlight_mask2]
            highlight_cy3_2 = subset_cy3_2[highlight_mask2]
            if len(highlight_txr_2) < 2:
                print(f"{base_name}: Not enough highlighted points for TxR→Cy3 regression. Using fallback slope.")
                slope_txr_to_cy3 = fallback_slope_cy3
            else:
                coeffs2 = np.polyfit(highlight_txr_2, highlight_cy3_2, 1)
                slope_txr_to_cy3 = coeffs2[0]
    
    # print(f"{base_name}: Computed slopes -> Cy3→TxR: {slope_cy3_to_txr:.4f}, TxR→Cy3: {slope_txr_to_cy3:.4f}")
    
    # ------------------------------
    # 4b) Compute additional slopes needed:
    #     (i) Slope from FAM into TxR (slope_FAM_to_TxR)
    #     (ii) Slope from FAM into Cy3 (slope_FAM_to_Cy3)
    #     (iii) Slope from TxR into Cy5 (slope_txr_to_Cy5)
    #     (iv) Slope from Cy3 into FAM (slope_Cy3_to_FAM)
    # ------------------------------
    # (i) Compute slope_FAM_to_TxR: Use pixels where FAM and TxR are >0.
    valid_idx_FAM_TxR = (masked_FAM > 0) & (masked_txr > 0)
    valid_FAM_for_TxR = masked_FAM[valid_idx_FAM_TxR]
    valid_txr_for_FAM = masked_txr[valid_idx_FAM_TxR]
    if len(valid_FAM_for_TxR) < 2:
        slope_FAM_to_TxR = 0.0
        print(f"{base_name}: Not enough valid pixels for FAM→TxR correction. Using fallback slope 0.0")
    else:
        fam_75th_TxR = np.percentile(valid_FAM_for_TxR, 75)
        subset_mask = valid_FAM_for_TxR >= fam_75th_TxR
        subset_FAM_TxR = valid_FAM_for_TxR[subset_mask]
        subset_txr_for_TxR = valid_txr_for_FAM[subset_mask]
        if len(subset_FAM_TxR) < 2:
            slope_FAM_to_TxR = 0.0
            print(f"{base_name}: Not enough points in FAM→TxR subset. Using fallback slope 0.0")
        else:
            ratio = subset_FAM_TxR / subset_txr_for_TxR
            ratio_thresh = np.percentile(ratio, 60)
            highlight_mask = ratio >= ratio_thresh
            highlight_FAM_TxR = subset_FAM_TxR[highlight_mask]
            highlight_txr_for_TxR = subset_txr_for_TxR[highlight_mask]
            if len(highlight_FAM_TxR) < 2:
                slope_FAM_to_TxR = 0.0
                print(f"{base_name}: Not enough highlighted points for FAM→TxR regression. Using fallback slope 0.0")
            else:
                coeffs = np.polyfit(highlight_FAM_TxR, highlight_txr_for_TxR, 1)
                slope_FAM_to_TxR = coeffs[0]
    
    # (ii) Compute slope_FAM_to_Cy3: Use pixels where FAM and Cy3 are >0.
    valid_idx_FAM_Cy3 = (masked_FAM > 0) & (masked_cy3 > 0)
    valid_FAM_for_Cy3 = masked_FAM[valid_idx_FAM_Cy3]
    valid_cy3_for_FAM = masked_cy3[valid_idx_FAM_Cy3]
    if len(valid_FAM_for_Cy3) < 2:
        slope_FAM_to_Cy3 = 0.0
        print(f"{base_name}: Not enough valid pixels for FAM→Cy3 correction. Using fallback slope 0.0")
    else:
        fam_75th_Cy3 = np.percentile(valid_FAM_for_Cy3, 75)
        subset_mask = valid_FAM_for_Cy3 >= fam_75th_Cy3
        subset_FAM_Cy3 = valid_FAM_for_Cy3[subset_mask]
        subset_cy3_for_FAM = valid_cy3_for_FAM[subset_mask]
        if len(subset_FAM_Cy3) < 2:
            slope_FAM_to_Cy3 = 0.0
            print(f"{base_name}: Not enough points in FAM→Cy3 subset. Using fallback slope 0.0")
        else:
            ratio = subset_FAM_Cy3 / subset_cy3_for_FAM
            ratio_thresh = np.percentile(ratio, 40)
            highlight_mask = ratio >= ratio_thresh
            highlight_FAM_Cy3 = subset_FAM_Cy3[highlight_mask]
            highlight_cy3_for_FAM = subset_cy3_for_FAM[highlight_mask]
            if len(highlight_FAM_Cy3) < 2:
                slope_FAM_to_Cy3 = 0.0
                print(f"{base_name}: Not enough highlighted points for FAM→Cy3 regression. Using fallback slope 0.0")
            else:
                coeffs = np.polyfit(highlight_FAM_Cy3, highlight_cy3_for_FAM, 1)
                slope_FAM_to_Cy3 = coeffs[0]
    
    # (iii) Compute slope_txr_to_Cy5: Use pixels where TxR and Cy5 are >0.
    valid_idx_Cy5_TxR = (masked_cy5 > 0) & (masked_txr > 0)
    valid_txr_for_Cy5 = masked_txr[valid_idx_Cy5_TxR]
    valid_cy5_for_corr = masked_cy5[valid_idx_Cy5_TxR]
    if len(valid_txr_for_Cy5) < 2:
        slope_txr_to_Cy5 = 0.0
        print(f"{base_name}: Not enough valid pixels for TxR→Cy5 correction. Using fallback slope 0.0")
    else:
        txr_75th_Cy5 = np.percentile(valid_txr_for_Cy5, 75)
        subset_mask = valid_txr_for_Cy5 >= txr_75th_Cy5
        subset_txr_Cy5 = valid_txr_for_Cy5[subset_mask]
        subset_cy5_for_corr = valid_cy5_for_corr[subset_mask]
        if len(subset_txr_Cy5) < 2:
            slope_txr_to_Cy5 = 0.0
            print(f"{base_name}: Not enough points in TxR→Cy5 subset. Using fallback slope 0.0")
        else:
            ratio = subset_txr_Cy5 / subset_cy5_for_corr
            ratio_thresh = np.percentile(ratio, 60)
            highlight_mask = ratio >= ratio_thresh
            highlight_txr_Cy5 = subset_txr_Cy5[highlight_mask]
            highlight_cy5_for_corr = subset_cy5_for_corr[highlight_mask]
            if len(highlight_txr_Cy5) < 2:
                slope_txr_to_Cy5 = 0.0
                print(f"{base_name}: Not enough highlighted points for TxR→Cy5 regression. Using fallback slope 0.0")
            else:
                coeffs = np.polyfit(highlight_txr_Cy5, highlight_cy5_for_corr, 1)
                slope_txr_to_Cy5 = coeffs[0]
    
    # (iv) Compute slope_Cy3_to_FAM: Use pixels where Cy3 and FAM are >0.
    valid_idx_Cy3_FAM = (masked_cy3 > 0) & (masked_FAM > 0)
    valid_cy3_for_FAM = masked_cy3[valid_idx_Cy3_FAM]
    valid_FAM_for_corr = masked_FAM[valid_idx_Cy3_FAM]
    if len(valid_cy3_for_FAM) < 2:
        slope_Cy3_to_FAM = 0.0
        print(f"{base_name}: Not enough valid pixels for Cy3→FAM correction. Using fallback slope 0.0")
    else:
        cy3_75th_FAM = np.percentile(valid_cy3_for_FAM, 75)
        subset_mask = valid_cy3_for_FAM >= cy3_75th_FAM
        subset_cy3_FAM = valid_cy3_for_FAM[subset_mask]
        subset_FAM_for_corr = valid_FAM_for_corr[subset_mask]
        if len(subset_cy3_FAM) < 2:
            slope_Cy3_to_FAM = 0.0
            print(f"{base_name}: Not enough points in Cy3→FAM subset. Using fallback slope 0.0")
        else:
            ratio = subset_cy3_FAM / subset_FAM_for_corr
            ratio_thresh = np.percentile(ratio, 60)
            highlight_mask = ratio >= ratio_thresh
            highlight_cy3_FAM = subset_cy3_FAM[highlight_mask]
            highlight_FAM_for_corr = subset_FAM_for_corr[highlight_mask]
            if len(highlight_cy3_FAM) < 2:
                slope_Cy3_to_FAM = 0.0
                print(f"{base_name}: Not enough highlighted points for Cy3→FAM regression. Using fallback slope 0.0")
            else:
                coeffs = np.polyfit(highlight_cy3_FAM, highlight_FAM_for_corr, 1)
                slope_Cy3_to_FAM = coeffs[0]
    
    # ------------------------------
    # 4) Density Plots + Slopes
    # ------------------------------   
    if show_plots:
        # Density plot 1: Cy3→TxR correction
        plt.figure(figsize=(8,6))
        hb = plt.hexbin(valid_cy3, valid_txr, gridsize=500, cmap='viridis', mincnt=1)
        plt.colorbar(hb, label='Counts')
        plt.axvline(x=cy3_75th, color='black', linestyle='--', linewidth=0.8,
                    label=f'TxR 75th ({cy3_75th:.2f})')
        if len(highlight_cy3) >= 2:
            x_fit = np.linspace(highlight_cy3.min(), highlight_cy3.max(), 100)
            intercept_val = np.polyfit(highlight_cy3, highlight_txr, 1)[1]
            y_fit = slope_cy3_to_txr * x_fit + intercept_val
            plt.plot(x_fit, y_fit, color='green', linewidth=1,
                     label=f'Regression (slope={slope_cy3_to_txr:.4f})')
        plt.xlabel('Cy3 Intensity')
        plt.ylabel('TxR Intensity')
        plt.title(f'{base_name}: Density Plot (Cy3 → TxRed Correction)')
        plt.legend(loc='lower right')
        plt.show()
         
        # Density plot 2: TxR→Cy3 correction
        plt.figure(figsize=(8,6))
        hb = plt.hexbin(valid_txr, valid_cy3, gridsize=500, cmap='viridis', mincnt=1)
        plt.colorbar(hb, label='Counts')
        plt.axvline(x=txr_75th, color='black', linestyle='--', linewidth=0.8,
                    label=f'TxRed 75th ({txr_75th:.2f})')
        if len(highlight_txr_2) >= 2:
            x_fit = np.linspace(highlight_txr_2.min(), highlight_txr_2.max(), 100)
            intercept_val = np.polyfit(highlight_txr_2, highlight_cy3_2, 1)[1]
            y_fit = slope_txr_to_cy3 * x_fit + intercept_val
            plt.plot(x_fit, y_fit, color='green', linewidth=1,
                      label=f'Regression (slope={slope_txr_to_cy3:.4f})')
        plt.xlabel('TxR Intensity')
        plt.ylabel('Cy3 Intensity')
        plt.title(f'{base_name}: Density Plot (TxR → Cy3 Correction)')
        plt.legend(loc='lower right')
        plt.show()
        
        # # Density plot 3: FAM→TxR correction
        # plt.figure(figsize=(8,6))
        # hb = plt.hexbin(valid_FAM_for_TxR, valid_txr_for_FAM, gridsize=500, cmap='viridis', mincnt=1)
        # plt.colorbar(hb, label='Counts')
        # plt.axvline(x=fam_75th_TxR, color='black', linestyle='--', linewidth=0.8,
        #             label=f'FAM 75th ({fam_75th_TxR:.2f})')
        # if len(highlight_FAM_TxR) >= 2:
        #     x_fit = np.linspace(highlight_FAM_TxR.min(), highlight_FAM_TxR.max(), 100)
        #     intercept_val = np.polyfit(highlight_FAM_TxR, highlight_txr_for_TxR, 1)[1]
        #     y_fit = slope_FAM_to_TxR * x_fit + intercept_val
        #     plt.plot(x_fit, y_fit, color='green', linewidth=1,
        #              label=f'Regression (slope={slope_FAM_to_TxR:.4f})')
        # plt.xlabel('FAM Intensity')
        # plt.ylabel('TxR Intensity')
        # plt.title(f'{base_name}: Density Plot (FAM → TxRed Correction)')
        # plt.legend(loc='lower right')
        # plt.show()
         
        # Density plot 4: FAM→Cy3 correction
        plt.figure(figsize=(8,6))
        hb = plt.hexbin(valid_FAM_for_Cy3, valid_cy3_for_FAM, gridsize=500, cmap='viridis', mincnt=1)
        plt.colorbar(hb, label='Counts')
        plt.axvline(x=fam_75th_Cy3, color='black', linestyle='--', linewidth=0.8,
                    label=f'FAM 75th ({fam_75th_Cy3:.2f})')
        if len(highlight_FAM_Cy3) >= 2:
            x_fit = np.linspace(highlight_FAM_Cy3.min(), highlight_FAM_Cy3.max(), 100)
            intercept_val = np.polyfit(highlight_FAM_Cy3, highlight_cy3_for_FAM, 1)[1]
            y_fit = slope_FAM_to_Cy3 * x_fit + intercept_val
            plt.plot(x_fit, y_fit, color='green', linewidth=1,
                      label=f'Regression (slope={slope_FAM_to_Cy3:.4f})')
        plt.xlabel('FAM Intensity')
        plt.ylabel('Cy3 Intensity')
        plt.title(f'{base_name}: Density Plot (FAM → Cy3 Correction)')
        plt.legend(loc='lower right')
        plt.show()

        # Density plot 5: TxR→Cy5 correction
        plt.figure(figsize=(8,6))
        hb = plt.hexbin(valid_txr_for_Cy5, valid_cy5_for_corr, gridsize=500, cmap='viridis', mincnt=1)
        plt.colorbar(hb, label='Counts')
        plt.axvline(x=txr_75th_Cy5, color='black', linestyle='--', linewidth=0.8,
                    label=f'TxR 75th ({txr_75th_Cy5:.2f})')
        if len(highlight_txr_Cy5) >= 2:
            x_fit = np.linspace(highlight_txr_Cy5.min(), highlight_txr_Cy5.max(), 100)
            intercept_val = np.polyfit(highlight_txr_Cy5, highlight_cy5_for_corr, 1)[1]
            y_fit = slope_txr_to_Cy5 * x_fit + intercept_val
            plt.plot(x_fit, y_fit, color='green', linewidth=1,
                     label=f'Regression (slope={slope_txr_to_Cy5:.4f})')
        plt.xlabel('TxR Intensity')
        plt.ylabel('Cy5 Intensity')
        plt.title(f'{base_name}: Density Plot (TxR → Cy5 Correction)')
        plt.legend(loc='lower right')
        plt.show()
         
        # Density plot 6: Cy3→FAM correction
        plt.figure(figsize=(8,6))
        hb = plt.hexbin(valid_cy3_for_FAM, valid_FAM_for_corr, gridsize=500, cmap='viridis', mincnt=1)
        plt.colorbar(hb, label='Counts')
        plt.axvline(x=cy3_75th_FAM, color='black', linestyle='--', linewidth=0.8,
                    label=f'Cy3 75th ({cy3_75th_FAM:.2f})')
        if len(highlight_cy3_FAM) >= 2:
            x_fit = np.linspace(highlight_cy3_FAM.min(), highlight_cy3_FAM.max(), 100)
            intercept_val = np.polyfit(highlight_cy3_FAM, highlight_FAM_for_corr, 1)[1]
            y_fit = slope_Cy3_to_FAM * x_fit + intercept_val
            plt.plot(x_fit, y_fit, color='green', linewidth=1,
                      label=f'Regression (slope={slope_Cy3_to_FAM:.4f})')
        plt.xlabel('Cy3 Intensity')
        plt.ylabel('FAM Intensity')
        plt.title(f'{base_name}: Density Plot (Cy3 → FAM Correction)')
        plt.legend(loc='lower right')
        plt.show()
    
        # ------------------------------
    # 5) Apply the corrections to the image channels
    # ------------------------------
    # For TxR: subtract bleedthrough from Cy3 and from FAM.
    corrected_txr = ch_txr - slope_cy3_to_txr * ch_cy3 #- slope_FAM_to_TxR * ch_FAM
    corrected_txr = np.clip(corrected_txr, 0, None)
    # For Cy3: subtract bleedthrough from TxR and from FAM.
    corrected_cy3 = ch_cy3 - slope_txr_to_cy3 * ch_txr - slope_FAM_to_Cy3 * ch_FAM
    corrected_cy3 = np.clip(corrected_cy3, 0, None)
    # For Cy5: subtract bleedthrough from TxR only.
    corrected_cy5 = ch_cy5 - slope_txr_to_Cy5 * ch_txr
    corrected_cy5 = np.clip(corrected_cy5, 0, None)
    # For FAM: subtract bleedthrough from Cy3 only.
    corrected_FAM = ch_FAM - slope_Cy3_to_FAM * ch_cy3
    corrected_FAM = np.clip(corrected_FAM, 0, None)
    
    # ------------------------------
    # 6) Final reconstruction corrected 4-channel image + scale to uint16
    # ------------------------------
    corrected_image = np.dstack([corrected_cy5, corrected_txr, corrected_cy3, corrected_FAM])
    min_val, max_val = corrected_image.min(), corrected_image.max()
    if max_val == min_val:
        corrected_image_uint16 = np.zeros_like(corrected_image, dtype=np.uint16)
    else:
        corrected_image_uint16 = ((65535 * (corrected_image - min_val) / (max_val - min_val)).astype(np.uint16))
    
    # ------------------------------
    # 7) Histograms for all channels before/after corrections
    # ------------------------------
    if show_plots:
        # Histograms for Cy5:
        xlim_before = np.percentile(masked_cy5, 99.9)
        xlim_after = np.percentile(corrected_cy5[binary_mask], 99.5)
        p70_cy5_before = np.percentile(masked_cy5, 70)
        p75_cy5_before = np.percentile(masked_cy5, 75)
        p80_cy5_before = np.percentile(masked_cy5, 80)
        p70_cy5_after = np.percentile(corrected_cy5[binary_mask], 70)
        p75_cy5_after = np.percentile(corrected_cy5[binary_mask], 75)
        p80_cy5_after = np.percentile(corrected_cy5[binary_mask], 80)
        
        plt.figure(figsize=(12,5))
        plt.subplot(1,2,1)
        plt.hist(masked_cy5, bins=500, color='green', alpha=0.3)
        plt.axvline(p70_cy5_before, color='red', linestyle=':', linewidth=1, label='70th Percentile')
        plt.axvline(p75_cy5_before, color='black', linestyle=':', linewidth=1, label='75th Percentile')
        plt.axvline(p80_cy5_before, color='blue', linestyle=':', linewidth=1, label='80th Percentile')
        #plt.xlim(100, xlim_before)
        plt.title(f'{base_name}: Masked Cy5 (Before Correction)')
        plt.gca().yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,.0f}'))
        plt.legend()
        
        plt.subplot(1,2,2)
        plt.hist(corrected_cy5[binary_mask], bins=500, color='green', alpha=0.3)
        plt.axvline(p70_cy5_after, color='red', linestyle=':', linewidth=1, label='70th Percentile')
        plt.axvline(p75_cy5_after, color='black', linestyle=':', linewidth=1, label='75th Percentile')
        plt.axvline(p80_cy5_after, color='blue', linestyle=':', linewidth=1, label='80th Percentile')
        #plt.xlim(100, xlim_after)
        plt.title(f'{base_name}: Corrected Cy5 (After Correction)')
        plt.gca().yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,.0f}'))
        plt.legend()
        plt.show()
        
        # Histograms for Cy3:
        xlim_before = np.percentile(masked_cy3, 99.9)
        xlim_after = np.percentile(corrected_cy3[binary_mask], 99.5)
        p70_cy3_before = np.percentile(masked_cy3, 70)
        p75_cy3_before = np.percentile(masked_cy3, 75)
        p80_cy3_before = np.percentile(masked_cy3, 80)
        p70_cy3_after = np.percentile(corrected_cy3[binary_mask], 70)
        p75_cy3_after = np.percentile(corrected_cy3[binary_mask], 75)
        p80_cy3_after = np.percentile(corrected_cy3[binary_mask], 80)
        
        plt.figure(figsize=(12,5))
        plt.subplot(1,2,1)
        plt.hist(masked_cy3, bins=500, color='orange', alpha=0.3)
        plt.axvline(p70_cy3_before, color='red', linestyle=':', linewidth=1, label='70th Percentile')
        plt.axvline(p75_cy3_before, color='black', linestyle=':', linewidth=1, label='75th Percentile')
        plt.axvline(p80_cy3_before, color='blue', linestyle=':', linewidth=1, label='80th Percentile')
        #plt.xlim(100, xlim_before)
        plt.title(f'{base_name}: Masked Cy3 (Before Correction)')
        plt.gca().yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,.0f}'))
        plt.legend()
        
        plt.subplot(1,2,2)
        plt.hist(corrected_cy3[binary_mask], bins=500, color='orange', alpha=0.3)
        plt.axvline(p70_cy3_after, color='red', linestyle=':', linewidth=1, label='70th Percentile')
        plt.axvline(p75_cy3_after, color='black', linestyle=':', linewidth=1, label='75th Percentile')
        plt.axvline(p80_cy3_after, color='blue', linestyle=':', linewidth=1, label='80th Percentile')
        #plt.xlim(200, xlim_after)
        plt.ylim(0, 350000)
        plt.title(f'{base_name}: Corrected Cy3 (After Correction)')
        plt.gca().yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,.0f}'))
        plt.legend()
        plt.show()
        
        # Histograms for txr:
        xlim_before = np.percentile(masked_txr, 99.9)
        xlim_after = np.percentile(corrected_txr[binary_mask], 99.5)
        p70_txr_before = np.percentile(masked_txr, 70)
        p75_txr_before = np.percentile(masked_txr, 75)
        p80_txr_before = np.percentile(masked_txr, 80)
        p70_txr_after = np.percentile(corrected_txr[binary_mask], 70)
        p75_txr_after = np.percentile(corrected_txr[binary_mask], 75)
        p80_txr_after = np.percentile(corrected_txr[binary_mask], 80)
        
        plt.figure(figsize=(12,5))
        plt.subplot(1,2,1)
        plt.hist(masked_txr, bins=500, color='purple', alpha=0.3)
        plt.axvline(p70_txr_before, color='red', linestyle=':', linewidth=1, label='70th Percentile')
        plt.axvline(p75_txr_before, color='black', linestyle=':', linewidth=1, label='75th Percentile')
        plt.axvline(p80_txr_before, color='blue', linestyle=':', linewidth=1, label='80th Percentile')
        #plt.xlim(100, xlim_before)
        plt.title(f'{base_name}: Masked TxR (Before Correction)')
        plt.gca().yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,.0f}'))
        plt.legend()
        
        plt.subplot(1,2,2)
        plt.hist(corrected_txr[binary_mask], bins=500, color='purple', alpha=0.3)
        plt.axvline(p70_txr_after, color='red', linestyle=':', linewidth=1, label='70th Percentile')
        plt.axvline(p75_txr_after, color='black', linestyle=':', linewidth=1, label='75th Percentile')
        plt.axvline(p80_txr_after, color='blue', linestyle=':', linewidth=1, label='80th Percentile')
        #plt.xlim(100, xlim_after)
        #plt.ylim(0,600000)
        plt.title(f'{base_name}: Corrected TxR (After Correction)')
        plt.gca().yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,.0f}'))
        plt.legend()
        plt.show()
        
        # Histograms for FAM:
        xlim_before = np.percentile(masked_FAM, 99.9)
        xlim_after = np.percentile(corrected_FAM[binary_mask], 99.5)
        p70_FAM_before = np.percentile(masked_FAM, 70)
        p75_FAM_before = np.percentile(masked_FAM, 75)
        p80_FAM_before = np.percentile(masked_FAM, 80)
        p70_FAM_after = np.percentile(corrected_FAM[binary_mask], 70)
        p75_FAM_after = np.percentile(corrected_FAM[binary_mask], 75)
        p80_FAM_after = np.percentile(corrected_FAM[binary_mask], 80)
        
        plt.figure(figsize=(12,5))
        plt.subplot(1,2,1)
        plt.hist(masked_FAM, bins=500, color='brown', alpha=0.3)
        plt.axvline(p70_FAM_before, color='red', linestyle=':', linewidth=1, label='70th Percentile')
        plt.axvline(p75_FAM_before, color='black', linestyle=':', linewidth=1, label='75th Percentile')
        plt.axvline(p80_FAM_before, color='blue', linestyle=':', linewidth=1, label='80th Percentile')
        #plt.xlim(100, xlim_before)
        plt.title(f'{base_name}: Masked FAM (Before Correction)')
        plt.gca().yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,.0f}'))
        plt.legend()
        
        plt.subplot(1,2,2)
        plt.hist(corrected_FAM[binary_mask], bins=500, color='brown', alpha=0.3)
        plt.axvline(p70_FAM_after, color='red', linestyle=':', linewidth=1, label='70th Percentile')
        plt.axvline(p75_FAM_after, color='black', linestyle=':', linewidth=1, label='75th Percentile')
        plt.axvline(p80_FAM_after, color='blue', linestyle=':', linewidth=1, label='80th Percentile')
        #plt.xlim(100, xlim_after)
        plt.title(f'{base_name}: Corrected FAM (After Correction)')
        plt.gca().yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,.0f}'))
        plt.legend()
        plt.show()
    
    
    return corrected_image_uint16, slope_cy3_to_txr, slope_FAM_to_TxR, slope_txr_to_cy3, slope_FAM_to_Cy3, slope_txr_to_Cy5, slope_Cy3_to_FAM

# ===============================================
# Loop over all bases / apply corrections / save corrected images to a dictionary
# ===============================================
corrected_base_images = {}
slopes_by_base = {}  # Stores tuple: (slope_txr_to_cy3, slope_cy3_to_txr) for each base

for base_name, image in base_tif_images.items():
    if base_name == 'align':  
        continue  # Skip non-image files
    else: 
        base_image = image
        eroded_binary_mask = eroded_binary_masks[base_name]
        print(f'{base_name}')
        


        corrected_img, slope_cy3_to_txr, slope_FAM_to_TxR, slope_txr_to_cy3, slope_FAM_to_Cy3, slope_txr_to_Cy5, slope_Cy3_to_FAM = correct_bleedthroughs(
            original_image=base_image,
            binary_mask=eroded_binary_mask,
            fallback_slope_txr=0.5,
            fallback_slope_cy3=0.5,
            show_plots=True,
            base_name=base_name
        )

        corrected_base_images[base_name] = corrected_img
        slopes_by_base[base_name] = (slope_cy3_to_txr, slope_txr_to_cy3, slope_FAM_to_TxR, slope_FAM_to_Cy3, slope_txr_to_Cy5, slope_Cy3_to_FAM)

        # ------------------------------------------------------------------
        # Zoomed in views for all channels:
        # TxR:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        axes[0].imshow(base_image[3000:3500, 3500:4000, 1], cmap='gray')
        axes[0].set_title(f'{base_name} - Original TxR (500x500)')
        axes[1].imshow(corrected_img[3000:3500, 3500:4000, 1], cmap='gray')
        axes[1].set_title(f'{base_name} - Corrected TxR (500x500)')
        plt.show()
        # Cy3:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        axes[0].imshow(base_image[3000:3500, 3500:4000, 2], cmap='gray')
        axes[0].set_title(f'{base_name} - Original Cy3 (500x500)')
        axes[1].imshow(corrected_img[3000:3500, 3500:4000, 2], cmap='gray')
        axes[1].set_title(f'{base_name} - Corrected Cy3 (500x500)')
        plt.show()
        # Cy5:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        axes[0].imshow(base_image[3000:3500, 3500:4000, 0], cmap='gray')
        axes[0].set_title(f'{base_name} - Original Cy5 (500x500)')
        axes[1].imshow(corrected_img[3000:3500, 3500:4000, 0], cmap='gray')
        axes[1].set_title(f'{base_name} - Corrected Cy5 (500x500)')
        plt.show()
        # FAM:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        axes[0].imshow(base_image[3000:3500, 3500:4000, 3], cmap='gray')
        axes[0].set_title(f'{base_name} - Original FAM (500x500)')
        axes[1].imshow(corrected_img[3000:3500, 3500:4000, 3], cmap='gray')
        axes[1].set_title(f'{base_name} - Corrected FAM (500x500)')
        plt.show()
        
        # ------------------------------------------------------------------
        # Print all slopes used in bleedthrough corrections:
        print(f"{base_name}: Corrected bleedthroughs with slopes -> Cy3→TxR: {slope_cy3_to_txr:.4f}, TxR→Cy3: {slope_txr_to_cy3:.4f}")
        print(f"{base_name}: Corrected bleedthroughs with slopes -> FAM→TxR: {slope_FAM_to_TxR:.4f}, FAM→Cy3: {slope_FAM_to_Cy3:.4f}")
        print(f"{base_name}: Corrected bleedthroughs with slopes -> TxR→Cy5: {slope_txr_to_Cy5:.4f}, Cy3→FAM: {slope_Cy3_to_FAM:.4f}")

    # else:
    #     print(f"No mask for {base_name}. Skipping.")


# Save corrected base images to a file
# import pickle
# with open('~/dir/corrected_base_images_ALLchannels.pkl', 'wb') as file:
#     pickle.dump(corrected_base_images, file)

# to load dictionary from file:
import pickle
with open('~/dir/corrected_base_images_ALLchannels.pkl', 'rb') as file:
    corrected_base_images = pickle.load(file)

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#%% STEP 5C: Save all bleedthrough correction outputs into a PDF
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.backends.backend_pdf import PdfPages

# Set working directory and output PDF file
directory = '~/dir'
pdf_path = os.path.join(directory, 'bleedthrough_corrections_outputs.pdf')

# Toggle saving mode: if True, figures are saved to the PDF instead of being shown interactively.
SAVE_TO_PDF = True

def save_or_show(pdf, fig):
    """Save the current figure to the PDF if enabled; otherwise, show it."""
    if SAVE_TO_PDF and pdf is not None:
        pdf.savefig(fig)
        plt.close(fig)
    else:
        plt.show()

def correct_bleedthroughs(
    original_image,
    binary_mask,
    fallback_slope_txr=0.5,  # fallback for Cy3→TxR
    fallback_slope_cy3=0.5,  # fallback for TxR→Cy3
    show_plots=True,
    base_name='Base',
    pdf=None  # PdfPages object; if provided, figures are saved here.
):
    """
    Corrects bleedthrough and produces a series of density plots and histograms.
    Instead of calling plt.show(), this version calls save_or_show() so that figures are either
    saved to the provided PdfPages object or displayed.
    """
    # 1) Ensure mask is boolean
    binary_mask = binary_mask.astype(bool)
    
    # 2) Extract channels
    ch_cy5 = original_image[:, :, 0]  # Cy5
    ch_txr = original_image[:, :, 1]  # TxR
    ch_cy3 = original_image[:, :, 2]  # Cy3
    ch_FAM = original_image[:, :, 3]  # FAM
    
    # 3) Apply the mask to each channel
    masked_cy3 = ch_cy3[binary_mask]
    masked_txr = ch_txr[binary_mask]
    masked_cy5 = ch_cy5[binary_mask]
    masked_FAM = ch_FAM[binary_mask]
    
    # 4) Calculate slopes for cross-channel bleedthrough
    valid_idx = (masked_cy3 > 0) & (masked_txr > 0)
    valid_cy3 = masked_cy3[valid_idx]
    valid_txr = masked_txr[valid_idx]
    
    if len(valid_cy3) < 2 or len(valid_txr) < 2:
        print(f"{base_name}: Not enough valid pixels. Using fallback slopes.")
        slope_cy3_to_txr = fallback_slope_txr
        slope_txr_to_cy3 = fallback_slope_cy3
    else:
        # --- Correction for bleedthrough from Cy3 → TxR ---
        cy3_75th = np.percentile(valid_cy3, 75)
        subset_mask = valid_cy3 >= cy3_75th
        subset_cy3 = valid_cy3[subset_mask]
        subset_txr = valid_txr[subset_mask]
        if len(subset_cy3) < 2:
            print(f"{base_name}: Not enough points in Cy3→TxR subset. Using fallback slope.")
            slope_cy3_to_txr = fallback_slope_txr
        else:
            subset_ratio = subset_cy3 / subset_txr
            ratio_threshold = np.percentile(subset_ratio, 60)
            highlight_mask = subset_ratio >= ratio_threshold
            highlight_cy3 = subset_cy3[highlight_mask]
            highlight_txr = subset_txr[highlight_mask]
            if len(highlight_cy3) < 2:
                print(f"{base_name}: Not enough highlighted points for Cy3→TxR regression. Using fallback slope.")
                slope_cy3_to_txr = fallback_slope_txr
            else:
                coeffs = np.polyfit(highlight_cy3, highlight_txr, 1)
                slope_cy3_to_txr = coeffs[0]
        
        # --- Correction for bleedthrough from TxR → Cy3 ---
        txr_75th = np.percentile(valid_txr, 75)
        subset_mask2 = valid_txr >= txr_75th
        subset_txr_2 = valid_txr[subset_mask2]
        subset_cy3_2 = valid_cy3[subset_mask2]
        if len(subset_txr_2) < 2:
            print(f"{base_name}: Not enough points in TxR→Cy3 subset. Using fallback slope.")
            slope_txr_to_cy3 = fallback_slope_cy3
        else:
            subset_ratio2 = subset_txr_2 / subset_cy3_2
            ratio_threshold2 = np.percentile(subset_ratio2, 60)
            highlight_mask2 = subset_ratio2 >= ratio_threshold2
            highlight_txr_2 = subset_txr_2[highlight_mask2]
            highlight_cy3_2 = subset_cy3_2[highlight_mask2]
            if len(highlight_txr_2) < 2:
                print(f"{base_name}: Not enough highlighted points for TxR→Cy3 regression. Using fallback slope.")
                slope_txr_to_cy3 = fallback_slope_cy3
            else:
                coeffs2 = np.polyfit(highlight_txr_2, highlight_cy3_2, 1)
                slope_txr_to_cy3 = coeffs2[0]
    
    # 5) Compute additional slopes:
    # (i) Slope from FAM into TxR
    valid_idx_FAM_TxR = (masked_FAM > 0) & (masked_txr > 0)
    valid_FAM_for_TxR = masked_FAM[valid_idx_FAM_TxR]
    valid_txr_for_TxR = masked_txr[valid_idx_FAM_TxR]
    if len(valid_FAM_for_TxR) < 2:
        slope_FAM_to_TxR = 0.0
        print(f"{base_name}: Not enough valid pixels for FAM→TxR correction. Using fallback slope 0.0")
    else:
        fam_75th_TxR = np.percentile(valid_FAM_for_TxR, 75)
        subset_mask = valid_FAM_for_TxR >= fam_75th_TxR
        subset_FAM_TxR = valid_FAM_for_TxR[subset_mask]
        subset_txr_for_TxR = valid_txr_for_TxR[subset_mask]
        if len(subset_FAM_TxR) < 2:
            slope_FAM_to_TxR = 0.0
            print(f"{base_name}: Not enough points in FAM→TxR subset. Using fallback slope 0.0")
        else:
            ratio = subset_FAM_TxR / subset_txr_for_TxR
            ratio_thresh = np.percentile(ratio, 60)
            highlight_mask = ratio >= ratio_thresh
            highlight_FAM_TxR = subset_FAM_TxR[highlight_mask]
            highlight_txr_for_TxR = subset_txr_for_TxR[highlight_mask]
            if len(highlight_FAM_TxR) < 2:
                slope_FAM_to_TxR = 0.0
                print(f"{base_name}: Not enough highlighted points for FAM→TxR regression. Using fallback slope 0.0")
            else:
                coeffs = np.polyfit(highlight_FAM_TxR, highlight_txr_for_TxR, 1)
                slope_FAM_to_TxR = coeffs[0]
    
    # (ii) Slope from FAM into Cy3
    valid_idx_FAM_Cy3 = (masked_FAM > 0) & (masked_cy3 > 0)
    valid_FAM_for_Cy3 = masked_FAM[valid_idx_FAM_Cy3]
    valid_cy3_for_FAM = masked_cy3[valid_idx_FAM_Cy3]
    if len(valid_FAM_for_Cy3) < 2:
        slope_FAM_to_Cy3 = 0.0
        print(f"{base_name}: Not enough valid pixels for FAM→Cy3 correction. Using fallback slope 0.0")
    else:
        fam_75th_Cy3 = np.percentile(valid_FAM_for_Cy3, 75)
        subset_mask = valid_FAM_for_Cy3 >= fam_75th_Cy3
        subset_FAM_Cy3 = valid_FAM_for_Cy3[subset_mask]
        subset_cy3_for_FAM = valid_cy3_for_FAM[subset_mask]
        if len(subset_FAM_Cy3) < 2:
            slope_FAM_to_Cy3 = 0.0
            print(f"{base_name}: Not enough points in FAM→Cy3 subset. Using fallback slope 0.0")
        else:
            ratio = subset_FAM_Cy3 / subset_cy3_for_FAM
            ratio_thresh = np.percentile(ratio, 40)
            highlight_mask = ratio >= ratio_thresh
            highlight_FAM_Cy3 = subset_FAM_Cy3[highlight_mask]
            highlight_cy3_for_FAM = subset_cy3_for_FAM[highlight_mask]
            if len(highlight_FAM_Cy3) < 2:
                slope_FAM_to_Cy3 = 0.0
                print(f"{base_name}: Not enough highlighted points for FAM→Cy3 regression. Using fallback slope 0.0")
            else:
                coeffs = np.polyfit(highlight_FAM_Cy3, highlight_cy3_for_FAM, 1)
                slope_FAM_to_Cy3 = coeffs[0]
    
    # (iii) Slope from TxR into Cy5
    valid_idx_Cy5_TxR = (masked_cy5 > 0) & (masked_txr > 0)
    valid_txr_for_Cy5 = masked_txr[valid_idx_Cy5_TxR]
    valid_cy5_for_corr = masked_cy5[valid_idx_Cy5_TxR]
    if len(valid_txr_for_Cy5) < 2:
        slope_txr_to_Cy5 = 0.0
        print(f"{base_name}: Not enough valid pixels for TxR→Cy5 correction. Using fallback slope 0.0")
    else:
        txr_75th_Cy5 = np.percentile(valid_txr_for_Cy5, 75)
        subset_mask = valid_txr_for_Cy5 >= txr_75th_Cy5
        subset_txr_Cy5 = valid_txr_for_Cy5[subset_mask]
        subset_cy5_for_corr = valid_cy5_for_corr[subset_mask]
        if len(subset_txr_Cy5) < 2:
            slope_txr_to_Cy5 = 0.0
            print(f"{base_name}: Not enough points in TxR→Cy5 subset. Using fallback slope 0.0")
        else:
            ratio = subset_txr_Cy5 / subset_cy5_for_corr
            ratio_thresh = np.percentile(ratio, 60)
            highlight_mask = ratio >= ratio_thresh
            highlight_txr_Cy5 = subset_txr_Cy5[highlight_mask]
            highlight_cy5_for_corr = subset_cy5_for_corr[highlight_mask]
            if len(highlight_txr_Cy5) < 2:
                slope_txr_to_Cy5 = 0.0
                print(f"{base_name}: Not enough highlighted points for TxR→Cy5 regression. Using fallback slope 0.0")
            else:
                coeffs = np.polyfit(highlight_txr_Cy5, highlight_cy5_for_corr, 1)
                slope_txr_to_Cy5 = coeffs[0]
    
    # (iv) Slope from Cy3 into FAM
    valid_idx_Cy3_FAM = (masked_cy3 > 0) & (masked_FAM > 0)
    valid_cy3_for_FAM = masked_cy3[valid_idx_Cy3_FAM]
    valid_FAM_for_corr = masked_FAM[valid_idx_Cy3_FAM]
    if len(valid_cy3_for_FAM) < 2:
        slope_Cy3_to_FAM = 0.0
        print(f"{base_name}: Not enough valid pixels for Cy3→FAM correction. Using fallback slope 0.0")
    else:
        cy3_75th_FAM = np.percentile(valid_cy3_for_FAM, 75)
        subset_mask = valid_cy3_for_FAM >= cy3_75th_FAM
        subset_cy3_FAM = valid_cy3_for_FAM[subset_mask]
        subset_FAM_for_corr = valid_FAM_for_corr[subset_mask]
        if len(subset_cy3_FAM) < 2:
            slope_Cy3_to_FAM = 0.0
            print(f"{base_name}: Not enough points in Cy3→FAM subset. Using fallback slope 0.0")
        else:
            ratio = subset_cy3_FAM / subset_FAM_for_corr
            ratio_thresh = np.percentile(ratio, 60)
            highlight_mask = ratio >= ratio_thresh
            highlight_cy3_FAM = subset_cy3_FAM[highlight_mask]
            highlight_FAM_for_corr = subset_FAM_for_corr[highlight_mask]
            if len(highlight_cy3_FAM) < 2:
                slope_Cy3_to_FAM = 0.0
                print(f"{base_name}: Not enough highlighted points for Cy3→FAM regression. Using fallback slope 0.0")
            else:
                coeffs = np.polyfit(highlight_cy3_FAM, highlight_FAM_for_corr, 1)
                slope_Cy3_to_FAM = coeffs[0]
    
    # 6) Generate density plots and regression overlays
    if show_plots:
        # Density Plot 1: Cy3 → TxR
        fig = plt.figure(figsize=(8,6))
        hb = plt.hexbin(valid_cy3, valid_txr, gridsize=500, cmap='viridis', mincnt=1)
        plt.colorbar(hb, label='Counts')
        plt.axvline(x=cy3_75th, color='black', linestyle='--', linewidth=0.8,
                    label=f'TxR 75th ({cy3_75th:.2f})')
        if len(highlight_cy3) >= 2:
            x_fit = np.linspace(highlight_cy3.min(), highlight_cy3.max(), 100)
            intercept_val = np.polyfit(highlight_cy3, highlight_txr, 1)[1]
            y_fit = slope_cy3_to_txr * x_fit + intercept_val
            plt.plot(x_fit, y_fit, color='green', linewidth=1,
                     label=f'Regression (slope={slope_cy3_to_txr:.4f})')
        plt.xlabel('Cy3 Intensity')
        plt.ylabel('TxR Intensity')
        plt.title(f'{base_name}: Density Plot (Cy3 → TxRed Correction)')
        plt.legend(loc='lower right')
        save_or_show(pdf, plt.gcf())
        
        # Density Plot 2: TxR → Cy3
        fig = plt.figure(figsize=(8,6))
        hb = plt.hexbin(valid_txr, valid_cy3, gridsize=500, cmap='viridis', mincnt=1)
        plt.colorbar(hb, label='Counts')
        plt.axvline(x=txr_75th, color='black', linestyle='--', linewidth=0.8,
                    label=f'TxRed 75th ({txr_75th:.2f})')
        if len(highlight_txr_2) >= 2:
            x_fit = np.linspace(highlight_txr_2.min(), highlight_txr_2.max(), 100)
            intercept_val = np.polyfit(highlight_txr_2, highlight_cy3_2, 1)[1]
            y_fit = slope_txr_to_cy3 * x_fit + intercept_val
            plt.plot(x_fit, y_fit, color='green', linewidth=1,
                     label=f'Regression (slope={slope_txr_to_cy3:.4f})')
        plt.xlabel('TxR Intensity')
        plt.ylabel('Cy3 Intensity')
        plt.title(f'{base_name}: Density Plot (TxR → Cy3 Correction)')
        plt.legend(loc='lower right')
        save_or_show(pdf, plt.gcf())
        
        # Density Plot 4: FAM → Cy3
        fig = plt.figure(figsize=(8,6))
        hb = plt.hexbin(valid_FAM_for_Cy3, valid_cy3_for_FAM, gridsize=500, cmap='viridis', mincnt=1)
        plt.colorbar(hb, label='Counts')
        plt.axvline(x=fam_75th_Cy3, color='black', linestyle='--', linewidth=0.8,
                    label=f'FAM 75th ({fam_75th_Cy3:.2f})')
        if len(highlight_FAM_Cy3) >= 2:
            x_fit = np.linspace(highlight_FAM_Cy3.min(), highlight_FAM_Cy3.max(), 100)
            intercept_val = np.polyfit(highlight_FAM_Cy3, highlight_cy3_for_FAM, 1)[1]
            y_fit = slope_FAM_to_Cy3 * x_fit + intercept_val
            plt.plot(x_fit, y_fit, color='green', linewidth=1,
                     label=f'Regression (slope={slope_FAM_to_Cy3:.4f})')
        plt.xlabel('FAM Intensity')
        plt.ylabel('Cy3 Intensity')
        plt.title(f'{base_name}: Density Plot (FAM → Cy3 Correction)')
        plt.legend(loc='lower right')
        save_or_show(pdf, plt.gcf())
        
        # Density Plot 5: TxR → Cy5
        fig = plt.figure(figsize=(8,6))
        hb = plt.hexbin(valid_txr_for_Cy5, valid_cy5_for_corr, gridsize=500, cmap='viridis', mincnt=1)
        plt.colorbar(hb, label='Counts')
        plt.axvline(x=txr_75th_Cy5, color='black', linestyle='--', linewidth=0.8,
                    label=f'TxR 75th ({txr_75th_Cy5:.2f})')
        if len(highlight_txr_Cy5) >= 2:
            x_fit = np.linspace(highlight_txr_Cy5.min(), highlight_txr_Cy5.max(), 100)
            intercept_val = np.polyfit(highlight_txr_Cy5, highlight_cy5_for_corr, 1)[1]
            y_fit = slope_txr_to_Cy5 * x_fit + intercept_val
            plt.plot(x_fit, y_fit, color='green', linewidth=1,
                     label=f'Regression (slope={slope_txr_to_Cy5:.4f})')
        plt.xlabel('TxR Intensity')
        plt.ylabel('Cy5 Intensity')
        plt.title(f'{base_name}: Density Plot (TxR → Cy5 Correction)')
        plt.legend(loc='lower right')
        save_or_show(pdf, plt.gcf())
        
        # Density Plot 6: Cy3 → FAM
        fig = plt.figure(figsize=(8,6))
        hb = plt.hexbin(valid_cy3_for_FAM, valid_FAM_for_corr, gridsize=500, cmap='viridis', mincnt=1)
        plt.colorbar(hb, label='Counts')
        plt.axvline(x=cy3_75th_FAM, color='black', linestyle='--', linewidth=0.8,
                    label=f'Cy3 75th ({cy3_75th_FAM:.2f})')
        if len(highlight_cy3_FAM) >= 2:
            x_fit = np.linspace(highlight_cy3_FAM.min(), highlight_cy3_FAM.max(), 100)
            intercept_val = np.polyfit(highlight_cy3_FAM, highlight_FAM_for_corr, 1)[1]
            y_fit = slope_Cy3_to_FAM * x_fit + intercept_val
            plt.plot(x_fit, y_fit, color='green', linewidth=1,
                     label=f'Regression (slope={slope_Cy3_to_FAM:.4f})')
        plt.xlabel('Cy3 Intensity')
        plt.ylabel('FAM Intensity')
        plt.title(f'{base_name}: Density Plot (Cy3 → FAM Correction)')
        plt.legend(loc='lower right')
        save_or_show(pdf, plt.gcf())
        
        # 7) Histograms for all channels
        # For Cy5:
        p70_cy5_before = np.percentile(masked_cy5, 70)
        p75_cy5_before = np.percentile(masked_cy5, 75)
        p80_cy5_before = np.percentile(masked_cy5, 80)
        corrected_cy5 = (ch_cy5 - slope_txr_to_Cy5 * ch_txr)[binary_mask]
        p70_cy5_after = np.percentile(corrected_cy5, 70)
        p75_cy5_after = np.percentile(corrected_cy5, 75)
        p80_cy5_after = np.percentile(corrected_cy5, 80)
        
        fig = plt.figure(figsize=(12,5))
        plt.subplot(1,2,1)
        plt.hist(masked_cy5, bins=500, color='green', alpha=0.3)
        plt.axvline(p70_cy5_before, color='red', linestyle=':', linewidth=1, label='70th Percentile')
        plt.axvline(p75_cy5_before, color='black', linestyle=':', linewidth=1, label='75th Percentile')
        plt.axvline(p80_cy5_before, color='blue', linestyle=':', linewidth=1, label='80th Percentile')
        plt.title(f'{base_name}: Masked Cy5 (Before Correction)')
        plt.gca().yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,.0f}'))
        plt.legend()
        plt.subplot(1,2,2)
        plt.hist(corrected_cy5, bins=500, color='green', alpha=0.3)
        plt.axvline(p70_cy5_after, color='red', linestyle=':', linewidth=1, label='70th Percentile')
        plt.axvline(p75_cy5_after, color='black', linestyle=':', linewidth=1, label='75th Percentile')
        plt.axvline(p80_cy5_after, color='blue', linestyle=':', linewidth=1, label='80th Percentile')
        plt.title(f'{base_name}: Corrected Cy5 (After Correction)')
        plt.gca().yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,.0f}'))
        plt.legend()
        save_or_show(pdf, plt.gcf())
        
        # For Cy3:
        p70_cy3_before = np.percentile(masked_cy3, 70)
        p75_cy3_before = np.percentile(masked_cy3, 75)
        p80_cy3_before = np.percentile(masked_cy3, 80)
        corrected_cy3 = (ch_cy3 - slope_txr_to_cy3 * ch_txr - slope_FAM_to_Cy3 * ch_FAM)[binary_mask]
        p70_cy3_after = np.percentile(corrected_cy3, 70)
        p75_cy3_after = np.percentile(corrected_cy3, 75)
        p80_cy3_after = np.percentile(corrected_cy3, 80)
        
        fig = plt.figure(figsize=(12,5))
        plt.subplot(1,2,1)
        plt.hist(masked_cy3, bins=500, color='orange', alpha=0.3)
        plt.axvline(p70_cy3_before, color='red', linestyle=':', linewidth=1, label='70th Percentile')
        plt.axvline(p75_cy3_before, color='black', linestyle=':', linewidth=1, label='75th Percentile')
        plt.axvline(p80_cy3_before, color='blue', linestyle=':', linewidth=1, label='80th Percentile')
        plt.title(f'{base_name}: Masked Cy3 (Before Correction)')
        plt.gca().yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,.0f}'))
        plt.legend()
        plt.subplot(1,2,2)
        plt.hist(corrected_cy3, bins=500, color='orange', alpha=0.3)
        plt.axvline(p70_cy3_after, color='red', linestyle=':', linewidth=1, label='70th Percentile')
        plt.axvline(p75_cy3_after, color='black', linestyle=':', linewidth=1, label='75th Percentile')
        plt.axvline(p80_cy3_after, color='blue', linestyle=':', linewidth=1, label='80th Percentile')
        plt.ylim(0, 350000)
        plt.title(f'{base_name}: Corrected Cy3 (After Correction)')
        plt.gca().yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,.0f}'))
        plt.legend()
        save_or_show(pdf, plt.gcf())
        
        # For TxR:
        p70_txr_before = np.percentile(masked_txr, 70)
        p75_txr_before = np.percentile(masked_txr, 75)
        p80_txr_before = np.percentile(masked_txr, 80)
        corrected_txr = (ch_txr - slope_cy3_to_txr * ch_cy3)[binary_mask]
        p70_txr_after = np.percentile(corrected_txr, 70)
        p75_txr_after = np.percentile(corrected_txr, 75)
        p80_txr_after = np.percentile(corrected_txr, 80)
        
        fig = plt.figure(figsize=(12,5))
        plt.subplot(1,2,1)
        plt.hist(masked_txr, bins=500, color='purple', alpha=0.3)
        plt.axvline(p70_txr_before, color='red', linestyle=':', linewidth=1, label='70th Percentile')
        plt.axvline(p75_txr_before, color='black', linestyle=':', linewidth=1, label='75th Percentile')
        plt.axvline(p80_txr_before, color='blue', linestyle=':', linewidth=1, label='80th Percentile')
        plt.title(f'{base_name}: Masked TxR (Before Correction)')
        plt.gca().yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,.0f}'))
        plt.legend()
        plt.subplot(1,2,2)
        plt.hist(corrected_txr, bins=500, color='purple', alpha=0.3)
        plt.axvline(p70_txr_after, color='red', linestyle=':', linewidth=1, label='70th Percentile')
        plt.axvline(p75_txr_after, color='black', linestyle=':', linewidth=1, label='75th Percentile')
        plt.axvline(p80_txr_after, color='blue', linestyle=':', linewidth=1, label='80th Percentile')
        plt.title(f'{base_name}: Corrected TxR (After Correction)')
        plt.gca().yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,.0f}'))
        plt.legend()
        save_or_show(pdf, plt.gcf())
        
        # For FAM:
        p70_FAM_before = np.percentile(masked_FAM, 70)
        p75_FAM_before = np.percentile(masked_FAM, 75)
        p80_FAM_before = np.percentile(masked_FAM, 80)
        corrected_FAM = (ch_FAM - slope_Cy3_to_FAM * ch_cy3)[binary_mask]
        p70_FAM_after = np.percentile(corrected_FAM, 70)
        p75_FAM_after = np.percentile(corrected_FAM, 75)
        p80_FAM_after = np.percentile(corrected_FAM, 80)
        
        fig = plt.figure(figsize=(12,5))
        plt.subplot(1,2,1)
        plt.hist(masked_FAM, bins=500, color='brown', alpha=0.3)
        plt.axvline(p70_FAM_before, color='red', linestyle=':', linewidth=1, label='70th Percentile')
        plt.axvline(p75_FAM_before, color='black', linestyle=':', linewidth=1, label='75th Percentile')
        plt.axvline(p80_FAM_before, color='blue', linestyle=':', linewidth=1, label='80th Percentile')
        plt.title(f'{base_name}: Masked FAM (Before Correction)')
        plt.gca().yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,.0f}'))
        plt.legend()
        plt.subplot(1,2,2)
        plt.hist(corrected_FAM, bins=500, color='brown', alpha=0.3)
        plt.axvline(p70_FAM_after, color='red', linestyle=':', linewidth=1, label='70th Percentile')
        plt.axvline(p75_FAM_after, color='black', linestyle=':', linewidth=1, label='75th Percentile')
        plt.axvline(p80_FAM_after, color='blue', linestyle=':', linewidth=1, label='80th Percentile')
        plt.title(f'{base_name}: Corrected FAM (After Correction)')
        plt.gca().yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,.0f}'))
        plt.legend()
        save_or_show(pdf, plt.gcf())
    
    # 8) Apply corrections to form final image and scale to uint16
    corrected_txr = np.clip(ch_txr - slope_cy3_to_txr * ch_cy3, 0, None)
    corrected_cy3 = np.clip(ch_cy3 - slope_txr_to_cy3 * ch_txr - slope_FAM_to_Cy3 * ch_FAM, 0, None)
    corrected_cy5 = np.clip(ch_cy5 - slope_txr_to_Cy5 * ch_txr, 0, None)
    corrected_FAM = np.clip(ch_FAM - slope_Cy3_to_FAM * ch_cy3, 0, None)
    
    corrected_image = np.dstack([corrected_cy5, corrected_txr, corrected_cy3, corrected_FAM])
    min_val, max_val = corrected_image.min(), corrected_image.max()
    if max_val == min_val:
        corrected_image_uint16 = np.zeros_like(corrected_image, dtype=np.uint16)
    else:
        corrected_image_uint16 = ((65535 * (corrected_image - min_val) / (max_val - min_val)).astype(np.uint16))
    
    return corrected_image_uint16, slope_cy3_to_txr, slope_FAM_to_TxR, slope_txr_to_cy3, slope_FAM_to_Cy3, slope_txr_to_Cy5, slope_Cy3_to_FAM

# ===============================================
# Main loop: process all bases and save all figures to a PDF
# ===============================================
corrected_base_images = {}
slopes_by_base = {}  # Dictionary to store slopes for each base

with PdfPages(pdf_path) as pdf:
    for base_name, image in base_tif_images.items():
        if base_name == 'align':  
            continue  # Skip non-image files
        else: 
            base_image = image
            eroded_binary_mask = eroded_binary_masks[base_name]
            print(f'Processing {base_name}...')
            corrected_img, slope_cy3_to_txr, slope_FAM_to_TxR, slope_txr_to_cy3, \
            slope_FAM_to_Cy3, slope_txr_to_Cy5, slope_Cy3_to_FAM = correct_bleedthroughs(
                original_image=base_image,
                binary_mask=eroded_binary_mask,
                fallback_slope_txr=0.5,
                fallback_slope_cy3=0.5,
                show_plots=True,
                base_name=base_name,
                pdf=pdf
            )
            corrected_base_images[base_name] = corrected_img
            slopes_by_base[base_name] = (slope_cy3_to_txr, slope_txr_to_cy3, slope_FAM_to_TxR,
                                          slope_FAM_to_Cy3, slope_txr_to_Cy5, slope_Cy3_to_FAM)
            
            # Zoomed-in views for each channel
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            axes[0].imshow(base_image[3000:3500, 3500:4000, 1], cmap='gray')
            axes[0].set_title(f'{base_name} - Original TxR (500x500)')
            axes[1].imshow(corrected_img[3000:3500, 3500:4000, 1], cmap='gray')
            axes[1].set_title(f'{base_name} - Corrected TxR (500x500)')
            save_or_show(pdf, plt.gcf())
            
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            axes[0].imshow(base_image[3000:3500, 3500:4000, 2], cmap='gray')
            axes[0].set_title(f'{base_name} - Original Cy3 (500x500)')
            axes[1].imshow(corrected_img[3000:3500, 3500:4000, 2], cmap='gray')
            axes[1].set_title(f'{base_name} - Corrected Cy3 (500x500)')
            save_or_show(pdf, plt.gcf())
            
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            axes[0].imshow(base_image[3000:3500, 3500:4000, 0], cmap='gray')
            axes[0].set_title(f'{base_name} - Original Cy5 (500x500)')
            axes[1].imshow(corrected_img[3000:3500, 3500:4000, 0], cmap='gray')
            axes[1].set_title(f'{base_name} - Corrected Cy5 (500x500)')
            save_or_show(pdf, plt.gcf())
            
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            axes[0].imshow(base_image[3000:3500, 3500:4000, 3], cmap='gray')
            axes[0].set_title(f'{base_name} - Original FAM (500x500)')
            axes[1].imshow(corrected_img[3000:3500, 3500:4000, 3], cmap='gray')
            axes[1].set_title(f'{base_name} - Corrected FAM (500x500)')
            save_or_show(pdf, plt.gcf())
            
            print(f"{base_name}: Slopes -> Cy3→TxR: {slope_cy3_to_txr:.4f}, TxR→Cy3: {slope_txr_to_cy3:.4f}")
            print(f"{base_name}: Slopes -> FAM→TxR: {slope_FAM_to_TxR:.4f}, FAM→Cy3: {slope_FAM_to_Cy3:.4f}")
            print(f"{base_name}: Slopes -> TxR→Cy5: {slope_txr_to_Cy5:.4f}, Cy3→FAM: {slope_Cy3_to_FAM:.4f}")
        # else:
        #     print(f"No mask for {base_name}. Skipping.")
    
    # Add a summary page with all slope outputs
    fig, ax = plt.subplots(figsize=(8.5, 11))
    ax.axis('off')
    summary_lines = ["Bleedthrough Correction Slopes:\n"]
    for base, slopes in slopes_by_base.items():
        summary_lines.append(f"{base}:")
        summary_lines.append(f"  Cy3 → TxR: {slopes[0]:.4f}")
        summary_lines.append(f"  TxR → Cy3: {slopes[1]:.4f}")
        summary_lines.append(f"  FAM → TxR: {slopes[2]:.4f}")
        summary_lines.append(f"  FAM → Cy3: {slopes[3]:.4f}")
        summary_lines.append(f"  TxR → Cy5: {slopes[4]:.4f}")
        summary_lines.append(f"  Cy3 → FAM: {slopes[5]:.4f}\n")
    summary_text = "\n".join(summary_lines)
    ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, fontsize=10,
            va='top', family='monospace')
    save_or_show(pdf, plt.gcf())

#%%
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#STEP 6A: get centroids for every bead
import pandas as pd
import numpy as np
from skimage.measure import regionprops_table

# Dictionary to store DataFrames for each base (image)
base_call_dfs = {}

# Iterate through each base image and corresponding binary mask
for base_var, binary_mask in binary_masks.items():
    # if base_var in corrected_base_images and base_var in segmentations:
        # Get the corrected image and segmentation
    print(f'{base_var} running')
    segmentation = segmentations[base_var]

    # Use regionprops_table to get centroids for each segmentation
    props = regionprops_table(segmentation, properties=['label', 'centroid'])

    # Convert to df
    df_centroids = pd.DataFrame.from_dict(props)

    # Rename columns for clarity
    df_centroids.columns = ['label', 'centroid_x', 'centroid_y']

    # Set index by 'label' for easy reference 
    df_centroids.set_index('label', inplace=True)

    # Store the final df for this base in the dictionary
    base_call_dfs[base_var] = df_centroids


#Save the dictionary to a file
# import pickle
# with open('~/dir/base_call_dfs.pkl', 'wb') as file:
#     pickle.dump(base_call_dfs, file)


# to load dictionary from file:
import pickle
with open('~/dir/base_call_dfs.pkl', 'rb') as file:
    base_call_dfs = pickle.load(file)

#%%
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#STEP 6B: Save centroids separately, for ICP
# Dictionary to store centroid data for each image
centroid_dict = {}

# Extract centroids from each merged DataFrame
for base_var, merged_df in base_call_dfs.items():
    # Extract only the 'centroid_x' and 'centroid_y' columns
    centroids = merged_df[['centroid_x', 'centroid_y']].values  # Convert to numpy array
    centroid_dict[base_var] = centroids  # Store the centroids in the dictionary
    

# To save dictionary to file:
import pickle
with open('~/dir/centroid_dict.pkl', 'wb') as file:
    pickle.dump(centroid_dict, file)

# to load dictionary from file:
import pickle
with open('~/dir/centroid_dict.pkl', 'rb') as file:
    centroid_dict = pickle.load(file)


#%%
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#STEP 7A: transforming images to be on the same coordinate system using ICP: ALL FUNCTIONS
import open3d as o3d
import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt

# Step 1: Filter for edge points
def filter_edge_points(centroids, image_shape, edge_threshold=650, edges=None):
    """
    Filters centroids to only include those near specified edges of the image.

    Parameters:
    - centroids: The array of centroid coordinates.
    - image_shape: The shape of the image (height, width).
    - edge_threshold: Distance from the edge to consider as 'near the edge'.
    - edges: List of edges to filter for ('top', 'bottom', 'left', 'right').
             If None, defaults to 'all' edges.

    Returns:
    - edge_centroids: Filtered array of edge centroids.
    """
    if edges is None:
        edges = ['top', 'bottom', 'left', 'right']  # Default to all edges
    
    height, width = image_shape
    x, y = centroids[:, 0], centroids[:, 1]

    # Initialize mask for selected edges
    edge_mask = np.zeros_like(x, dtype=bool)

    # Apply masks for specified edges
    if 'top' in edges:
        edge_mask |= (y < edge_threshold)
    if 'bottom' in edges:
        edge_mask |= (y > height - edge_threshold)
    if 'left' in edges:
        edge_mask |= (x < edge_threshold)
    if 'right' in edges:
        edge_mask |= (x > width - edge_threshold)

    # Apply the mask to filter centroids near the specified edges
    edge_centroids = centroids[edge_mask]

    return edge_centroids



# Step 2: Apply ICP on the matched centroids with a given threshold

def apply_icp(source_points, target_points, threshold=300, trans_init=None):
    # Convert centroids into Open3D point cloud format
    sourcepoints_3d = np.hstack((source_points, np.zeros((source_points.shape[0], 1))))
    targetpoints_3d = np.hstack((target_points, np.zeros((target_points.shape[0], 1)))) 

    source_pcd = o3d.geometry.PointCloud()
    source_pcd.points = o3d.utility.Vector3dVector(sourcepoints_3d)

    target_pcd = o3d.geometry.PointCloud()
    target_pcd.points = o3d.utility.Vector3dVector(targetpoints_3d)

    # Apply trans_init if provided, directly transforming source_pcd
    if trans_init is not None:
        source_pcd.transform(trans_init)
        print("Applied initial manual transformation:")
        print(trans_init)
    else:
        trans_init = np.eye(4)  # Default to identity (no initial transformation)

    # Set up convergence criteria with higher max iterations
    convergence_criteria = o3d.pipelines.registration.ICPConvergenceCriteria(
        max_iteration=20000,
        relative_fitness=1e-6,
        relative_rmse=1e-6
    )

    # Apply ICP with identity transformation as the starting point for iterative refinement
    reg_p2p = o3d.pipelines.registration.registration_icp(
        source_pcd, target_pcd, threshold, np.eye(4),
        o3d.pipelines.registration.TransformationEstimationPointToPoint(),
        convergence_criteria
    )

    transformation_matrix = reg_p2p.transformation

    # QC: Extract fitness and inlier RMSE
    fitness = reg_p2p.fitness
    inlier_rmse = reg_p2p.inlier_rmse

    print(f"ICP Transformation Matrix: \n{transformation_matrix}")
    print(f"Fitness: {fitness}")
    print(f"Inlier RMSE: {inlier_rmse}")

    # Combine the manual trans_init with the resulting transformation matrix from ICP
    final_transformation = transformation_matrix @ trans_init

    return final_transformation, fitness, inlier_rmse


# Step 3: function to visualize source and target points
import matplotlib.pyplot as plt

def plot_centroids(source_centroids, target_centroids):
    # First figure and axes for overall and initial zoomed-in plots
    fig, axes = plt.subplots(1, 2, figsize=(15, 7))

    # Plot overall view
    axes[0].scatter(source_centroids[:, 0], source_centroids[:, 1], color='orange', label='Source Centroids', alpha=0.3, s=0.2)
    axes[0].scatter(target_centroids[:, 0], target_centroids[:, 1], color='blue', label='Target Centroids', alpha=0.3, s=0.2)
    axes[0].set_title(f'{base_var}: Overall Centroids')
    axes[0].legend(loc='lower left')

    # First zoomed-in view
    zoom_x_min, zoom_x_max = 7300, 7800
    zoom_y_min, zoom_y_max = 900, 1400
    axes[1].scatter(source_centroids[:, 0], source_centroids[:, 1], color='orange', alpha=0.5)
    axes[1].scatter(target_centroids[:, 0], target_centroids[:, 1], color='blue', alpha=0.3)
    axes[1].set_xlim(zoom_x_min, zoom_x_max)
    axes[1].set_ylim(zoom_y_min, zoom_y_max)
    axes[1].set_title(f'{base_var}: Zoomed-in Centroids ({zoom_x_min}x{zoom_y_min} to {zoom_x_max}x{zoom_y_max})')

    plt.tight_layout()
    plt.show()

    # Second figure for two more zoomed-in plot
    fig, ax = plt.subplots(1, 2, figsize=(15, 7))

    zoom_x_min, zoom_x_max = 200, 1100
    zoom_y_min, zoom_y_max = 7400, 8300
    ax[0].scatter(source_centroids[:, 0], source_centroids[:, 1], color='orange', alpha=0.5)
    ax[0].scatter(target_centroids[:, 0], target_centroids[:, 1], color='blue', alpha=0.3)
    ax[0].set_xlim(zoom_x_min, zoom_x_max)
    ax[0].set_ylim(zoom_y_min, zoom_y_max)
    ax[0].set_title(f'{base_var}: Zoomed-in Centroids ({zoom_x_min}x{zoom_y_min} to {zoom_x_max}x{zoom_y_max})')

    zoom_x_min, zoom_x_max = 7000, 7800
    zoom_y_min, zoom_y_max = 7400, 8200
    ax[1].scatter(source_centroids[:, 0], source_centroids[:, 1], color='orange', alpha=0.5)
    ax[1].scatter(target_centroids[:, 0], target_centroids[:, 1], color='blue', alpha=0.3)
    ax[1].set_xlim(zoom_x_min, zoom_x_max)
    ax[1].set_ylim(zoom_y_min, zoom_y_max)
    ax[1].set_title(f'{base_var}: Zoomed-in Centroids ({zoom_x_min}x{zoom_y_min} to {zoom_x_max}x{zoom_y_max})')

    plt.tight_layout()
    plt.show()

    # Third figure for more zoomed-in plots
    fig, ax = plt.subplots(1, 2, figsize=(15, 7))

    zoom_x_min, zoom_x_max = 3000, 3500
    zoom_y_min, zoom_y_max = 3000, 3500
    ax[0].scatter(source_centroids[:, 0], source_centroids[:, 1], color='orange', alpha=0.5)
    ax[0].scatter(target_centroids[:, 0], target_centroids[:, 1], color='blue', alpha=0.3)
    ax[0].set_xlim(zoom_x_min, zoom_x_max)
    ax[0].set_ylim(zoom_y_min, zoom_y_max)
    ax[0].set_title(f'{base_var}: Zoomed-in Centroids ({zoom_x_min}x{zoom_y_min} to {zoom_x_max}x{zoom_y_max})')

    zoom_x_min, zoom_x_max = 6000, 6500
    zoom_y_min, zoom_y_max = 6000, 6500
    ax[1].scatter(source_centroids[:, 0], source_centroids[:, 1], color='orange', alpha=0.5)
    ax[1].scatter(target_centroids[:, 0], target_centroids[:, 1], color='blue', alpha=0.3)
    ax[1].set_xlim(zoom_x_min, zoom_x_max)
    ax[1].set_ylim(zoom_y_min, zoom_y_max)
    ax[1].set_title(f'{base_var}: Zoomed-in Centroids ({zoom_x_min}x{zoom_y_min} to {zoom_x_max}x{zoom_y_max})')

    plt.tight_layout()
    plt.show()
    
    # Fourth figure for more zoomed-in plots
    fig, ax = plt.subplots(1, 2, figsize=(15, 7))

    zoom_x_min, zoom_x_max = 2000, 2500
    zoom_y_min, zoom_y_max = 2000, 2500
    ax[0].scatter(source_centroids[:, 0], source_centroids[:, 1], color='orange', alpha=0.5)
    ax[0].scatter(target_centroids[:, 0], target_centroids[:, 1], color='blue', alpha=0.3)
    ax[0].set_xlim(zoom_x_min, zoom_x_max)
    ax[0].set_ylim(zoom_y_min, zoom_y_max)
    ax[0].set_title(f'{base_var}: Zoomed-in Centroids ({zoom_x_min}x{zoom_y_min} to {zoom_x_max}x{zoom_y_max})')

    zoom_x_min, zoom_x_max = 1000, 1500
    zoom_y_min, zoom_y_max = 1000, 1500
    ax[1].scatter(source_centroids[:, 0], source_centroids[:, 1], color='orange', alpha=0.5)
    ax[1].scatter(target_centroids[:, 0], target_centroids[:, 1], color='blue', alpha=0.3)
    ax[1].set_xlim(zoom_x_min, zoom_x_max)
    ax[1].set_ylim(zoom_y_min, zoom_y_max)
    ax[1].set_title(f'{base_var}: Zoomed-in Centroids ({zoom_x_min}x{zoom_y_min} to {zoom_x_max}x{zoom_y_max})')

    plt.tight_layout()
    plt.show()

    # Fifth figure for more zoomed-in plots
    fig, ax = plt.subplots(1, 2, figsize=(15, 7))

    zoom_x_min, zoom_x_max = 200, 800
    zoom_y_min, zoom_y_max = 3900, 4500
    ax[0].scatter(source_centroids[:, 0], source_centroids[:, 1], color='orange', alpha=0.5)
    ax[0].scatter(target_centroids[:, 0], target_centroids[:, 1], color='blue', alpha=0.3)
    ax[0].set_xlim(zoom_x_min, zoom_x_max)
    ax[0].set_ylim(zoom_y_min, zoom_y_max)
    ax[0].set_title(f'{base_var}: Zoomed-in Centroids ({zoom_x_min}x{zoom_y_min} to {zoom_x_max}x{zoom_y_max})')

    zoom_x_min, zoom_x_max = 900, 1500
    zoom_y_min, zoom_y_max = 3900, 4500
    ax[1].scatter(source_centroids[:, 0], source_centroids[:, 1], color='orange', alpha=0.5)
    ax[1].scatter(target_centroids[:, 0], target_centroids[:, 1], color='blue', alpha=0.3)
    ax[1].set_xlim(zoom_x_min, zoom_x_max)
    ax[1].set_ylim(zoom_y_min, zoom_y_max)
    ax[1].set_title(f'{base_var}: Zoomed-in Centroids ({zoom_x_min}x{zoom_y_min} to {zoom_x_max}x{zoom_y_max})')

    plt.tight_layout()
    plt.show()

# Step 4: Threshold sweep function
def find_best_threshold(source_centroids, target_centroids, threshold_range=np.arange(0, 600, 2), trans_init=None): #sweep from 0 to 300 in increments of 1
    best_threshold = None
    best_rmse = float('inf')
    best_fitness = 0
    best_transformation = None

    # Sweep across different thresholds
    for threshold in threshold_range:
        if threshold == 0:  # Skip threshold 0 since it is invalid
            continue

        print(f"Testing threshold: {threshold}")
        transformation_matrix, fitness, inlier_rmse = apply_icp(
            source_centroids, target_centroids, threshold, trans_init=trans_init
        ) # Pass the trans_init to apply_icp

        # Check if fitness is above 0.7 and RMSE is lower than the current best
        if fitness >= 0.6 and inlier_rmse < best_rmse:
            best_rmse = inlier_rmse
            best_fitness = fitness
            best_threshold = threshold
            best_transformation = transformation_matrix

        #print(f"Threshold: {threshold}, Fitness: {fitness}, RMSE: {inlier_rmse}")

    print(f"Best threshold: {best_threshold}, Best Fitness: {best_fitness}, Best RMSE: {best_rmse}")
    return best_threshold, best_transformation, best_fitness, best_rmse


# Step 5: Test the function with example centroids

def apply_full_transformation(source_centroids, target_centroids, edge_threshold=900, trans_init=None):
    image_shape = (8800, 8800)  # Example image size

    # Filter for edge points only
    source_edge_centroids = filter_edge_points(source_centroids, image_shape, edge_threshold=edge_threshold, edges=['right', 'bottom'])
    target_edge_centroids = filter_edge_points(target_centroids, image_shape, edge_threshold=edge_threshold, edges=['right', 'bottom'])
    
    # Call the threshold sweep function, passing in trans_init
    best_threshold, best_transformation, best_fitness, best_rmse = find_best_threshold(
        source_edge_centroids, target_edge_centroids, trans_init=trans_init
    )
    
    # Check if `best_transformation` is None or not set properly
    if best_transformation is None:
        print("Warning: ICP did not return a valid transformation matrix. Using `trans_init` as fallback.")
        best_transformation = trans_init if trans_init is not None else np.eye(4)

    # Apply the transformation to all source centroids
    source_centroids_3d = np.hstack((source_centroids, np.zeros((source_centroids.shape[0], 1))))  
    transformed_source_points_3d = (best_transformation @ np.hstack((source_centroids_3d, np.ones((source_centroids_3d.shape[0], 1)))).T).T[:, :2]
    
    # Plot centroids before and after ICP
    plot_centroids(source_centroids, target_centroids)  # Before ICP (All points)
    plot_centroids(transformed_source_points_3d, target_centroids)  # After ICP
    
    
    # Create a side-by-side comparison of original and transformed points
    import pandas as pd
    comparison = np.hstack((source_centroids, transformed_source_points_3d))
    df_comparison = pd.DataFrame(comparison, columns=["original_x", "original_y", "transformed_x", "transformed_y"])

    
    return df_comparison, best_transformation


# Step 5: Test the function with example centroids

# source_centroids = centroid_dict['base02']  # Source image centroids
# target_centroids = centroid_dict['base01']  # Target image centroids

# apply_full_transformation(source_centroids, target_centroids, edge_threshold=300)


#%%
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
# STEP 7B: apply transformation from every image to the coordinate system of 'base02'
# Function to visualize all aligned centroids
def plot_aligned_centroids(centroid_dict, base_var, transformed_dict):
    plt.figure(figsize=(10, 10))
    
    # Use the updated colormap with the correct usage
    cmap = plt.colormaps.get_cmap('tab10')
    colors = [cmap(i / len(centroid_dict)) for i in range(len(centroid_dict))]  # Get distinct colors

    # Plot each base's centroids
    for i, (base_name, transformed_centroids_df) in enumerate(transformed_dict.items()):
        # Extract transformed x and y from the DataFrame
        plt.scatter(transformed_centroids_df['transformed_x'], transformed_centroids_df['transformed_y'], 
                    color=colors[i], label=f'{base_name} (aligned)', alpha=0.2, s=0.4)
    
    # Also plot the target base centroids (base_var)
    plt.scatter(centroid_dict[base_var][:, 0], centroid_dict[base_var][:, 1], 
                color='gray', label=f'{base_var} (reference)', alpha=0.2, s=0.4)

    plt.title(f'Aligned Centroids to {base_var}')
    plt.legend(loc='lower left')  # Move legend to bottom-left
    plt.grid(True)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()


# Main Process: Align all centroids to 'base02'
target_base = 'base02'
target_centroids = centroid_dict[target_base]

# Create an empty dictionary to store DataFrames for each base
df_transformed_centroids = {}
ICP_transformation_matrices = {}

for base_var, source_centroids in centroid_dict.items():
    if base_var != target_base:

        # USE THIS ONLY IF NEEDED -- in cases where base transformed very poorly
        # Provide custom transformation matrix for certain bases (example below)
        # if base_var == 'base02':
        #     trans_init = np.array([
        #         [np.cos(np.deg2rad(-3)), -np.sin(np.deg2rad(-3)), 0, 0],
        #         [np.sin(np.deg2rad(-3)), np.cos(np.deg2rad(-3)), 0, 0],  # Translate 50 units in up direction
        #         [0, 0, 1, 0],
        #         [0, 0, 0, 1]
        #     ])
        # elif base_var == 'base03':
        #     trans_init = np.array([
        #         [np.cos(np.deg2rad(-3)), -np.sin(np.deg2rad(-3)), 0, 0],
        #         [np.sin(np.deg2rad(-3)), np.cos(np.deg2rad(-3)), 0, 0],  # Translate 100 units UPWARD
        #         [0, 0, 1, 0],
        #         [0, 0, 0, 1]
        #     ])
        # else:
        #     trans_init = np.eye(4)  # Use default behavior for others (identity transformation)
        
        # Apply the transformation (returns a tuple with DataFrame and transformation matrix)
        df_comparison, best_transformation = apply_full_transformation(source_centroids, target_centroids, edge_threshold=500, trans_init=trans_init)

        # Store the best transformation matrix in the dictionary
        ICP_transformation_matrices[base_var] = best_transformation

        # Save the transformed DataFrame in the dictionary
        df_transformed_centroids[base_var] = df_comparison

    else:
        # For the target centroids, store them directly (no transformation)
        df_target = pd.DataFrame({
            'original_x': target_centroids[:, 0],  # Original is the same as transformed
            'original_y': target_centroids[:, 1],
            'transformed_x': target_centroids[:, 0],
            'transformed_y': target_centroids[:, 1]
        })
        df_transformed_centroids[base_var] = df_target


# Visualize the aligned centroids
plot_aligned_centroids(centroid_dict, target_base, df_transformed_centroids)


# To save dictionary to file:
import pickle
with open('~/dir/ICP_transformation_matrices.pkl', 'wb') as file:
    pickle.dump(ICP_transformation_matrices, file)

# to load dictionary from file:
import pickle
with open('~/dir/ICP_transformation_matrices.pkl', 'rb') as file:
    ICP_transformation_matrices = pickle.load(file)


#%%
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#STEP 7C: add transformed coordinates to each corresponding df in 'base_call_dfs'

base_calls_with_coords = {}

for base_var, base_call_df in base_call_dfs.items():
    if base_var == 'align':
        continue
    else:
        # Add the 'transformed_x' and 'transformed_y' columns
        base_calls_with_coords[base_var] = base_call_df.merge(df_transformed_centroids[base_var], right_index=True, left_index=True)
    
for base_var, base_call_df in base_call_dfs.items():
    if base_var == 'align':
        continue
    else:
        # Perform the merge based on the 'centroid_x' column in base_call_df and 'original_x' in df_transformed_centroids
        base_calls_with_coords[base_var] = base_call_df.merge(
                                    df_transformed_centroids[base_var], 
                                    left_on='centroid_x', 
                                    right_on='original_x', 
                                    how='inner'  # Choose 'inner', 'left', 'right', or 'outer' depending on your desired behavior
                                )
#%%
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#STEP 8: Find nearest centroids / assign indices
    #Steps:
    # Assign index to base01: Assign each centroid in base01 a unique index
    # Compare base02 to base01: Compare centroids in base02 to base01. 
        #   If a centroid in base02 is within 13 units of a base01 centroid 
        #   in both transformed_x and transformed_y, assign it the same index. 
        #   If it doesn't match any base01 centroid, assign it a new index.
    # Repeat for base03 through base14: 
        #   Compare each subsequent base to all previously indexed centroids 
        #   from base01, base02, etc., following the same logic.
    # Track base calls: Track which base image each base call came from 
        #   and include that information in the final table.

import pandas as pd
import numpy as np
from scipy.spatial import cKDTree

# Define threshold for matching centroids within 10 units
MATCHING_THRESHOLD = 12

# Step 1: Assign unique indices to base02 centroids
def initialize_indices(base_df):
    base_df['index'] = np.arange(len(base_df))
    return base_df

# Step 2: Compare and assign index from base01 to base02 (and further)
def assign_indices(df_to_match, reference_df, matching_threshold):
    # Create KDTree for the reference_df to find nearest neighbors
    ref_tree = cKDTree(reference_df[['transformed_x', 'transformed_y']].values)

    # Query for nearest neighbor in the reference base within the threshold
    distances, indices = ref_tree.query(df_to_match[['transformed_x', 'transformed_y']].values, distance_upper_bound=matching_threshold)

    # Filter out the invalid matches (indices that are out of bounds)
    valid_mask = (distances < matching_threshold) & (indices < len(reference_df))

    # Create a new index column in df_to_match
    df_to_match['index'] = np.nan  # Initialize with NaN
    df_to_match.loc[valid_mask, 'index'] = reference_df.iloc[indices[valid_mask]]['index'].values

    # Assign new index to unmatched centroids (where index is still NaN)
    unmatched_mask = df_to_match['index'].isna()
    max_index = reference_df['index'].max()

    # Directly assign new indices to unmatched rows
    df_to_match.loc[unmatched_mask, 'index'] = np.arange(max_index + 1, max_index + 1 + unmatched_mask.sum())

    return df_to_match

# Step 3: Repeat for all subsequent bases, comparing with previously matched bases
def match_all_bases_to_base02(base_calls_with_coords, threshold):
    """
    Matches all other bases to the reference base 'base02'.

    Parameters:
    - base_calls_with_coords: Dictionary containing dataframes for all bases (e.g., base01, base02, ..., base14, base06).
    - threshold: Distance threshold for nearest neighbor matching.

    Returns:
    - Updated base_calls_with_coords with aligned indices.
    """
    # Step 1: Initialize indices for base02 (reference base)
    base02_df = initialize_indices(base_calls_with_coords['base02'])

    # Store the updated dataframe for base02
    base_calls_with_coords['base02'] = base02_df

    # Step 2: Iterate through all bases except base06
    for i in range(1, 15):  # base01 to base14
        if i == 2:
            continue  # Skip base02 since it is the reference base

        base_var = f'base{i:02}'
        df_to_match = base_calls_with_coords[base_var]

        # Compare each base directly to base06
        df_to_match = assign_indices(df_to_match, base02_df, threshold)

        # Store the updated dataframe back into the dictionary
        base_calls_with_coords[base_var] = df_to_match

    return base_calls_with_coords

# Step 4: Track which base each index/base call came from
def create_final_output(base_calls_with_coords):
    output_list = []

    # Collect the relevant data for each base and append to the output list
    for base_var, df in base_calls_with_coords.items():
        for _, row in df.iterrows():
            output_list.append({
                'index': int(row['index']),
                'x': row['transformed_x'],
                'y': row['transformed_y'],
                'base': base_var#,
                #'base_call': row['base_call']
            })

    # Create a DataFrame with the final output
    final_output_df = pd.DataFrame(output_list)
    return final_output_df

# Example usage:
# Assuming base_calls_with_coords is the dictionary with 14 DataFrames

# Step 5: Match all bases and assign indices
base_calls_with_coords = match_all_bases_to_base02(base_calls_with_coords, MATCHING_THRESHOLD)

## Step 6: Create the final output
final_output = create_final_output(base_calls_with_coords)

# Show the output
print(final_output.head())


# to save the dictionary to a file
import pickle
with open('~/dir/base_calls_with_coords.pkl', 'wb') as file:
    pickle.dump(base_calls_with_coords, file)

# to load dictionary from file:
import pickle
with open('~/dir/base_calls_with_coords.pkl', 'rb') as file:
    base_calls_with_coords = pickle.load(file)


#%%
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#STEP 9: Generate a master mask for composite puck coords 
#base_calls_with_coords: has original and transformed coordinates + indices
import pandas as pd

# Step 1: Combine all DataFrames from the dictionary into one giant DataFrame
combined_df = pd.concat(base_calls_with_coords.values(), ignore_index=True)

# Step 2: Select only the columns 'transformed_x', 'transformed_y', and 'index'
selected_columns = combined_df[['transformed_x', 'transformed_y', 'index']]


# Step 3: Collapse by 'index' and calculate avg of 'transformed_x' and 'transformed_y' to become one centroid
collapsed_df = selected_columns.groupby('index', as_index=False).agg({
    'transformed_x': 'mean',
    'transformed_y': 'mean'
})

# Step 4: Turn into integers
collapsed_df['transformed_x']=collapsed_df['transformed_x'].astype(int)
collapsed_df['transformed_y']=collapsed_df['transformed_y'].astype(int)
collapsed_df['index'] = collapsed_df['index'].astype(int)
print(collapsed_df)


# Step 5: Get rid of indices (i.e. rows) where transformed_x or transformed_y are <0 or >8800
    #but keep current indices
filtered_df = collapsed_df[(collapsed_df['transformed_x'] >= 0) & (collapsed_df['transformed_x'] < 8800) &
                           (collapsed_df['transformed_y'] >= 0) & (collapsed_df['transformed_y'] < 8800)]
print(filtered_df)

# Step 6: Initialize a mask + create/store the Master Mask
all_centroids_mask = np.zeros_like(base_tif_images['base01'][:,:,0])

#extract centroids
coords = np.array(filtered_df[['transformed_x', 'transformed_y']])
#labels=np.array(filtered_df['index']) # labels of each bead

# Assign a label of 1 to all centroid locations in the mask
all_centroids_mask[coords[:,0].astype(int), coords[:,1].astype(int)] = 1 #labels#+1

# Generate sequential labels for connected regions
from skimage.measure import label
labeled_mask = (label(all_centroids_mask)).astype(np.uint32)

 # Expand labeled regions by distance of 4
from skimage.segmentation import expand_labels
expanded_labeled_mask = (expand_labels(labeled_mask, distance=4)).astype(np.uint32)
print(expanded_labeled_mask.shape)
print(f'Max label for base01 after expansion: {np.amax(expanded_labeled_mask)}')

 # Store expanded labeled mask + check for int32 data type  (which can store up to 4,294,967,295 unique values)
all_centroids_mask = expanded_labeled_mask
print(all_centroids_mask.dtype, all_centroids_mask.shape, all_centroids_mask.nbytes)


# Step 7: Visualize master mask
plt.imshow(all_centroids_mask)
plt.show()

# zoomed in version
plt.imshow(all_centroids_mask[3000:3300, 3000:3300])
plt.show()

# different color visualization
from matplotlib.colors import ListedColormap
custom_cmap = ListedColormap(['black', '#cccccc']) # #define light grey for 1 and 'black' for 0
plt.imshow(all_centroids_mask, cmap=custom_cmap, vmin=0, vmax=1, interpolation='nearest')
plt.show()
plt.imshow(all_centroids_mask[2500:3000, 2500:3000], cmap=custom_cmap, vmin=0, vmax=1, interpolation='nearest')
plt.show()

# Step 8: Save expanded labeled mask

import pickle
with open('~/dir/expanded_labeled_mask.pkl', 'wb') as file:
    pickle.dump(expanded_labeled_mask, file)