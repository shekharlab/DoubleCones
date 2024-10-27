from cellpose import models
import tifffile as tiff
import numpy as np
import sys
import matplotlib.pyplot as plt
from skimage import measure
from skimage import io
from skimage.measure import label, regionprops
import os
from cellpose import utils
import cv2

# Function to load the z-stack
def load_z_stack(file_path):
    try:
        z_stack = tiff.imread(file_path)
        # Check if the loaded image is a 3D array
        # if z_stack.ndim != 3:
        # raise ValueError("The loaded image is not a 3D z-stack.")
        print(f"Z-stack loaded successfully with shape {z_stack.shape}")
        return z_stack
    except Exception as e:
        print(f"Error loading z-stack: {e}")
        sys.exit(1)  # Exit if the image can't be loaded

def quick_segment(image, **kwargs):
    model = models.Cellpose(gpu=True, model_type='nuclei')
    masks, flows, styles, diams = model.eval(image, **kwargs)
    plot_two_images(image, masks)
    return masks

# Function to run Cellpose segmentation
def run_cellpose_segmentation(z_stack, gpu = True, **kwargs):
    # try:
        # Check if z-stack is non-empty
        if z_stack is None or not np.any(z_stack):
            raise ValueError("The z-stack is empty or invalid.")
        
        # Initialize the Cellpose model
        model = models.Cellpose(gpu=gpu, model_type='nuclei')
        
        # Run 3D segmentation
        masks, flows, styles, diams = model.eval(z_stack, **kwargs)
                                                 # diameter=20,
                                                 # channels=channels,
                                                 # do_3D=True)
                                                 # cellprob_threshold= -2, 
                                                 # normalize = True, 
                                                 # rescale = False
                                                 # flow_threshold=0.3 # unused in 3D
                                                 # mask_threshold=0.1 # doesn't exist?
                                                # )
        
        # Check if any masks were returned
        if masks is None or np.max(masks) == 0:
            raise ValueError("No segmentation masks were created.")
        
        print(f"Segmentation successful. Found {np.max(masks)} objects.")
        return masks
    # except Exception as e:
    #     print(f"Error during segmentation: {e}")
    #     sys.exit(1)  # Exit if segmentation fails

# Function to save the segmented z-stack
def save_segmented_z_stack(masks, output_file):
    try:
        tiff.imwrite(output_file, masks.astype('uint16'))
        print(f"Segmented z-stack saved successfully as {output_file}")
    except Exception as e:
        print(f"Error saving segmented z-stack: {e}")
        sys.exit(1)  # Exit if saving fails

def run_and_save_cellpose(z_stack_file, gpu = True, normalize_each = False, **kwargs):
    
    filename_parts = []
    
    # Iterate through the kwargs and create parts for the filename
    for key, value in kwargs.items():
        filename_parts.append(f"{key}_{value}")

    # Join the parts with underscores
    filename = "-".join(filename_parts)

    # Add a file extension (e.g., .txt)
    base_name = z_stack_file.split('.')[0]
    output_file = base_name + "-" + f"{filename}.tif"

    # Step 1: Load the z-stack
    z_stack = load_z_stack(z_stack_file)
    if(len(z_stack.shape) == 4):
        z_stack = z_stack[:,0,:,:] # DAPI channel

    # If normalize
    if(normalize_each):
        
        # Normalize each slice in the z-stack
        # If z-stack is 3D (z, height, width), we iterate over z
        # If z-stack is 4D (z, height, width, channels), you need to normalize per channel
        normalized_stack = np.zeros_like(z_stack)  # Create an empty array to store normalized images
        
        for z in range(z_stack.shape[0]):  # Iterate through each z-slice
            normalized_stack[z] = normalize_1_99(z_stack[z])

        # Check that it worked
        plot_two_images(z_stack[0], normalized_stack[0])
        breakpoint()
        z_stack = normalized_stack

    # Step 2: Run Cellpose segmentation
    segmented_masks = run_cellpose_segmentation(z_stack, gpu = gpu, **kwargs)
    
    # Step 3: Save the segmented result
    print(f"Writing output file to: {output_file}")
    save_segmented_z_stack(segmented_masks, output_file)
    return output_file

# Normalize each image slice between 0 and 1
def normalize_image(img):
    min_val = np.min(img)
    max_val = np.max(img)
    if max_val > min_val:
        return (img - min_val) / (max_val - min_val)  # Normalize to [0, 1]
    else:
        return img  # If max == min, the image is constant, no need to normalize

def normalize_1_99(image):
    """
    Normalize the image based on the 1st and 99th percentiles.

    Args:
    - image (numpy array): The input image to be normalized.

    Returns:
    - normalized_image (numpy array): The normalized image.
    """
    # Calculate the 1st and 99th percentiles
    p1 = np.percentile(image, 1)
    p99 = np.percentile(image, 99)

    # Clip the image
    clipped_image = np.clip(image, p1, p99)

    # Normalize to the range 0-255 or 0-1
    normalized_image = (clipped_image - p1) / (p99 - p1)  # Scale to [0, 1]
    # normalized_image = normalized_image.astype(np.float32)
    normalized_image = (normalized_image * 255).astype(np.uint8)  # Scale to [0, 255] if needed

    return normalized_image
    
def plot_two_images(image1, image2, title1=None, title2=None):
    """
    Function to plot two images side by side using matplotlib.

    Args:
    - image1 (numpy array): The first image to be plotted.
    - image2 (numpy array): The second image to be plotted.
    - title1 (str): Optional. Title for the first image.
    - title2 (str): Optional. Title for the second image.

    Returns:
    None
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))  # Create 1x2 subplot
    
    # Plot first image
    axes[0].imshow(image1, cmap='gray', vmin = 0, vmax = 5000)
    axes[0].axis('off')  # Turn off axis for the first image
    if title1:
        axes[0].set_title(title1)
    
    # Plot second image
    axes[1].imshow(image2, cmap='gray')
    axes[1].axis('off')  # Turn off axis for the second image
    if title2:
        axes[1].set_title(title2)
    
    plt.tight_layout()  # Adjust layout so titles don't overlap
    plt.show()  # Display the plot


def check_masks(z_stack_file, mask_file, vmin = None, vmax = None, binarize_masks = False):
    
    # Load the z-stack and masks
    if isinstance(z_stack_file, str):
        z_stack = load_z_stack(z_stack_file)
        if(len(z_stack.shape) == 4):
            z_stack = z_stack[:,0,:,:] # DAPI channel
    else: 
        z_stack = z_stack_file

    if isinstance(mask_file, str):
        masks = tiff.imread(mask_file)  # Load the corresponding masks from Cellpose
    else: 
        masks = mask_file
    
    # Select the number of slices to view
    num_slices_to_view = 10  # Adjust as needed
    total_slices = z_stack.shape[0]  # Get total number of slices
    
    # Generate a list of slice indices to view
    slice_indices = np.linspace(0, total_slices - 1, num_slices_to_view, dtype=int)
    
    # Set up the plot
    fig, axes = plt.subplots(num_slices_to_view, 2, figsize=(10, 4 * num_slices_to_view))
    
    global_min = min(np.min(img) for img in z_stack)
    global_max = max(np.max(img) for img in z_stack)
    if vmin is None:
        vmin = global_min
    if vmax is None:
        vmax = global_max

    # if binarize_masks:
    #     masks2 = np.zeros_like(masks)
    #     masks2[np.where(masks != 0)] = 1
    #     masks = masks2.copy()
    
    # Loop through selected slices and display the original image and mask
    for i, slice_index in enumerate(slice_indices):
        # Display original z-stack slice
        axes[i, 0].imshow(z_stack[slice_index], cmap='gray', vmin = vmin, vmax = vmax)
        axes[i, 0].set_title(f'Original Slice {slice_index}')
        # axes[i, 0].axis('off')  # Hide axes for clarity
    
        # Display corresponding mask slice
        gray_image = z_stack[slice_index]
        scaled_image = cv2.normalize(gray_image, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)
        rgb_image = cv2.cvtColor(scaled_image, cv2.COLOR_GRAY2RGB)
        contour_image = rgb_image.copy()
        outlines = utils.masks_to_outlines(masks[slice_index])
        outX, outY = np.nonzero(outlines)
        contour_image[outX,outY] = np.array([255,75,75]) # Red
        # contour_image = cv2.cvtColor(z_stack[slice_index], cv2.COLOR_GRAY2RGB)
        # outline = utils.masks_to_outlines(masks[slice_index])
        # contour_image[np.where(outline[slice_index] != 0)] = [255,0,0]
        axes[i, 1].imshow(contour_image)  # cmap='jet', alpha=0.5)
        axes[i, 1].set_title(f'Mask Slice {slice_index}')
        # axes[i, 1].axis('off')  # Hide axes for clarity
    
    # Show the plots
    plt.tight_layout()
    plt.show()


def histogram(data_col, nbins = 10):
    
    # Create a histogram of the 'values' column
    plt.figure(figsize=(10, 6))  # Set the figure size
    data_col.hist(bins=nbins, color='skyblue', edgecolor='black')  # Create histogram
    
    # Customize the plot
    plt.title('Histogram of Values')
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.grid(axis='y', alpha=0.75)  # Add gridlines for better readability
    # plt.xticks(range(1, 8))  # Set x-ticks based on the data range
    
    # Show the plot
    plt.show()

def scatterplot(x, y, title="Scatter Plot", xlabel="X-axis", ylabel="Y-axis"):
    """
    Creates a basic scatter plot.

    Parameters:
    x (list or array-like): Data for the X-axis.
    y (list or array-like): Data for the Y-axis.
    title (str): Title of the scatter plot (default is "Scatter Plot").
    xlabel (str): Label for the X-axis (default is "X-axis").
    ylabel (str): Label for the Y-axis (default is "Y-axis").
    """
    plt.figure(figsize=(8, 6))
    plt.scatter(x, y, color='blue', alpha=0.7)
    
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    plt.grid(True)
    plt.show()

def plot_image(image, title="Image"):
    """
    Plots an image using matplotlib.

    Parameters:
    image (ndarray): The image data (can be grayscale or RGB).
    title (str): The title of the image plot.
    """
    plt.imshow(image, cmap='gray')
    plt.title(title)
    # plt.axis('off')  # Hide axes ticks and labels
    plt.show()