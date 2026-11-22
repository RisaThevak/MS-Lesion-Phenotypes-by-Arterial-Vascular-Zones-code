import numpy as np
import scipy.ndimage as ndimage
import nibabel as nib

# Step 1: Load the MRI volume from a NIfTI file
nifti_file_path = '/export02/data/risa/Atlas_MNI152/Atlas_MNI152/ArterialAtlas_level2.nii'  # Update with the path to your NIfTI file
nifti_image = nib.load(nifti_file_path)

# Get the volume data (the MRI image with integer labels for MCA, PCA, ACA, etc.)
volume = nifti_image.get_fdata()
volume = np.squeeze(volume)

# Step 2: Define voxel resolution (example: isotropic resolution of 1mm)
affine = nifti_image.affine
voxel_spacing = np.abs(affine[:3, :3].diagonal())  # Extract voxel dimensions from affine matrix
r_x, r_y, r_z = voxel_spacing  # in mm

# Convert the 10mm boundary into voxel units
max_distance_voxels_x = 10 / r_x
max_distance_voxels_y = 10 / r_y
max_distance_voxels_z = 10 / r_z

# Step 3: Define labels for ACA, MCA, PCA, and Basilar Artery
aca_labels = [1.000076294876635, 1.9999999995343387]  # ACA
mca_labels = [3.000076294410974, 3.9999999990686774]  # MCA
pca_labels = [5.0000762939453125, 5.999999998603016]  # PCA
#basilar_labels = [7.000076293479651, 7.999999998137355]  # Basilar Artery

# Combine all target labels for easier boundary computation
target_labels = aca_labels + mca_labels + pca_labels #+ basilar_labels

# Define labels for each arterial territory
territories = {
    "ACA": aca_labels,
    "MCA": mca_labels,
    "PCA": pca_labels,
    #"Basilar": basilar_labels
}

# Step 4: Define labels to exclude (background and ventricles)
exclude_labels = [9.00007629301399, 9.999999997671694]  # Ventricles labels
exclude_mask = np.isin(volume, exclude_labels)  # True where the volume matches any of the excluded labels

# Create an empty mask for the final watershed regions
watershed_mask = np.zeros_like(volume, dtype=bool)

# Step 5: Remove small clusters (1 to 5 voxels) within each territory
min_cluster_size = 16  # Minimum cluster size to keep (1 to 5 voxels will be removed)

filtered_territories = {}

for territory, labels in territories.items():
    # Create a binary mask for the current territory
    territory_mask = np.isin(volume, labels)
    territory_mask &= ~exclude_mask  # Exclude ventricles

    # Label connected components in the territory mask
    labeled_territory, num_labels = ndimage.label(territory_mask)
    
    # Remove small clusters (1 to 5 voxels)
    filtered_territory_mask = np.zeros_like(territory_mask, dtype=bool)
    
    for label in range(1, num_labels + 1):
        # Get the size of the current connected component (cluster)
        cluster_mask = (labeled_territory == label)
        cluster_size = np.sum(cluster_mask)
        
        # If the cluster size is larger than or equal to the minimum size, keep it
        if cluster_size >= min_cluster_size:
            filtered_territory_mask |= cluster_mask
    
    # Store the filtered territory
    filtered_territories[territory] = filtered_territory_mask

# Step 6: Expand each filtered territory (dilation)
structuring_element = np.ones((int(max_distance_voxels_x), int(max_distance_voxels_y), int(max_distance_voxels_z)), dtype=bool)

expanded_territories = {}

for territory, filtered_mask in filtered_territories.items():
    # Dilate the filtered territory mask to expand it
    expanded_territories[territory] = ndimage.binary_dilation(filtered_mask, structure=structuring_element)

# Step 7: Define overlap between different territories as watershed regions
for territory1, mask1 in expanded_territories.items():
    for territory2, mask2 in expanded_territories.items():
        # Skip if both territories are the same
        if territory1 == territory2:
            continue
        
        # Define the overlap (intersection) of two expanded territories
        overlap = mask1 & mask2
        
        if np.any(overlap):  # If there is overlap
            # Label connected components in the overlap
            labeled_overlap, num_labels = ndimage.label(overlap)
            
            # Loop through each connected component (cluster)
            for label in range(1, num_labels + 1):
                # Get the size of the current connected component (cluster)
                cluster_mask = (labeled_overlap == label)
                cluster_size = np.sum(cluster_mask)
                
                # If the cluster size is larger than or equal to the minimum size, keep it; otherwise, remove it
                if cluster_size >= min_cluster_size:
                    watershed_mask |= cluster_mask

# Step 8: Mask out the background (intensity 0 in the original volume) and ventricles from the watershed mask
watershed_mask &= (volume != 0)  # Exclude background
#watershed_mask &= ~exclude_mask  # Exclude ventricles

# Step 9: Save the resulting watershed mask as a NIfTI file
watershed_img = nib.Nifti1Image(watershed_mask.astype(np.uint8), affine=nifti_image.affine)
watershed_output_path = '/export02/data/risa/Atlas_MNI152/Atlas_MNI152/watershed_regions.nii'
nib.save(watershed_img, watershed_output_path)

