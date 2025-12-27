import os
import nibabel as nib  
import re
from skimage.filters import threshold_yen
from skimage.morphology import remove_small_objects

dataset_path = r'C:\Users\Kavin Ramadoss\Images for CNN\test_data\test_images'
output_path = r'C:\Users\Kavin Ramadoss\Images for CNN\test_data\test_masks'

file_names = os.listdir(dataset_path)
file_names.sort(key=lambda f: int(re.sub('\D', '', f)))

for i, image_file in enumerate(file_names):
    image_path = os.path.join(dataset_path, image_file)  
    image = nib.load(image_path)
    image_data = image.get_fdata()

    thresh = threshold_yen(image_data)
    bw = image_data > thresh

    cleared = remove_small_objects(bw, 750)

    mask_file = f'mask_{i+1}.nii'
    mask_nifti = nib.Nifti1Image(cleared.astype(int), image.affine, image.header)
    nib.save(mask_nifti , os.path.join(output_path, mask_file))