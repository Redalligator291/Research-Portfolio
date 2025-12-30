import cv2
import numpy as np
import matplotlib.pyplot as plt
from skimage import io
from skimage.filters import threshold_otsu
from skimage.filters.rank import entropy
from skimage.morphology import disk
from skimage.color import rgb2gray

# Read the .tif image 
img_path = r"C:\Users\vaithi\Tissue work\ND13326Hi_Run01_BSED_slice_0087.dm3.tif"
img = cv2.imread(img_path)  #  Read as-is

# Convert to grayscale
gray_img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

# Find the contours in the grayscale image
contours, _ = cv2.findContours(gray_img, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

# Find the largest contour (assuming it's the actual image)
largest_contour = max(contours, key=cv2.contourArea)

# Get the bounding rectangle of the largest contour
x, y, w, h = cv2.boundingRect(largest_contour)

# Crop the image to the bounding rectangle
cropped_img = img[y:y+h, x:x+w]

# Convert the cropped image to grayscale
cropped_img_gray = rgb2gray(cropped_img)

# Apply entropy filtering and Otsu thresholding
entropy_img = entropy(cropped_img_gray, disk(25))
plt.hist(entropy_img.flat, bins=100, range=(0,5))
thresh = threshold_otsu(entropy_img)
binary = entropy_img <= thresh

plt.imshow(binary)
plt.show()

# Save the segmented image as a .tif file
output_path = "segmented_image.tif"
cv2.imwrite(output_path, np.uint8(binary * 255))