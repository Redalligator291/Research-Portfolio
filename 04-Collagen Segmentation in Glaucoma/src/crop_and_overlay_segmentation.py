import cv2
import numpy as np
import matplotlib.pyplot as plt
from skimage import io
from skimage.filters import threshold_otsu
from skimage.filters.rank import entropy
from skimage.morphology import disk
from imutils.perspective import four_point_transform

# This line of code reads the image
img = r"C:\Users\vaithi\Tissue work\Data\ND13326Hi_Run01_BSED_slice_0087.dm3.tif"

# Reads the image in using OpenCV
image = cv2.imread(img)

# This shows the original image using a python library called matplotlib
plt.imshow(image)
plt.title("Original Image")
plt.show()

# These next lines crop the images using Gaussian Blur and other functions
gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
blur = cv2.GaussianBlur(gray, (5, 5), 0)
thresh = cv2.threshold(blur, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)[1]
contours, hierarchy = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
cnts = sorted(contours, key=cv2.contourArea, reverse=True)
cnt = cnts[0]
x, y, w, h = cv2.boundingRect(cnt)
crop = image[y:y+h, x:x+w]

# This shows the cropped image
plt.imshow(crop)
plt.title("Cropped Image")
plt.show()

seg_img=r"C:\Users\vaithi\Tissue work\segmented_image.tif"
org_img=r"C:\Users\vaithi\Tissue work\Data\ND13326Hi_Run01_BSED_slice_0087.dm3.tif"
# Read the original image again
original_image = cv2.imread(org_img)
segmented_image = cv2.imread(seg_img)

# Get the dimensions of the segmented image
segmented_image_height, segmented_image_width = segmented_image.shape[:2]

# Create a new image with the same dimensions as the original image
result = original_image.copy()

# Paste the segmented image onto the cropped area of the original image
result[y:y+h, x:x+w] = segmented_image

plt.imshow(result)
plt.show()
# Save the result as a new image
output_path = "segmented_img_with_border.tif"
cv2.imwrite(output_path, result)