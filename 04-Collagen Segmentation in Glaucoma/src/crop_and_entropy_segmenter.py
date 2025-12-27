import cv2
import numpy as np
import matplotlib.pyplot as plt
from skimage import io
from skimage.filters import threshold_otsu
from skimage.filters.rank import entropy
from skimage.morphology import disk
from imutils. perspective import four_point_transform

#This line of code reads the image
img = r"C:\Users\vaithi\Tissue work\Data\ND13326Hi_Run01_BSED_slice_0087.dm3.tif"

#Reads the image in using OpenCV
image = cv2.imread(img)

#This shows the original image using a python library called matplotlib
plt.imshow(image)
plt.title("Original Image")
plt.show()

#These next lines crop the images using Gaussin Blur and other functions
gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
blur = cv2.GaussianBlur(gray, (5,5), 0)
thresh = cv2.threshold(blur, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)[1]
contours,hierarchy = cv2.findContours(thresh,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
cnts = sorted(contours, key=cv2.contourArea, reverse=True)
cnt = cnts[0]
x,y,w,h = cv2.boundingRect(cnt)
crop = image[y:y+h,x:x+w]

#This shows the cropped image
plt.imshow(crop)
plt.title("Cropped Image")
plt.show()

#These lines of code below save the cropped image
output_path = "cropped_img.tif"
cv2.imwrite(output_path, crop)

#Reads the cropped image using OpenCV
img_path = r"C:\Users\vaithi\Tissue work\cropped_img.tif"
img = cv2.imread(img_path, cv2.IMREAD_GRAYSCALE)
#Runs the entropy function on the image and performs a texture anaylsis
entropy_img = entropy(img, disk(100))
plt.hist(entropy_img.flat, bins=35, range=(0,5))
thresh = threshold_otsu(entropy_img)
binary = entropy_img <= thresh
inverted_binary = ~binary

#Shows the fully segmented image
plt.imshow(inverted_binary, cmap='binary')
plt.title("Fully Segmented Image")
plt.show()

#Saves the segmented image as a .tif file
output_path = "segmented_img.tif"
cv2.imwrite(output_path, np.uint8(binary * 255))