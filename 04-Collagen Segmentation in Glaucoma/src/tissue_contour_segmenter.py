import cv2
import numpy as np

# Read the image
image = cv2.imread(r'C:\Users\vaithi\Tissue work\tissue image better.png')

# Convert the image to grayscale
gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

# Apply Canny edge detection
edges = cv2.Canny(gray_image, 100, 100)

# Find contours of the edges
contours, hierarchy = cv2.findContours(edges, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

# Draw contours on the original image
for contour in contours:
   # Get the contour's bounding box
   x, y, w, h = cv2.boundingRect(contour)

   # Check if the contour is a tissue or fluid contour
   if w > h:
       # Tissue contour is wider than taller
       cv2.drawContours(image, [contour], 0, (0, 0, 0), 2)
   else:
       # Fluid contour is taller than wider
       cv2.drawContours(image, [contour], 0, (255, 255, 255), -1)

# Display the result
cv2.imshow('Result', image)
cv2.waitKey(0)
cv2.destroyAllWindows()