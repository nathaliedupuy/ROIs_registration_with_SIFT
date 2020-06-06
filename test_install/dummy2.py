import cv2

cv2.__version__

print('All good so far with open cv!')


# Initiate SIFT detector
sift = cv2.xfeatures2d.SIFT_create()

print('Creation SIFT object succesful!')