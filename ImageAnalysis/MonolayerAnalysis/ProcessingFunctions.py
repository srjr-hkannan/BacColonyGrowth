# Dependencies
import os
import cv2
import numpy as np
import pandas as pd
import datetime
import json
import scipy.io as sio
import seaborn as sns

import tifffile
import glob
import imageio
from PIL import Image
from matplotlib import pyplot as plt
from scipy.ndimage import gaussian_filter
import skg

def processImage(stackFile):
    
    #Sum Z-Projection
    stackImage = tifffile.imread(stackFile)
    projectionImage = stackImage.sum(axis=0)
    projectionImage[projectionImage < 0] = 0
    projectionImage = 255*projectionImage/np.amax(projectionImage)
    
    # Convert to 8-bit
    projectionImage = (projectionImage).astype(np.uint8)
    smoothImage = cv2.GaussianBlur(projectionImage, (0,0), sigmaX=1, sigmaY=1)

    _, threshImage = cv2.threshold(projectionImage, 20, 255, cv2.THRESH_BINARY)
    cont = cv2.findContours(threshImage, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    # get largest contours
    contours = cv2.findContours(threshImage, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    contours = contours[0] if len(contours) == 2 else contours[1]
    big_contour = max(contours, key=cv2.contourArea)

    # smooth contour
    peri = cv2.arcLength(big_contour, True)
    big_contour = cv2.approxPolyDP(big_contour, 0.001 * peri, True)

    # draw white filled contour on black background
    contour_img = np.zeros_like(threshImage)
    cv2.drawContours(contour_img, [big_contour], 0, 255, -1)

    # apply dilate to connect the white areas in the alpha channel
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (15,15))
    dilate = cv2.morphologyEx(contour_img, cv2.MORPH_DILATE, kernel)

    # make edge outline
    edge = cv2.Canny(dilate, 0, 200)

    # thicken edge
    #edge = cv2.GaussianBlur(edge, (0,0), sigmaX=0.3, sigmaY=0.3)
    _, outlineImage = cv2.threshold(edge, 100, 255, cv2.THRESH_BINARY)  
    
    return projectionImage, smoothImage, outlineImage

def offsetCalc(distance, m):
    return(distance*np.sqrt(1+m*m))
    