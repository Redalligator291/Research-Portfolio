import os
import nibabel as nib  
import numpy as np
import keras
from keras.models import Model
from keras.layers import Input, Conv3D, MaxPooling3D, UpSampling3D, Concatenate

# Set path variables
train_img_dir = 'training_images/images/'
train_mask_dir = 'training_images/masks/'

# Load training data 
train_imgs = []
for img_file in os.listdir(train_img_dir):
    img_path = os.path.join(train_img_dir, img_file)  
    img = nib.load(img_path)
    train_imgs.append(img)
    
train_masks = []
for mask_file in os.listdir(train_mask_dir):  
    mask_path = os.path.join(train_mask_dir, mask_file)
    mask = nib.load(mask_path)
    train_masks.append(mask)

assert len(train_imgs) == len(train_masks)

# Define model
img_input = Input(shape=(None, None, None, 1))
mask_input = Input(shape=(None, None, None, 1))

x = Conv3D(32, (3,3,3), activation='relu')(img_input) 
x = Conv3D(32, (3,3,3), activation='relu')(x)
x = MaxPooling3D((2,2,2))(x)  


x = Conv3D(64, (3,3,3), activation='relu')(x) 
x = Conv3D(64, (3,3,3), activation='relu')(x)
x = MaxPooling3D((2,2,2))(x)

x = Conv3D(128, (3,3,3), activation='relu')(x)
x = Conv3D(128, (3,3,3), activation='relu')(x)
x = MaxPooling3D((2,2,2))(x)

x = Conv3D(256, (3,3,3), activation='relu')(x)
x = Conv3D(256, (3,3,3), activation='relu')(x)
x = MaxPooling3D((2,2,2))(x)

# Decoder pathway
x = UpSampling3D((2,2,2))(x)
x = Conv3D(128, (3,3,3), activation='relu')(x)
x = Concatenate()([x, conv3_2]) # Skip connection

x = UpSampling3D((2,2,2))(x) 
x = Conv3D(64, (3,3,3), activation='relu')(x)
x = Concatenate()([x, conv2_2]) # Skip connection

x = UpSampling3D((2,2,2))(x)
x = Conv3D(32, (3,3,3), activation='relu')(x) 
x = Concatenate()([x, conv1_2]) # Skip connection

output = Conv3D(num_labels, (1,1,1), activation='softmax')(x)   

model = Model([img_input, mask_input], output)  

# Compile and train model
model.compile(optimizer='adam', loss='categorical_crossentropy')
model.fit(
    [train_imgs, train_masks],
    train_masks,
    epochs=50,  
    batch_size=8)