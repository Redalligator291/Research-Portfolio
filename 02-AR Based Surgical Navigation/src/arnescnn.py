import os
import numpy as np
import nibabel as nib
from skimage.transform import resize
from tensorflow.keras.utils import to_categorical
from tensorflow.keras import layers, models

# Function to load medical data
def load_medical_data(image_folder, mask_folder):
    image_paths = [os.path.join(image_folder, filename) for filename in os.listdir(image_folder) if filename.endswith('.nii.gz')]
    mask_paths = [os.path.join(mask_folder, filename) for filename in os.listdir(mask_folder) if filename.endswith('.nii.gz')]

    # Assuming the data is paired (image and mask have the same filename)
    image_paths.sort()
    mask_paths.sort()

    images = [nib.load(path).get_fdata() for path in image_paths]
    masks = [nib.load(path).get_fdata() for path in mask_paths]

    return images, masks

# Function to preprocess data
def preprocess_data(images, masks):
    processed_images = np.array([resize_image(image) for image in images])
    processed_masks = np.array([resize_mask(mask) for mask in masks])

    return processed_images, processed_masks

def resize_image(image, target_shape=(128, 128, 128)):
    resized_image = resize(image, target_shape, mode='constant', anti_aliasing=True)
    return resized_image

def resize_mask(mask, target_shape=(128, 128, 128)):
    resized_mask = resize(mask, target_shape, mode='constant', anti_aliasing=False)
    return resized_mask

# Define the 3D U-Net model
def unet_3d(input_shape=(128, 128, 128, 1), num_classes=2):
    inputs = layers.Input(shape=input_shape)

    conv1 = layers.Conv3D(32, (3, 3, 3), activation='relu', padding='same')(inputs)
    conv1 = layers.Conv3D(32, (3, 3, 3), activation='relu', padding='same')(conv1)
    pool1 = layers.MaxPooling3D(pool_size=(2, 2, 2))(conv1)

    conv2 = layers.Conv3D(64, (3, 3, 3), activation='relu', padding='same')(pool1)
    conv2 = layers.Conv3D(64, (3, 3, 3), activation='relu', padding='same')(conv2)
    pool2 = layers.MaxPooling3D(pool_size=(2, 2, 2))(conv2)

    conv3 = layers.Conv3D(128, (3, 3, 3), activation='relu', padding='same')(pool2)
    conv3 = layers.Conv3D(128, (3, 3, 3), activation='relu', padding='same')(conv3)
    pool3 = layers.MaxPooling3D(pool_size=(2, 2, 2))(conv3)

    conv4 = layers.Conv3D(256, (3, 3, 3), activation='relu', padding='same')(pool3)
    conv4 = layers.Conv3D(256, (3, 3, 3), activation='relu', padding='same')(conv4)

    up5 = layers.Concatenate()([layers.UpSampling3D(size=(2, 2, 2))(conv4), conv3])
    conv5 = layers.Conv3D(128, (3, 3, 3), activation='relu', padding='same')(up5)
    conv5 = layers.Conv3D(128, (3, 3, 3), activation='relu', padding='same')(conv5)

    up6 = layers.Concatenate()([layers.UpSampling3D(size=(2, 2, 2))(conv5), conv2])
    conv6 = layers.Conv3D(64, (3, 3, 3), activation='relu', padding='same')(up6)
    conv6 = layers.Conv3D(64, (3, 3, 3), activation='relu', padding='same')(conv6)

    up7 = layers.Concatenate()([layers.UpSampling3D(size=(2, 2, 2))(conv6), conv1])
    conv7 = layers.Conv3D(32, (3, 3, 3), activation='relu', padding='same')(up7)
    conv7 = layers.Conv3D(32, (3, 3, 3), activation='relu', padding='same')(conv7)

    outputs = layers.Conv3D(num_classes, (1, 1, 1), activation='softmax')(conv7)

    model = models.Model(inputs=inputs, outputs=outputs)
    model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])

    return model

# Define folder paths for training data
training_data_folder = r'C:\Users\Kavin Ramadoss\Images for CNN\training_data'

# Load and preprocess training data
training_images, training_masks = load_medical_data(os.path.join(training_data_folder, 'images'), os.path.join(training_data_folder, 'masks'))
processed_training_images, processed_training_masks = preprocess_data(training_images, training_masks)

# Convert masks to categorical format (assuming binary segmentation)
num_classes = 2  # Assuming binary segmentation
y_train_categorical = to_categorical(processed_training_masks, num_classes=num_classes)

# Define the 3D U-Net model
model = unet_3d(input_shape=processed_training_images.shape[1:], num_classes=num_classes)
model.summary()

# Train the model on the training data
model.fit(processed_training_images, y_train_categorical, epochs=50, batch_size=1)

# Save the trained model
model.save(r'C:\Users\Kavin Ramadoss\Images for CNN\unet_3d_model.h5')

