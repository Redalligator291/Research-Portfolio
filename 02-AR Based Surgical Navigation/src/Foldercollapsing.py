import os
import shutil

def flatten_directory(src_folder, dest_folder):
    for root, dirs, files in os.walk(src_folder):
        for file in files:
            src_path = os.path.join(root, file)
            dest_path = os.path.join(dest_folder, file)
            shutil.move(src_path, dest_path)

# Specify source and destination folders
src_folder = r"C:\Users\Kavin Ramadoss\Images for CNN\Testing Images for CNN"
dest_folder = r"C:\Users\Kavin Ramadoss\Images for CNN\Testing Images for CNN"

# Call the function to flatten the directory
flatten_directory(src_folder, dest_folder)
