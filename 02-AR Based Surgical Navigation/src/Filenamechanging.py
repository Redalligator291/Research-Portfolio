import os

def rename_files(folder_path, prefix="test_image", extension=".nii.gz"):
    files = os.listdir(folder_path)

    for i, file_name in enumerate(files, start=1):
        old_path = os.path.join(folder_path, file_name)
        new_name = f"{prefix}{i}{extension}"
        new_path = os.path.join(folder_path, new_name)

        os.rename(old_path, new_path)
        print(f"Renamed: {file_name} to {new_name}")

# Example usage
folder_to_rename = r"C:\Users\Kavin Ramadoss\Images for CNN\Testing Images for CNN"
rename_files(folder_to_rename)
