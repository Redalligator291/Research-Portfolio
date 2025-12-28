# Collagen Segmenatation in Glaucoma 

This folder contains  image-processing scripts used for cropping and entropy-based segmentation of tissue images.

## Files

- `crop_and_entropy_segment.py` — Crop an image to the largest contour, compute entropy texture, apply Otsu thresholding, and save a segmentation mask (`segmented_img.tif`).
- `crop_and_entropy_segmenter.py` — Similar to `crop_and_entropy_segment.py`; produces a cropped image (`cropped_img.tif`) and an entropy-based segmentation (`segmented_img.tif`).
- `crop_and_overlay_segmentation.py` — Crop the largest region, load an externally computed segmented image, and overlay it back onto the original image; saves `segmented_img_with_border.tif`.
- `tissue_contour_segmenter.py` — Simple contour-based tissue vs fluid classification using Canny edge detection and contour aspect ratio; displays the result in a window.

## Quick Notes

- Input paths in scripts are currently hard-coded. Update `img` / `img_path` variables to point to your local data or modify scripts to accept command-line arguments.
- Scripts use OpenCV, scikit-image, NumPy, Matplotlib, and imutils. Example minimal requirements:

```
opencv-python
scikit-image
numpy
matplotlib
imutils
```

- Output files are TIFF images saved in the working directory by default.

## Suggested usage

1. Edit the image path at the top of the script or add argument parsing.
2. Run the script, for example:

```bash
python crop_and_entropy_segment.py
```

3. Inspect outputs (`cropped_img.tif`, `segmented_img.tif`, `segmented_img_with_border.tif`).

## Next steps (optional)

- Replace hard-coded paths with CLI args or a small runner script.
- Add a small `requirements.txt` or `pyproject.toml`.
- Convert image I/O to use relative paths or a `data/` folder.
