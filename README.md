# Liver and Tumor Segmentation Pipeline

This project implements a MATLAB-based pipeline for segmenting the liver and tumors from NIfTI medical images. It includes stages for data loading, extensive preprocessing, segmentation, evaluation, and visualization. All image processing and segmentation tasks are performed using classical image processing algorithms and techniques; no machine learning models are employed.

## Key Features & Pipeline Stages

The pipeline processes 3D medical image volumes through several stages:

1.  **Data Loading**:
    *   Loads NIfTI image and label files (e.g., `*.nii.gz`) using functions like [`data_io.loadNiiVolumes`](+data_io/loadNiiVolumes.m).
    *   Handles `.gz` decompression.

2.  **Preprocessing**:
    *   **ROI-based Normalization**: Normalizes volume intensity based on a defined Region of Interest using [`preprocessing.normalizeVolumeByRoi`](+preprocessing/normalizeVolumeByRoi.m).
    *   **Histogram Analysis**:
        *   Calculates volume histograms ([`histogram_tools.calculateVolumeHistogram`](+histogram_tools/calculateVolumeHistogram.m)).
        *   Groups histogram bins ([`histogram_tools.groupHistogram`](+histogram_tools/groupHistogram.m)).
        *   Detects intensity bands within histograms ([`histogram_tools.detectHistogramBand`](+histogram_tools/detectHistogramBand.m)).
        *   Removes histogram outliers ([`histogram_tools.removeHistogramOutliers`](+histogram_tools/removeHistogramOutliers.m)).
    *   **Contrast Stretching**: Adjusts volume contrast, for example, using [`preprocessing.stretchVolumeContrast`](+preprocessing/stretchVolumeContrast.m).
    *   **Artifact/Structure Removal**: Removes high-intensity structures (e.g., bones) using functions like [`preprocessing.removeHighIntensityArtifacts`](+preprocessing/removeHighIntensityArtifacts.m) and [`preprocessing.createStructureMask`](+preprocessing/createStructureMask.m).
    *   **Histogram Matching**: Matches volume histograms to a reference slice ([`preprocessing.matchHistogramToReference`](+preprocessing/matchHistogramToReference.m)).
    *   **Bilateral Filtering**: Applies bilateral filters for noise reduction while preserving edges (e.g., [`preprocessing.applyBilateralFilter`](+preprocessing/applyBilateralFilter.m)).
    *   **Mask Application**: Applies masks to isolate regions of interest ([`preprocessing.applyMask`](+preprocessing/applyMask.m)).

3.  **Segmentation**:
    *   **Liver Segmentation**:
        *   Initial 2D segmentation based on intensity thresholds and morphological operations (e.g., [`segmentation.segmentLiverByThreshold2D`](run_liver_analysis_pipeline.m)).
        *   2D refinement using morphological operations ([`segmentation.refineLiverMask2D`](run_liver_analysis_pipeline.m)).
        *   3D refinement, often by selecting the largest connected component ([`segmentation.refineMask3DLargestComponent`](run_liver_analysis_pipeline.m)).
        *   Post-processing for connectivity ([`segmentation.postProcessLiverMaskConnectivity`](run_liver_analysis_pipeline.m)).
    *   **Tumor Segmentation**: Segments tumors within the liver region, typically involving thresholding and morphological operations (e.g., [`segmentation.segmentTumors`](run_liver_analysis_pipeline.m)).

4.  **Evaluation**:
    *   Calculates segmentation performance metrics such as Dice score and Recall for both liver and tumor against ground truth labels ([`evaluation.calculateSegmentationMetrics`](+evaluation/calculateSegmentationMetrics.m)).
    *   Displays these metrics in the command window ([`evaluation.displayMetrics`](+evaluation/displayMetrics.m)).

5.  **Visualization**:
    *   Displays 2D slices of volumes at various processing stages, often with overlays ([`visualization.viewVolumeSlices`](run_liver_analysis_pipeline.m)).
    *   Plots histograms and grouped histograms with detected bands ([`visualization.plotHistogramComparison`](run_liver_analysis_pipeline.m), [`visualization.plotGroupedHistogramWithBand`](run_liver_analysis_pipeline.m)).
    *   Provides a GUI for comparing segmentation results with ground truth slice by slice ([`visualization.displaySegmentationComparisonGUI`](+visualization/displaySegmentationComparisonGUI.m)).
    *   Renders 3D surfaces of the segmented liver, tumors, and optionally bones ([`visualization.render3DSegmentation`](run_liver_analysis_pipeline.m)).

## Directory Structure

The project is organized using MATLAB packages:

*   `+data_io/`: Functions for data input/output.
*   `+evaluation/`: Scripts for evaluating segmentation performance.
*   `+histogram_tools/`: Utilities for histogram manipulation.
*   `+morph_ops/`: Functions for morphological operations.
*   `+preprocessing/`: Scripts for image preprocessing steps.
*   `+segmentation/`: Algorithms for liver and tumor segmentation.
*   `+visualization/`: Tools for displaying images, plots, and 3D renderings.
*   `data/`: Intended for storing NIfTI image and label files (Note: this directory is ignored by Git as per [`.gitignore`](.gitignore)).
*   `old/`: Contains older versions or auxiliary scripts.
*   `run_liver_analysis_pipeline.m`: Main script to execute the entire pipeline.

## How to Run

1.  Ensure you have MATLAB installed with the Image Processing Toolbox. You will also need a NIfTI I/O library (the project uses `load_nii`, which is common in the neuroimaging community).
2.  Place your NIfTI image and label files (e.g., `liver_XX.nii.gz`) into a `data/imagesTr` and `data/labelsTr` subdirectory respectively, or update the paths in the main configuration section of the script.
3.  Open MATLAB and navigate to the project's root directory.
4.  Run the main pipeline script: `run_liver_analysis_pipeline.m`.

The script [`run_liver_analysis_pipeline.m`](run_liver_analysis_pipeline.m) contains a `config` structure at the beginning where parameters like file paths, algorithm settings, and verbosity can be adjusted.

## Results and Methodology

This project achieves liver and tumor segmentation using a sequence of **classical image processing techniques**. These include:

*   Intensity-based thresholding
*   Histogram analysis and manipulation (normalization, stretching, matching)
*   Morphological operations (erosion, dilation, opening, closing, connected component analysis)
*   Spatial filtering (e.g., bilateral filtering)

**No machine learning or deep learning models are used in this pipeline.** The segmentation is entirely driven by algorithmic processing of image intensities and spatial relationships. The performance of the pipeline is evaluated using standard metrics like Dice coefficient and Recall, which are reported at the end of the script execution. The specific results will vary depending on the input data and the fine-tuning of parameters within the `config` structure.

For reference, the algorithm with the current configuration was tested on the following dataset:
*   **Dataset**: [Task03_LiverTumor - Medical Decathlon Challenge](http://medicaldecathlon.com/)

#### Results for Task03_LiverTumor_liver_80

```
--- Segmentation Metrics ---
Dice Score Liver: 0.9496
Dice Score Tumor: 0.8680
Recall Liver:     0.9533
Recall Tumor:     0.8288
---------------------------
```


## Dependencies

*   MATLAB
*   Image Processing Toolbox (for functions like `imhist`, `imfill`, `strel`, morphological operations, etc.)
*   NIfTI file I/O tools (e.g., "Tools for NIfTI and ANALYZE image" by Jimmy Shen, which provides `load_nii`).