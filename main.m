clear; clc; clf;
disp("Code Run")

% Option to visualize slices
visualizeSlicesFlag = false;


% Load the slices from files
[trainVolume, labelVolume] = loadNiiFile(...
    'data/imagesTr/liver_80.nii.gz', ...
    'data/labelsTr/liver_80.nii.gz'  ...
);


% Visualize the slices
if visualizeSlicesFlag
    visualizeSlices(trainVolume, labelVolume);
end

%%
% Define the region of interest (ROI) parameters for cross normalization
roiParams.x = 42;
roiParams.y = 13;
roiParams.r = 10;

% Maximum value for uint16
maxValue = 65536;

% ----------- START REFACTOR FOLLOWING CODE ----------- %

% Number of bins for histogram
nBins = 65536; % Number of bins for histogram
[meanValue, normalizedSlice] = normalizingSlices(trainVolume, roiParams, maxValue);
[hMean, hMean_clean] = histogramOnAllSlices(normalizedSlice, nBins);

%%
fprintf('Mean value in ROI: %.2f\n', meanValue);

% Visualize mean histogram of all slices
plotHistograms(hMean, hMean_clean);

% Group histogram data
groupSize = 2000;

[grouped_hMean, grouped_hMean_clean, binCenters] = groupHistogramData(hMean, hMean_clean, groupSize, nBins);
[lowerIntensity, upperIntensity] = bandDetection(grouped_hMean_clean);

plotGroupedHistograms(binCenters, grouped_hMean, grouped_hMean_clean, nBins, lowerIntensity, upperIntensity);

%%

% ------------ STOP REFACTOR FOLLOWING CODE ----------- %

% Stretch the slices using the lower and upper intensity values
stretchedSlice = stretchSlices(normalizedSlice, lowerIntensity, upperIntensity, 4);

slice_idx = 154;

% Remove high value (in this case mainly areas representing the bones)
mask = stretchedSlice(:,:,slice_idx) > maxValue / 2.5;
mask = imdilate(double(mask), strel("disk", 4)); % Dilation to remove small holes
mask = imerode(mask, strel("disk", 2)); % Erosion to remove small holes
mask = imfill(mask, "holes"); % Fill holes
mask = imfilter(mask, fspecial("gaussian", 10)); % Gaussian filter to smooth the mask
% ------------ START OF TEST CODE ----------- %

tempSlice = uint16(double(stretchedSlice(:,:,slice_idx)) .* double(1 - mask));
stretchedSlice(:,:,slice_idx) = tempSlice;

% Visualizzazione della slice modificata
figure;
imshow(stretchedSlice(:,:,slice_idx));
doubleStretchedSlice = stretchSlices(stretchedSlice, 0.15, 0.3, 0.4);

figure;
% Visualizzazione della slice modificata
imshow(doubleStretchedSlice(:,:,slice_idx));
se_open = strel("disk", 8);

figure;
imshow(imopen(imfill(imfilter(doubleStretchedSlice(:,:,slice_idx), fspecial("disk", 3)) > nBins / 3, "holes"), se_open))

if visualizeSlicesFlag
    visualizeSlices(normalizedSlice, stretchedSlice);
end


disp("Code End")