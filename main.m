clear; clc; clf;
disp("Code Run")

% Option to visualize slices
visualizeSlicesFlag = false;




% Load the slices from files
[trainVolume, labelVolume] = loadNiiFile('/home/santal/Documents/University/1-anno/ImageProcessing/LiverProject/data/imagesTr/liver_80.nii.gz', '/home/santal/Documents/University/1-anno/ImageProcessing/LiverProject/data/labelsTr/liver_80.nii.gz');


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
stretchedSlice = stretchSlices(normalizedSlice, lowerIntensity, upperIntensity, 5);

slice_idx = 154;

% Remove high value (in this case mainly areas representing the bones)
mask = stretchedSlice(:,:,slice_idx) > maxValue / 2.8;
mask = imdilate(double(mask), strel("disk", 3)); % Dilation to remove small holes
mask = imfill(mask, "holes"); % Fill holes
mask = imfilter(mask, fspecial("gaussian", 10)); % Gaussian filter to smooth the mask
%%
% ------------ START OF TEST CODE ----------- %

tempSlice = uint16(double(stretchedSlice(:,:,slice_idx)) .* double(1 - mask));
stretchedSlice(:,:,slice_idx) = tempSlice;

% Visualizzazione della slice modificata
figure;
imshow(stretchedSlice(:,:,slice_idx));
stretchedSlice = stretchSlices(stretchedSlice, 0.15, 0.3, 0.5);

% Visualizzazione della slice modificata
figure;
imshow(stretchedSlice(:,:,slice_idx));
se = strel("arbitrary", 2);

figure;
imshow(imfill(imopen(imfilter(stretchedSlice(:,:,slice_idx), fspecial("gaussian", 10)), se) > nBins / 3, "holes"))
figure;
imshow(imclose(stretchedSlice(:,:,slice_idx), se))

if visualizeSlicesFlag
    visualizeSlices(normalizedSlice, stretchedSlice);
end


disp("Code End")