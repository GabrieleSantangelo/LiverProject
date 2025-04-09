clear; clc; clf;
disp("Code Run")

% Option to visualize slices
visualizeSlicesFlag = true;


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
%[text] Define the region of interest (ROI) parameters for cross normalization
roiParams.x = 35;
roiParams.y = 38;
roiParams.r = 30;


% Maximum value for uint16
maxValue = 65536;

% ----------- START REFACTOR FOLLOWING CODE ----------- %

% Number of bins for histogram
nBins = 65536; % Number of bins for histogram
[meanValue, normalizedSlice_temp] = normalizingSlices(trainVolume, roiParams, maxValue);

normalizedSlice = histogramMachingAllSlice(normalizedSlice_temp, 150);
[hMean, hMean_clean] = histogramOnAllSlices(normalizedSlice, nBins);

visualizeSlices(normalizedSlice, trainVolume)
%%
fprintf('Mean value in ROI: %.2f\n', meanValue);

% Visualize mean histogram of all slices
plotHistograms(hMean, hMean_clean);

% Group histogram data
groupSize = 2000;

[grouped_hMean, grouped_hMean_clean, binCenters] = groupHistogramData(hMean, hMean_clean, groupSize, nBins);
[lowerIntensity, upperIntensity] = bandDetection(grouped_hMean, 3000000); % 30000 con grouped_hMean_clean | 3000000 con grouped_hMean

plotGroupedHistograms(binCenters, grouped_hMean, grouped_hMean_clean, nBins, lowerIntensity, upperIntensity);
%%
% ------------ STOP REFACTOR FOLLOWING CODE ----------- %
% Stretch the slices using the lower and upper intensity values
stretchedSlice = stretchSlices(normalizedSlice, lowerIntensity, upperIntensity, 8);

dims = size(trainVolume);
nSlice = dims(3);

no_bones_slice = zeros(size(trainVolume,1), size(trainVolume,2), nSlice,'uint16');


   
figure;
imshow(stretchedSlice(:,:,195))
title("Porco caneeeeee");
for slice_idx=1:nSlice


    % Remove high value (in this case mainly areas representing the bones)
    mask = stretchedSlice(:,:,slice_idx) > maxValue * 0.9;
    mask = imdilate(double(mask), strel("disk", 4)); % Dilation to remove small holes
    mask = imerode(mask, strel("disk", 2)); % Erosion to remove small holes
    mask = imfill(mask, "holes"); % Fill holes
    mask = imfilter(mask, fspecial("gaussian", 10)); % Gaussian filter to smooth the mask

    tempSlice = uint16(double(normalizedSlice(:,:,slice_idx)) .* double(1 - mask));
    no_bones_slice(:,:,slice_idx) = tempSlice;
    figure(1);clf;
    imshow(no_bones_slice(:,:,slice_idx))
    pause(0.001)
end


% Visualizzazione della slice modificata
figure;
imshow(no_bones_slice(:,:,slice_idx));
%%

[hMean_no_bones, hMean_clean_no_bones] = histogramOnAllSlices(no_bones_slice, maxValue);

[grouped_hMean_no_bones, grouped_hMean_clean_no_bones, binCenters] = groupHistogramData(hMean_no_bones, hMean_clean_no_bones, groupSize, nBins);


[lowerIntensity_no_bones, upperIntensity_no_bones] = bandDetection(grouped_hMean_no_bones);

plotGroupedHistograms(binCenters, grouped_hMean_no_bones, grouped_hMean_clean_no_bones, nBins, lowerIntensity_no_bones, upperIntensity_no_bones);


%%
flag = false

doubleStretchedSlice = stretchSlices(no_bones_slice, lowerIntensity_no_bones, upperIntensity_no_bones, 5.8);
doubleStretchedSlice(1:150,: ,:) = 0;

if flag
    for slice_idx=1:nSlice
    
        figure(1); clf;
        subplot(1, 2, 1);
        imshow(doubleStretchedSlice(:,:,slice_idx))
        title(['Label Slice', num2str(slice_idx)]);
    
        subplot(1, 2, 2);
        imshow(no_bones_slice(:,:,slice_idx))
        title(['Label Slice', num2str(slice_idx)]);
        pause(0.0001);
    end
end

sliceDouble = double(doubleStretchedSlice);
mu = mean(sliceDouble(:));
sigma = std(sliceDouble(:));

lower = mu + 0.7*sigma;
upper = mu + 3*sigma;

% temp_1 = doubleStretchedSlice > maxValue / 14 & doubleStretchedSlice < maxValue / 2; %14 e 2
% temp_1 = sliceDouble > lower & sliceDouble < upper;
for slice_idx=1:nSlice

    level = graythresh(doubleStretchedSlice(:,:,slice_idx));
    temp_1(:,:,slice_idx) = imbinarize(doubleStretchedSlice(:,:,slice_idx), level);

    figure(1); clf;
    imshow(temp_1(:,:,slice_idx))
    title("temp_1")
    pause(0.0001)
end


%%


temp_2 = imfill(temp_1, 18, "holes"); % 18
for slice_idx=1:nSlice
    figure(1); clf;
    imshow(temp_2(:,:,slice_idx))
    title("temp_2")
    pause(0.0001)
end


temp_3 = imerode(BW_final, strel('diamond', 10));
for slice_idx=1:nSlice
    figure(1); clf;
    imshow(temp_3(:,:,slice_idx))
    title("temp_3")
    pause(0.0001)
end

mask = imfill(temp_3, 26, "holes");


for slice_idx=1:nSlice
    figure(1); clf;
    imshow(mask(:,:,slice_idx))
    title("mask")
    pause(0.0001)
end



for slice_idx=1:nSlice
    
    figure(1); clf;
    subplot(1, 2, 1);
    imshow(doubleStretchedSlice(:,:,slice_idx))
    title(['Label Slice', num2str(slice_idx)]);
    
    se = strel('disk', 7); % Prova con diversi raggi
    mask = imopen(mask, se);
    subplot(1, 2, 2);
    imshow(mask(:,:,slice_idx),[])
    title(['Label Slice', num2str(slice_idx)]);
    pause(0.0001);
end
%%
% ------------ START OF TEST CODE ----------- %
%slice_idx = 154
mod = zeros(size(trainVolume,1), size(trainVolume,2), nSlice, 'uint16');

for slice_idx=1:nSlice

    figure(2); clf;
    tempSlice = doubleStretchedSlice(:,:,slice_idx);
    
    [~,threshold] = edge(tempSlice,'sobel');
    fudgeFactor = 0.5;
    BWs = edge(tempSlice,'sobel',threshold * fudgeFactor);

    subplot(1, 2, 1);
    imshow(tempSlice);
    
    se90 = strel('line',3,90);
    se0 = strel('line',3,0);
    
    BWsdil = imdilate(BWs,[se90 se0]);
    %imshow(BWsdil)
    %title('Dilated Gradient Mask')
    
    BWdfill = imfill(BWsdil,'holes');
    %imshow(BWdfill)
    %title('Binary Image with Filled Holes')
    
    
    BWnobord = imclearborder(BWdfill,4);
    %imshow(BWnobord)
    %title('Cleared Border Image')
    
    seD = strel('diamond',1);
    BWfinal = imerode(BWnobord,seD);
    BWfinal = imerode(BWfinal,seD);
    %imshow(BWfinal)
    %title('Segmented Image');
    
    % seD5 = strel('diamond',8);
    % BWfinal = imopen(BWfinal, seD5);
    %imshow(BWfinal)
    %title("AntoMetod")

    subplot(1, 2, 2);
    imshow(BWfinal)
    title(['Label Slice', num2str(slice_idx)]);

    pause(0.01);
    %mod(:,:,slice_idx) = BWfinal;

end

%visualizeSlices(mod, labelVolume)
%%

figure;
% Visualizzazione della slice modificata
imshow(doubleStretchedSlice(:,:,slice_idx));
se_open = strel("disk", 8);

figure;
imshow(imopen(imfill(imfilter(doubleStretchedSlice(:,:,slice_idx), fspecial("disk", 3)) > nBins / 6, "holes"), se_open))

if visualizeSlicesFlag
    visualizeSlices(normalizedSlice, stretchedSlice);
end


disp("Code End")

%[appendix]
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":39.7}
%---
