function [normalizedVolume, roiMeanValue] = normalizeVolumeByRoi(inputVolume, roiParams, outputMaxValue)
% normalizeVolumeByRoi Normalizes a 3D volume based on ROI statistics.
%
% The process involves:
%   1. Per-slice initial intensity adjustment (shift and stretch).
%   2. Calculation of mean and max intensity within a circular ROI for each adjusted slice.
%   3. Calculation of global mean and max ROI values across all slices.
%   4. Per-slice gamma correction to match slice ROI statistics to global ROI statistics.
%
% Args:
%   inputVolume (uint16 array): The input 3D volume.
%   roiParams (struct): Parameters for the ROI.
%                       Fields: x (center x), y (center y), r (radius).
%   outputMaxValue (double): The maximum value for the uint16 output (e.g., 65535).
%
% Returns:
%   normalizedVolume (uint16 array): The 3D normalized volume.
%   roiMeanValue (double): The average mean intensity value within the ROI
%                          across all slices *before* the final gamma correction.

    [nRows, nCols, nSlices] = size(inputVolume);

    % Preallocate arrays
    processedSlices = zeros(nRows, nCols, nSlices, 'uint16');
    partialRoiMeans = zeros(nSlices, 1);
    partialRoiMaxs = zeros(nSlices, 1);

    % --- Pass 1: Initial per-slice processing and ROI statistics ---
    for i = 1:nSlices
        currentSlice = inputVolume(:, :, i);

        % Shift and normalize intensity (mimicking original logic)
        sliceShifted = currentSlice + 1024; % Original shift value
        sliceNormDouble = double(sliceShifted) ./ double(outputMaxValue);
        
        % Stretch to [0, 1]
        minVal = min(sliceNormDouble(:));
        maxVal = max(sliceNormDouble(:));
        if minVal == maxVal % Avoid division by zero if slice is flat
            sliceStretched = zeros(size(sliceNormDouble));
        else
            sliceStretched = imadjust(sliceNormDouble, [minVal maxVal], [0 1]);
        end
        processedSlices(:, :, i) = uint16(sliceStretched .* double(outputMaxValue));

        % Compute mean and max in ROI for the processed slice
        [X, Y] = meshgrid(1:nCols, 1:nRows);
        mask = (X - roiParams.x).^2 + (Y - roiParams.y).^2 <= roiParams.r.^2;
        roiPixels = processedSlices(mask & (processedSlices(:,:,i) > 0)); % Consider only non-zero pixels in ROI
        if ~isempty(roiPixels)
            partialRoiMeans(i) = mean(roiPixels);
            partialRoiMaxs(i) = max(roiPixels);
        else
            partialRoiMeans(i) = 0; % Or some other default / NaN
            partialRoiMaxs(i) = 0;  % Or some other default / NaN
        end
    end

    % Calculate global mean of ROI means and global mean of ROI maxs
    roiMeanValue = mean(partialRoiMeans(partialRoiMeans > 0)); % Average of valid means
    meanOfRoiMaxs = mean(partialRoiMaxs(partialRoiMaxs > 0)); % Average of valid maxs
    if isnan(roiMeanValue), roiMeanValue = 0; end
    if isnan(meanOfRoiMaxs), meanOfRoiMaxs = double(outputMaxValue); end % Fallback

    % --- Pass 2: Final normalization using gamma correction ---
    normalizedVolume = zeros(size(processedSlices), 'uint16');
    for i = 1:nSlices
        if partialRoiMaxs(i) == 0 || partialRoiMeans(i) == 0 % Avoid division by zero or log(0)
            normalizedVolume(:, :, i) = processedSlices(:, :, i); % Keep as is or set to 0
            continue;
        end
        ratio = meanOfRoiMaxs / partialRoiMaxs(i);
        adjustedPartialMean = partialRoiMeans(i) * ratio;
        gamma = log(roiMeanValue) / log(adjustedPartialMean); % log_base(X, base) = log(X)/log(base)
        
        normalizedSliceValues = (double(processedSlices(:, :, i)) .* ratio) .^ gamma;
        normalizedVolume(:, :, i) = uint16(min(max(normalizedSliceValues,0), double(outputMaxValue))); % Clip and cast
    end
end