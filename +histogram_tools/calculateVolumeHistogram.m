function [hMean, hMeanClean] = calculateVolumeHistogram(volume, nBins, ignoreZeros)
% calculateVolumeHistogram Calculates the mean histogram over all slices in a volume.
%
% Args:
%   volume (numeric array): The input 3D volume.
%   nBins (integer): The number of bins for the histogram.
%   ignoreZeros (logical, optional): If true, zeros in the volume are
%                                    ignored during histogram calculation.
%                                    Default is false.
%
% Returns:
%   hMean (double array): The mean histogram across all slices.
%   hMeanClean (double array): The mean histogram with outliers removed.

    if nargin < 3
        ignoreZeros = false;
    end

    nSlices = size(volume, 3);
    hTotal = zeros(nBins, 1);

    for i = 1:nSlices
        currentSlice = volume(:, :, i);
        if ignoreZeros
            pixels = currentSlice(currentSlice > 0);
            if isempty(pixels) % Handle case where slice is all zeros or empty after filtering
                h = zeros(nBins,1);
            else
                h = imhist(pixels, nBins);
            end
        else
            h = imhist(currentSlice, nBins);
        end
        hTotal = hTotal + h;
    end

    hMean = hTotal / nSlices;
    hMeanClean = histogram_tools.removeHistogramOutliers(hMean); % Assumes removeHistogramOutliers is in +histogram_tools
end