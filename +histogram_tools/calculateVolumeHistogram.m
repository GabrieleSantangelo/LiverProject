function [hMean, hMeanClean] = calculateVolumeHistogram(volume, nBins)
% calculateVolumeHistogram Calculates the mean histogram over all slices in a volume.
%
% Args:
%   volume (numeric array): The input 3D volume.
%   nBins (integer): The number of bins for the histogram.
%
% Returns:
%   hMean (double array): The mean histogram across all slices.
%   hMeanClean (double array): The mean histogram with outliers removed.

    nSlices = size(volume, 3);
    hTotal = zeros(nBins, 1);

    for i = 1:nSlices
        currentSlice = volume(:, :, i);

        h = imhist(currentSlice, nBins);
        hTotal = hTotal + h;
    end

    hMean = hTotal / nSlices;
    hMeanClean = histogram_tools.removeHistogramOutliers(hMean); % Assumes removeHistogramOutliers is in +histogram_tools
end