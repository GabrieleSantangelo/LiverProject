function [grouped_hMean, grouped_hMean_clean, binCenters] = groupHistogram(hMean, hMean_clean, groupSize, nOriginalBins)
% groupHistogram Groups histogram bins.
%
% Args:
%   hMean (double array): The original mean histogram.
%   hMean_clean (double array): The original mean histogram (cleaned version).
%   groupSize (integer): The number of original bins to combine into one new bin.
%   nOriginalBins (integer): The total number of bins in the original histograms.
%
% Returns:
%   grouped_hMean (double array): The grouped version of hMean.
%   grouped_hMean_clean (double array): The grouped version of hMean_clean.
%   binCenters (double array): The center positions of the new grouped bins.

    nGroups = ceil(nOriginalBins / groupSize);
    grouped_hMean = zeros(nGroups, 1);
    grouped_hMean_clean = zeros(nGroups, 1);
    binCenters = zeros(nGroups, 1);

    for i = 1:nGroups
        idx_start = (i - 1) * groupSize + 1;
        idx_end = min(i * groupSize, nOriginalBins);

        grouped_hMean(i) = sum(hMean(idx_start:idx_end));
        grouped_hMean_clean(i) = sum(hMean_clean(idx_start:idx_end));
        % Calculate bin center based on original bin indices (1-based)
        binCenters(i) = (idx_start + idx_end) / 2; % Center of the group range
    end
end