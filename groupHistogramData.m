function [grouped_hMean, grouped_hMean_clean, binCenters] = groupHistogramData(hMean, hMean_clean, groupSize, nBins)
    nGroups = ceil(nBins / groupSize);
    grouped_hMean = zeros(nGroups,1);
    grouped_hMean_clean = zeros(nGroups,1);
    binCenters = zeros(nGroups,1);
    
    for i = 1:nGroups
        idx_start = (i-1)*groupSize + 1;
        idx_end = min(i*groupSize, nBins);
        grouped_hMean(i) = sum(hMean(idx_start:idx_end));
        grouped_hMean_clean(i) = sum(hMean_clean(idx_start:idx_end));
        binCenters(i) = (idx_start + idx_end - 1) / 2;
    end
end