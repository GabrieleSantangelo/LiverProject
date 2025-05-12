function [meanValue, hMeanNorm, hMean_clean] = analyzeSlices(normalizedSlice, allSlices, partialMeanValues, nBins)
    nSlice = size(normalizedSlice, 3);

    meanValue = mean(partialMeanValues);

    % Istogramma medio non normalizzato
    hTotal = zeros(nBins,1);
    for i = 1:nSlice
        h = imhist(allSlices(:, :, i), nBins);
        hTotal = hTotal + h;
    end
    hMean = hTotal / nSlice;
    hMean_clean = removeHistogramOutliers(hMean);

    % Istogramma medio normalizzato
    hTotalNorm = zeros(nBins,1);
    for i = 1:nSlice
        h_norm = imhist(normalizedSlice(:, :, i), nBins);
        hTotalNorm = hTotalNorm + h_norm;
    end
    hMeanNorm = hTotalNorm / nSlice;
end
