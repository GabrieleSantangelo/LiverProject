function cleanedHistogram = removeHistogramOutliers(histogramVector, percentiles)
% removeHistogramOutliers Replaces outliers in a histogram with interpolated values.
%
% Args:
%   histogramVector (double array): The input histogram (1D vector).
%   percentiles (double array, optional): A 2-element vector specifying the
%       lower and upper percentiles to define outliers. Default is [2 98].
%
% Returns:
%   cleanedHistogram (double array): The histogram with outliers interpolated.

    if nargin < 2
        percentiles = [2 98]; % Default percentiles
    end

    cleanedHistogram = histogramVector;
    outlierIndices = isoutlier(histogramVector, "percentiles", percentiles);

    if any(outlierIndices)
        x = (1:length(histogramVector))';
        x_good = x(~outlierIndices);
        y_good = histogramVector(~outlierIndices);
        cleanedHistogram(outlierIndices) = interp1(x_good, y_good, x(outlierIndices), 'linear', 'extrap');
    end
end