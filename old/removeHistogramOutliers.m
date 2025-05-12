function hMean_clean = removeHistogramOutliers(hMean)
    outliers = isoutlier(hMean, "percentile",[2 98]);
    hMean_clean = hMean;
    
    x = (1:length(hMean))';
    x_good = x(~outliers);
    y_good = hMean(~outliers);
    x_bad = x(outliers);
    
    hMean_clean(outliers) = interp1(x_good, y_good, x_bad, 'linear', 'extrap');
end
