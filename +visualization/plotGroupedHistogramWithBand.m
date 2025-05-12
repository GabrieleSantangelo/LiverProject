function plotGroupedHistogramWithBand(binCenters, groupedHist1, groupedHist2, nOriginalBins, lowerIntensityNorm, upperIntensityNorm, titleStr)
% plotGroupedHistogramWithBand Plots grouped histograms and detected intensity band.
%
% Args:
%   binCenters (double array): Center positions of the grouped bins.
%   groupedHist1 (double array): First grouped histogram (e.g., original).
%   groupedHist2 (double array): Second grouped histogram (e.g., cleaned).
%   nOriginalBins (integer): Total number of bins in the original (ungrouped) histogram.
%                            Used to scale normalized intensity bands.
%   lowerIntensityNorm (double): Normalized lower intensity bound (0 to 1).
%   upperIntensityNorm (double): Normalized upper intensity bound (0 to 1).
%   titleStr (char, optional): Overall title for the figure.

    if nargin < 7 || isempty(titleStr)
        titleStr = 'Grouped Histograms with Detected Band';
    end

    binCenters_vec = binCenters(:); % Ensure binCenters is a column vector
    groupedHist1_vec = groupedHist1(:); % Ensure groupedHist1 is a column vector
    groupedHist2_vec = groupedHist2(:); % Ensure groupedHist2 is a column vector

    % Handle potential non-unique binCenters for bar plotting by aggregating counts
    [uniqueBinCenters, ~, groupIdx] = unique(binCenters_vec);
    
    aggregatedHist1 = groupedHist1_vec;
    aggregatedHist2 = groupedHist2_vec;

    if length(uniqueBinCenters) < length(binCenters_vec)
        warningMessage = sprintf('%s: Duplicate binCenters detected. Aggregating histogram counts.', titleStr);
        warning(warningMessage);
        
        % Aggregate groupedHist1
        if ~isempty(groupedHist1_vec)
            aggregatedHist1 = accumarray(groupIdx, groupedHist1_vec, [length(uniqueBinCenters) 1], @sum, 0);
        end

        % Aggregate groupedHist2
        if ~isempty(groupedHist2_vec)
            aggregatedHist2 = accumarray(groupIdx, groupedHist2_vec, [length(uniqueBinCenters) 1], @sum, 0);
        end
    end

    figure('Name', titleStr, 'NumberTitle', 'off');
    subplot(2,1,1);
    if ~isempty(aggregatedHist1)
        bar(uniqueBinCenters, aggregatedHist1);
    end
    title('Grouped Histogram (Original)');
    xlabel('Bin Center (Original Scale)');
    ylabel('Frequency');
    
    subplot(2,1,2);
    if ~isempty(aggregatedHist2)
        bar(uniqueBinCenters, aggregatedHist2);
    end
    title('Grouped Histogram (Cleaned/Processed) with Detected Band');
    xlabel('Bin Center (Original Scale)');
    ylabel('Frequency');
    hold on;
    xline(lowerIntensityNorm * nOriginalBins, 'r--', 'LineWidth', 2, 'Label', 'Lower Bound');
    xline(upperIntensityNorm * nOriginalBins, 'r--', 'LineWidth', 2, 'Label', 'Upper Bound');
    hold off;
end