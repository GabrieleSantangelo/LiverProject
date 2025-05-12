function plotHistogramComparison(hist1, hist2, titleStr)
% plotHistogramComparison Plots two histograms side-by-side for comparison.
%
% Args:
%   hist1 (double array): The first histogram vector.
%   hist2 (double array): The second histogram vector (e.g., cleaned version).
%   titleStr (char, optional): Overall title for the figure.

    if nargin < 3 || isempty(titleStr)
        titleStr = 'Histogram Comparison';
    end

    figure('Name', titleStr, 'NumberTitle', 'off');
    subplot(2,1,1);
    bar(hist1);
    title('Original Histogram');
    xlabel('Intensity');
    ylabel('Frequency');
    
    subplot(2,1,2);
    bar(hist2);
    title('Processed/Cleaned Histogram');
    xlabel('Intensity');
    ylabel('Frequency');
end