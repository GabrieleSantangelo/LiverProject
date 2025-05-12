function [lowerThreshold, upperThreshold] = getThresholdsFromRoiInteractive(slice, titleMessage)
% getThresholdsFromRoiInteractive Allows user to draw a polygon ROI on a slice
% and suggests lower/upper thresholds based on min/max intensity in the ROI.
%
% Args:
%   slice (numeric array): 2D image slice.
%   titleMessage (char): Message to display as title for the ROI drawing figure.
%
% Returns:
%   lowerThreshold (double): Minimum intensity value in the ROI (normalized if slice is integer).
%   upperThreshold (double): Maximum intensity value in the ROI (normalized if slice is integer).

    if nargin < 2 || isempty(titleMessage)
        titleMessage = 'Draw ROI to determine thresholds';
    end

    hFig = figure;
    imshow(slice, []);
    title(titleMessage);
    
    roi = drawpolygon; % Or drawfreehand, drawellipse
    wait(roi); % Wait for user to finish drawing

    mask = createMask(roi);
    roiPixels = slice(mask);
    
    close(hFig); % Close the drawing figure

    if isempty(roiPixels)
        error('No pixels selected in ROI. Cannot determine thresholds.');
    end

    lowerThreshold = min(roiPixels(:));
    upperThreshold = max(roiPixels(:));
    
    % If slice is integer (e.g. uint16), normalize thresholds to [0,1]
    if isinteger(slice)
        maxVal = double(intmax(class(slice)));
        lowerThreshold = double(lowerThreshold) / maxVal;
        upperThreshold = double(upperThreshold) / maxVal;
    end
    
    fprintf('Interactive ROI selection: Min=%.4f, Max=%.4f (normalized if applicable)\n', lowerThreshold, upperThreshold);
end