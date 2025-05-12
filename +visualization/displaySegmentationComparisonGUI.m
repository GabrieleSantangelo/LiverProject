function displaySegmentationComparisonGUI(baseImageVolume, groundTruthVolume, predictionVolume)
% displaySegmentationComparisonGUI Creates a simple GUI to scroll through and
% compare segmentation results against ground truth.
%
% Args:
%   baseImageVolume (numeric array): The grayscale image volume for background.
%   groundTruthVolume (numeric array): Ground truth labels (0, 1:liver, 2:tumor).
%   predictionVolume (numeric array): Predicted labels (0, 1:liver, 2:tumor).

    [nRows, nCols, nSlices] = size(baseImageVolume);

    % Create figure and UI controls
    hFig = figure('Name', 'Segmentation Comparison', 'NumberTitle', 'off', ...
                  'Position', [100, 100, 900, 450]);

    hAx1 = subplot(1, 2, 1, 'Parent', hFig);
    hAx2 = subplot(1, 2, 2, 'Parent', hFig);

    % Slider for slice navigation
    uicontrol('Style', 'slider', ...
              'Min', 1, 'Max', nSlices, 'Value', 1, ...
              'Position', [150, 20, 600, 20], ...
              'Callback', @slider_callback, ...
              'SliderStep', [1/(nSlices-1), 10/(nSlices-1)]); % Minor, Major step

    hSliceText = uicontrol('Style', 'text', ...
                           'Position', [400, 0, 100, 20], ...
                           'String', 'Slice: 1');

    % Initial display
    currentSlice = 1;
    updateDisplay(currentSlice);

    function slider_callback(source, ~)
        currentSlice = round(source.Value);
        source.Value = currentSlice; % Snap to integer
        set(hSliceText, 'String', ['Slice: ', num2str(currentSlice)]);
        updateDisplay(currentSlice);
    end

    function updateDisplay(sliceIdx)
        % Ground Truth Display
        axes(hAx1); %#ok<LAXES>
        imshow(baseImageVolume(:,:,sliceIdx), []);
        hold on;
        visboundaries(groundTruthVolume(:,:,sliceIdx) == 1, 'Color', 'g', 'LineWidth', 0.5, 'LineStyle', ':'); % GT Liver
        visboundaries(groundTruthVolume(:,:,sliceIdx) == 2, 'Color', 'b', 'LineWidth', 1);      % GT Tumor
        hold off;
        title(hAx1, sprintf('Ground Truth - Slice %d', sliceIdx));
        axis(hAx1, 'image', 'off');

        % Prediction Display
        axes(hAx2); %#ok<LAXES>
        imshow(baseImageVolume(:,:,sliceIdx), []);
        hold on;
        visboundaries(predictionVolume(:,:,sliceIdx) == 1, 'Color', 'y', 'LineWidth', 0.5, 'LineStyle', ':'); % Pred Liver
        visboundaries(predictionVolume(:,:,sliceIdx) == 2, 'Color', 'r', 'LineWidth', 1);      % Pred Tumor
        hold off;
        title(hAx2, sprintf('Prediction - Slice %d', sliceIdx));
        axis(hAx2, 'image', 'off');
        
        drawnow;
    end

    % Ensure GUI remains responsive
    if ~isdeployed
        waitfor(hFig); % Keep figure open until closed by user in MATLAB environment
    end

end