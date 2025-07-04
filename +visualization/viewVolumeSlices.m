function viewVolumeSlices(volume1, volume2, titleStr, pauseDuration, adaptToRange1, adaptToRange2)
% viewVolumeSlices Displays slices from one or two 3D volumes side-by-side.
%
% Args:
%   volume1 (numeric array): The first 3D volume.
%   volume2 (numeric array, optional): The second 3D volume. If empty or not
%                                      provided, only volume1 is shown.
%   titleStr (char, optional): Title for the figure window.
%   pauseDuration (double, optional): Pause duration between slices.
%                                     Default is 0.01.

    if nargin < 2, volume2 = []; end
    if nargin < 3 || isempty(titleStr), titleStr = 'Slice Viewer'; end
    if nargin < 4 || isempty(pauseDuration), pauseDuration = 0.001; end
    if nargin < 5 || isempty(adaptToRange1), adaptToRange1 = false; end
    if nargin < 6 || isempty(adaptToRange2), adaptToRange2 = false; end

    nSlices = size(volume1, 3);
    
    figure('Name', titleStr, 'NumberTitle', 'off');
    for i = 1:nSlices
        clf; % Clear current figure
        if ~isempty(volume2)
            subplot(1, 2, 1);
            if adaptToRange1
                imshow(volume1(:, :, i), []);
            else
                imshow(volume1(:, :, i));
            end
            title(sprintf('Volume 1 - Slice %d/%d', i, nSlices));
            subplot(1, 2, 2);
            if adaptToRange2
                imshow(volume2(:, :, i), []);
            else
                imshow(volume2(:, :, i));
            end
            title(sprintf('Volume 2 - Slice %d/%d', i, nSlices));
        else
            if adaptToRange1
                imshow(volume1(:, :, i), []);
            else
                imshow(volume1(:, :, i));
            end
            title(sprintf('Slice %d/%d', i, nSlices));
        end
        drawnow;
        pause(pauseDuration);
    end
end