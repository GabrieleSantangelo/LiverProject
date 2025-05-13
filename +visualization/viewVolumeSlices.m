function viewVolumeSlices(volume1, volume2, titleStr, pauseDuration)
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

    nSlices = size(volume1, 3);
    
    figure('Name', titleStr, 'NumberTitle', 'off');
    for i = 1:nSlices
        clf; % Clear current figure
        if ~isempty(volume2)
            subplot(1, 2, 1);
            imshow(volume1(:, :, i));
            title(sprintf('Volume 1 - Slice %d/%d', i, nSlices));
            subplot(1, 2, 2);
            imshow(volume2(:, :, i));
            title(sprintf('Volume 2 - Slice %d/%d', i, nSlices));
        else
            imshow(volume1(:, :, i));
            title(sprintf('Slice %d/%d', i, nSlices));
        end
        drawnow;
        pause(pauseDuration);
    end
end