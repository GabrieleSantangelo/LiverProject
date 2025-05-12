function initialLiverMask = segmentLiverByThreshold2D(volume, lowerThreshold, upperThreshold, imfillConn, openStrelDiskSize, verbose, pauseDuration)
% segmentLiverByThreshold2D Segments liver initially using 2D thresholding and morphological operations per slice.
%
% Args:
%   volume (double array): Input volume, assumed to be normalized [0,1].
%   lowerThreshold (double): Lower intensity threshold for liver.
%   upperThreshold (double): Upper intensity threshold for liver.
%   imfillConn (integer): Connectivity for imfill (e.g., 4, 8, or 18 for 3D-like behavior on slices).
%   openStrelDiskSize (integer): Radius of disk for morphological opening.
%   verbose (logical): If true, display intermediate steps.
%   pauseDuration (double): Pause for visualization.
%
% Returns:
%   initialLiverMask (logical array): Initial 3D binary mask of the liver.

    [nRows, nCols, nSlices] = size(volume);
    initialLiverMask = false(nRows, nCols, nSlices);
    se_open = strel('disk', openStrelDiskSize);

    for slice_idx = 1:nSlices
        currentSlice = volume(:,:,slice_idx);
        mask_slice = (currentSlice >= lowerThreshold) & (currentSlice <= upperThreshold);
        
        % Using imfillConn for 2D slices, typically 4 or 8.
        % If 'holes' is meant for 3D context across slices, this needs a 3D approach.
        % Assuming 2D fill for now.
        mask_fill = imfill(mask_slice, imfillConn, "holes"); % imfillConn for 2D is usually 4 or 8
        mask_opened = imopen(mask_fill, se_open);
        initialLiverMask(:,:,slice_idx) = mask_opened;
        
        if verbose && mod(slice_idx, 20) == 0
            figure(200); clf;
            subplot(1,2,1); imshow(currentSlice); title(sprintf('Filtered Slice %d', slice_idx));
            subplot(1,2,2); imshow(mask_opened); title('Initial Liver Mask');
            drawnow; pause(pauseDuration);
        end
    end
end