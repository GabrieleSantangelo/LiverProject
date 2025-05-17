function structureMaskArray = createStructureMask(inputVolume, maxValue, params, verbose, pauseDuration)
% createStructureMask Creates a mask for unwanted structures (e.g., remaining bones/high density).
% The mask identifies regions ABOVE a certain threshold in the inputVolume.
% The returned mask is TRUE for structures to be KEPT (i.e., 1 - original mask).
%
% Args:
%   inputVolume (uint16 array): Volume to process (e.g., doubleStretchedSlice).
%   maxValue (double): Maximum possible intensity value.
%   params (struct): Parameters for mask creation.
%       .thresholdFactor (double): Factor of maxValue to define structure.
%       .dilateDiskRadius (int): Radius for dilation.
%       .erodeDiskRadius (int): Radius for erosion.
%       .gaussianSize (int): Size of Gaussian filter for smoothing.
%   verbose (logical): If true, display intermediate steps.
%   pauseDuration (double): Pause for visualization.
%
% Returns:
%   structureMaskArray (logical array): Mask where TRUE indicates regions to KEEP.

    [nRows, nCols, nSlices] = size(inputVolume);
    structureMaskArray = false(nRows, nCols, nSlices);

    for slice_idx = 1:nSlices
        mask = inputVolume(:,:,slice_idx) > (maxValue * params.thresholdFactor);
        mask = imdilate(double(mask), strel("disk", params.dilateDiskRadius));
        mask = imerode(mask, strel("disk", params.erodeDiskRadius));
        mask = imfill(mask, "holes");
        smoothed_mask = imfilter(mask, fspecial("gaussian", params.gaussianSize));
        
        structureMaskArray(:,:,slice_idx) = ~(smoothed_mask > 0.5); % Invert: true for regions to KEEP

        if verbose && mod(slice_idx, 1) == 0
            figure(101); clf;
            imshow(~structureMaskArray(:,:,slice_idx)); % Show regions being masked out
            title(sprintf('Secondary Structure Mask (masked out) - Slice %d', slice_idx));
            drawnow; pause(pauseDuration);
        end
    end
end