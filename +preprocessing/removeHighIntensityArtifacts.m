function [volumeWithoutArtifacts, artifactMaskCombined] = removeHighIntensityArtifacts(baseVolume, intensityRefVolume, maxValue, params, verbose, pauseDuration)
% removeHighIntensityArtifacts Removes high intensity artifacts (e.g., bones)
% by creating a mask from an intensity reference volume and applying it to a base volume.
%
% Args:
%   baseVolume (uint16 array): The volume to clean (e.g., normalizedSlice).
%   intensityRefVolume (uint16 array): Volume used to identify high intensity
%                                      regions (e.g., stretchedSlice1).
%   maxValue (double): Maximum possible intensity value (e.g., 65535).
%   params (struct): Parameters for artifact removal.
%       .window (int): Number of neighboring slices to consider for mask generation.
%       .thresholdFactor (double): Factor of maxValue to define high intensity.
%       .dilateDiskRadius (int): Radius for dilation.
%       .erodeDiskRadius (int): Radius for erosion.
%       .gaussianSize (int): Size of Gaussian filter for smoothing the mask.
%   verbose (logical): If true, display intermediate steps.
%   pauseDuration (double): Pause for visualization.
%
% Returns:
%   volumeWithoutArtifacts (uint16 array): Base volume with artifacts removed.
%   artifactMaskCombined (logical array): The combined 3D mask of detected artifacts.

    [nRows, nCols, nSlices] = size(baseVolume);
    volumeWithoutArtifacts = zeros(size(baseVolume), 'like', baseVolume);
    artifactMaskCombined = false(nRows, nCols, nSlices);

    for slice_idx = 1:nSlices
        tempSlice = baseVolume(:, :, slice_idx);
        cumulativeMaskForSlice = false(nRows, nCols);

        % Consider a window of slices from intensityRefVolume
        for j = 0:params.window-1
            ref_slice_idx = slice_idx + j;
            if ref_slice_idx > nSlices, continue; end % Boundary check

            mask = intensityRefVolume(:,:,ref_slice_idx) > (maxValue * params.thresholdFactor);
            mask = imdilate(double(mask), strel("disk", params.dilateDiskRadius));
            mask = imerode(mask, strel("disk", params.erodeDiskRadius));
            mask = imfill(mask, "holes");
            mask = imfilter(mask, fspecial("gaussian", params.gaussianSize)) > 0.5; % Threshold smoothed mask
            cumulativeMaskForSlice = cumulativeMaskForSlice | mask;
        end
        volumeWithoutArtifacts(:,:,slice_idx) = uint16(double(tempSlice) .* double(~cumulativeMaskForSlice));
        artifactMaskCombined(:,:,slice_idx) = cumulativeMaskForSlice;
        
        if verbose && mod(slice_idx, 1) == 0 % Display occasionally
            figure(100); clf;
            subplot(1,2,1); imshow(baseVolume(:,:,slice_idx)); title(sprintf('Base Slice %d', slice_idx));
            subplot(1,2,2); imshow(volumeWithoutArtifacts(:,:,slice_idx)); title(sprintf('Cleaned Slice %d', slice_idx));
            sgtitle('Bone Removal Pass 1');
            drawnow; pause(pauseDuration);
        end
    end
end