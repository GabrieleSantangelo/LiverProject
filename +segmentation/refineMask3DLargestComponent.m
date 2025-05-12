function refinedMask3D = refineMask3DLargestComponent(initialMask, sphereStrelRadius, minVolumeVoxels, connectivity, verbose, pauseDuration, displayVolume)
% refineMask3DLargestComponent Refines a 3D binary mask using 3D morphological operations
% and keeps only the largest connected component.
%
% Args:
%   initialMask (logical array): The initial 3D binary mask.
%   sphereStrelRadius (integer): Radius of the spherical structuring element for 3D opening.
%   minVolumeVoxels (integer): Minimum volume in voxels to keep a 3D connected component.
%   connectivity (integer): Connectivity for 3D operations (e.g., 6, 18, 26).
%   verbose (logical): If true, display intermediate steps.
%   pauseDuration (double): Pause for visualization.
%   displayVolume (numeric array, optional): Volume to display alongside mask for context.
%
% Returns:
%   refinedMask3D (logical array): The 3D refined binary mask.

    filledVolume_3D = imfill(initialMask, 'holes'); % Works in 3D
    
    se_3D = strel('sphere', sphereStrelRadius);
    openedVolume_3D = imopen(filledVolume_3D, se_3D);
    
    cleanedVolume_3D = bwareaopen(openedVolume_3D, minVolumeVoxels, connectivity);
    
    CC_3D = bwconncomp(cleanedVolume_3D, connectivity);
    refinedMask3D = false(size(initialMask));

    if CC_3D.NumObjects > 0
        volumes = cellfun(@numel, CC_3D.PixelIdxList);
        [~, idxLargest] = max(volumes);
        refinedMask3D(CC_3D.PixelIdxList{idxLargest}) = true;
    end

    if verbose
        slice_to_show = round(size(initialMask,3)/2);
        figure(202); clf;
        if nargin > 6 && ~isempty(displayVolume)
            subplot(1,2,1); imshow(displayVolume(:,:,slice_to_show),[]); title(sprintf('Base Im Slice %d', slice_to_show));
            hold on; visboundaries(initialMask(:,:,slice_to_show),'Color','y'); hold off;
            subplot(1,2,2); imshow(displayVolume(:,:,slice_to_show),[]); title('3D Refined Mask');
            hold on; visboundaries(refinedMask3D(:,:,slice_to_show),'Color','g'); hold off;
        else
            % Fallback if displayVolume is not provided
            subplot(1,2,1); imshow(initialMask(:,:,slice_to_show),[]); title(sprintf('Initial Mask Slice %d', slice_to_show));
            subplot(1,2,2); imshow(refinedMask3D(:,:,slice_to_show),[]); title('3D Refined Mask');
        end
        drawnow; pause(pauseDuration);
    end
end