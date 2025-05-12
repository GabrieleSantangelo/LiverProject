function refinedMask2D = refineLiverMask2D(initialMask, minArea, seCloseDisk, seOpenDisk, verbose, pauseDuration)
% refineLiverMask2D Refines a 2D liver mask slice-by-slice using morphological operations.
%
% Args:
%   initialMask (logical array): The initial 3D binary liver mask.
%   minArea (integer): Minimum area in pixels to keep a connected component.
%   seCloseDisk (integer): Radius of disk for morphological closing.
%   seOpenDisk (integer): Radius of disk for morphological opening.
%   verbose (logical): If true, display intermediate steps.
%   pauseDuration (double): Pause for visualization.
%
% Returns:
%   refinedMask2D (logical array): The 2D refined 3D binary mask.

    [nRows, nCols, nSlices] = size(initialMask);
    refinedMask2D = false(nRows, nCols, nSlices);
    
    se_close = strel('disk', seCloseDisk);
    se_open = strel('disk', seOpenDisk);

    for slice_idx = 1:nSlices
        binarySlice = initialMask(:,:,slice_idx);
        filledSlice = imfill(binarySlice, 'holes');
        closedSlice = imclose(filledSlice, se_close);
        openedSlice = imopen(closedSlice, se_open);
        cleanedSlice = bwareaopen(openedSlice, minArea);
        
        CC = bwconncomp(cleanedSlice);
        finalSlice = false(size(binarySlice));
        if CC.NumObjects > 0
            stats = regionprops(CC, 'Area');
            [~, idx] = max([stats.Area]);
            finalSlice(CC.PixelIdxList{idx}) = true;
        end
        refinedMask2D(:,:,slice_idx) = finalSlice;

        if verbose && mod(slice_idx, 20) == 0
            figure(201); clf;
            subplot(1,2,1); imshow(initialMask(:,:,slice_idx)); title(sprintf('Initial Mask Slice %d', slice_idx));
            subplot(1,2,2); imshow(finalSlice); title('2D Refined Mask');
            drawnow; pause(pauseDuration);
        end
    end
end