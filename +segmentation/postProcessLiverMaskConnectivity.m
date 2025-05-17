function finalProcessedMask = postProcessLiverMaskConnectivity(inputMask3D, nRowsFullVolume, nRowsFractionForBottom, connectivity3D, verbose, pauseDuration, displayVolume)
% postProcessLiverMaskConnectivity Post-processes a 3D liver mask by:
%   1. Per slice: Keeping only the component that extends furthest down,
%      if any component is in the 'bottom half' of the slice.
%   2. Globally: Keeping only the largest 3D connected component from the result of step 1.
%
% Args:
%   inputMask3D (logical array): The 3D binary liver mask to process.
%   nRowsFullVolume (integer): Total number of rows in the original volume.
%   nRowsFractionForBottom (double): Factor to define the 'bottom half' (e.g., 1.5 means lower 1/1.5 part).
%   connectivity3D (integer): Connectivity for 3D bwconncomp (e.g., 26).
%   verbose (logical): If true, display intermediate steps.
%   pauseDuration (double): Pause for visualization.
%   displayVolume (numeric array, optional): Volume to display alongside mask for context.
%
% Returns:
%   finalProcessedMask (logical array): The post-processed 3D binary liver mask.

    [~, nCols, nSlices] = size(inputMask3D);
    maskAfterStage1 = false(size(inputMask3D));
    startRowBottomHalf = floor(nRowsFullVolume / nRowsFractionForBottom) + 1;

    % Stage 1: Per-slice processing
    for slice_idx = 1:nSlices
        currentSliceMask = inputMask3D(:,:,slice_idx);
        if ~any(currentSliceMask(:)), continue; end

        bottomHalfMaskPortion = currentSliceMask(startRowBottomHalf:end, :);
        if ~any(bottomHalfMaskPortion(:)), continue; end
        
        CC_2D = bwconncomp(currentSliceMask, 8);
        if CC_2D.NumObjects <= 1
            maskAfterStage1(:,:,slice_idx) = currentSliceMask;
        else
            maxRowOverall = -1; idxLowestComponent = -1;
            for k = 1:CC_2D.NumObjects
                [rows, ~] = ind2sub([nRowsFullVolume, nCols], CC_2D.PixelIdxList{k});
                if isempty(rows), continue; end
                if max(rows) > maxRowOverall
                    maxRowOverall = max(rows);
                    idxLowestComponent = k;
                end
            end
            if idxLowestComponent > 0
                tempSliceMask = false(nRowsFullVolume, nCols);
                tempSliceMask(CC_2D.PixelIdxList{idxLowestComponent}) = true;
                maskAfterStage1(:,:,slice_idx) = tempSliceMask;
            end
        end
    end

    % Stage 2: Select largest 3D component from Stage 1 result
    finalProcessedMask = false(size(inputMask3D));
    if any(maskAfterStage1(:))
        CC_3D_Stage2 = bwconncomp(maskAfterStage1, connectivity3D);
        if CC_3D_Stage2.NumObjects > 0
            numVoxels = cellfun(@numel, CC_3D_Stage2.PixelIdxList);
            [~, idxLargestComponent] = max(numVoxels);
            finalProcessedMask(CC_3D_Stage2.PixelIdxList{idxLargestComponent}) = true;
        end
    end

    if verbose
        for slice_to_show=1:size(inputMask3D, 3)
            figure(203); clf;
            subplot(1,2,1); imshow(displayVolume(:,:,slice_to_show),[]); title(sprintf('Base Im Slice %d', slice_to_show));
            hold on; visboundaries(inputMask3D(:,:,slice_to_show),'Color','y'); hold off;
            subplot(1,2,2); imshow(displayVolume(:,:,slice_to_show),[]); title('Post-processed Liver Mask (Connectivity)');
            hold on; visboundaries(finalProcessedMask(:,:,slice_to_show),'Color','m'); hold off;
            drawnow; pause(pauseDuration);
        end
    end
end