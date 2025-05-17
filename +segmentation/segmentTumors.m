function finalTumorMask = segmentTumors(liverIntensityVolume, liverMask, lowerThresh, upperThresh, minTumorVol, openSphereRad, closeDiskRad2D, conn3D, maxValue, verbose, pauseDur, labelVolumeForViz)
% segmentTumors Segments tumors within a given liver mask.
%
% Args:
%   liverIntensityVolume (double array): Intensity volume (e.g., bilaterally filtered), normalized [0,1].
%   liverMask (logical array): 3D binary mask of the liver.
%   lowerThresh (double): Lower intensity threshold for tumor candidates within liver.
%   upperThresh (double): Upper intensity threshold for tumor candidates.
%   minTumorVol (integer): Minimum voxel volume for a tumor.
%   openSphereRad (integer): Radius for 3D morphological opening (strel 'sphere').
%   closeDiskRad2D (integer): Radius for 2D morphological closing per slice (strel 'disk').
%   conn3D (integer): Connectivity for 3D bwareaopen.
%   verbose (logical): If true, display intermediate steps.
%   pauseDur (double): Pause for visualization.
%   labelVolumeForViz (numeric array, optional): Ground truth label volume for side-by-side visualization.
%
% Returns:
%   finalTumorMask (logical array): 3D binary mask of detected tumors.

    liverIntensityVolume = double(liverIntensityVolume) / double(maxValue);

    initialTumorMask = false(size(liverIntensityVolume));
    validLiverPixels = liverMask & (liverIntensityVolume > 0); % Consider only non-zero pixels within liver
    
    % Apply thresholding only within the liver mask
    candidatePixels = liverIntensityVolume(validLiverPixels);
    tumorPixelsInLiver = (candidatePixels >= lowerThresh) & (candidatePixels <= upperThresh);
    initialTumorMask(validLiverPixels) = tumorPixelsInLiver;

    % 3D Morphological Opening
    openedMask_3D = imopen(initialTumorMask, strel('sphere', openSphereRad));

    % Remove Small Objects (3D)
    cleanedTumorMask_3D = bwareaopen(openedMask_3D, minTumorVol, conn3D);

    % 2D Morphological Closing (slice-by-slice) and Fill Holes
    closedMask_final = false(size(cleanedTumorMask_3D));
    se_close_2D = strel('disk', closeDiskRad2D);
    for k=1:size(cleanedTumorMask_3D, 3)
         slice_closed = imclose(cleanedTumorMask_3D(:,:,k), se_close_2D);
         closedMask_final(:,:,k) = imfill(slice_closed, "holes");
    end
    finalTumorMask = closedMask_final;
    finalTumorMask = imfill(finalTumorMask, 'holes'); % Final 3D fill

    if verbose
        for slice_to_show=1:size(liverIntensityVolume, 3)
            figure(204); clf;
            subplot(1,2,1); imshow(liverIntensityVolume(:,:,slice_to_show),[]); title(sprintf('Input Intensity Slice %d', slice_to_show));
            hold on; visboundaries(liverMask(:,:,slice_to_show),'Color','g', 'LineStyle', ':'); hold off;
            subplot(1,2,2); imshow(liverIntensityVolume(:,:,slice_to_show),[]); title('Detected Tumors (Red)');
            hold on; visboundaries(liverMask(:,:,slice_to_show),'Color','g', 'LineStyle', ':');
            visboundaries(finalTumorMask(:,:,slice_to_show),'Color','r'); hold off;
            if nargin > 10 && ~isempty(labelVolumeForViz)
                subplot(1,2,1); hold on; visboundaries(labelVolumeForViz(:,:,slice_to_show)==2, 'Color', 'b', 'LineStyle', '--'); hold off; title('Input + GT Tumor (Blue)');
            end
            drawnow; pause(pauseDur);
        end
    end
end