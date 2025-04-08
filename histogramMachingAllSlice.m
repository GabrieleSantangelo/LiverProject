function [matchedVolume] = histogramMachingAllSlice(slices, reference_idx)
    matchedVolume = slices;  
    

    if nargin < 2
        for i = 2:size(slices, 3)
            currentSlice = slices(:, :, i);
            referenceSlice = matchedVolume(:, :, i - 1);  
            matchedVolume(:, :, i) = imhistmatch(currentSlice, referenceSlice);
        end
    else
        referenceSlice = matchedVolume(:, :, reference_idx); 
        for i = 1:size(slices, 3)
            currentSlice = slices(:, :, i);
            matchedVolume(:, :, i) = imhistmatch(currentSlice, referenceSlice);
        end
    end

end