function [matchedVolume] = histogramMachingAllSlice(slices, reference_idx)
    matchedVolume = slices;  
    

    if nargin < 2
        for i = 2:size(slices, 3)
            currentSlice = slices(:, :, i-1);
            referenceSlice = matchedVolume(:, :, i);  
            nBins = length(unique(referenceSlice(:)));
            disp(nBins)
            if nBins < 2
                nBins = 2;
            end
            matchedVolume(:, :, i) = imhistmatch(currentSlice, referenceSlice, nBins);
        end
    else
        referenceSlice = matchedVolume(:, :, reference_idx); 
        nBins = length(unique(referenceSlice(:)));
        for i = 1:size(slices, 3)
            currentSlice = slices(:, :, i);
            matchedVolume(:, :, i) = imhistmatch(currentSlice, referenceSlice, 128);
        end
    end

end