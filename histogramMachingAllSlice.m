function [matchedVolume] = histogramMachingAllSlice(slices, reference_idx)
    matchedVolume = slices;  
    dim = size(slices, 3);
    

    if nargin < 2
        for i = 2:dim
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
        h = waitbar(0);
        for i = 1:dim
            waitbar(i/dim, h);
            currentSlice = slices(:, :, i);
            matchedVolume(:, :, i) = imhistmatch(currentSlice, referenceSlice, nBins);
        end
    end    

end