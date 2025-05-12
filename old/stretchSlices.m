function [stretchedSlice] = stretchSlices(normalizedSlice, lowerIntensity, upperIntensity, gamma)
    if nargin < 4
        gamma = 1;
    end
        
    stretchedSlice = zeros(size(normalizedSlice), 'uint16');
    dims = size(normalizedSlice);
    nSlice = dims(3);

    for i=1:nSlice
        stretchedSlice(:,:,i) = imadjust(normalizedSlice(:,:,i), [lowerIntensity upperIntensity], [0 1], gamma);
    end
end

