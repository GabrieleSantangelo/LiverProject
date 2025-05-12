function [lowerIntensity, upperIntensity] = bandDetection(histogram, k)
    meanFreq = mean(histogram);

    lowerIntensity=0;
    upperIntensity=1;

    dim = length(histogram);

    slope_arr = (circshift(histogram, -1) - histogram) .* dim;
    slope_arr(dim) = [];
    
    if nargin < 2
        k = median(abs(slope_arr));
    end

    for i=1:dim-1
        freq = histogram(i);

        if freq > meanFreq
            continue
        end

        if (slope_arr(i) > k) && (lowerIntensity == 0)
            lowerIntensity = i / dim;
        end

        if (slope_arr(dim-i) < - k) && (upperIntensity == 1)
            upperIntensity = i / dim;
        end

        if (lowerIntensity ~= 0) && (upperIntensity ~= 1)
            break
        end


    end
end

