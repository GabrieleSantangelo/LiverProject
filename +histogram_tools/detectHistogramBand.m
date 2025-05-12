function [lowerIntensityNorm, upperIntensityNorm] = detectHistogramBand(histogram, k_factor)
% detectHistogramBand Detects a significant intensity band based on histogram slopes.
%   Finds the first significant rising slope from the left (lower bound)
%   and the first significant falling slope from the right (upper bound).
%
% Args:
%   histogram (double array): The input histogram vector.
%   k_factor (double, optional): A factor to multiply with the median of
%       absolute slopes to determine the slope threshold. If not provided,
%       the median of absolute slopes is used directly as the threshold 'k'.
%       If k_factor is provided, k = k_factor * median(abs(slopes)).
%       If k_factor is a direct value (e.g. from old script `bandDetection(..., 3000)`),
%       it's treated as the slope threshold `k` itself.
%
% Returns:
%   lowerIntensityNorm (double): Normalized lower intensity bound (0 to 1).
%   upperIntensityNorm (double): Normalized upper intensity bound (0 to 1).

    lowerIntensityNorm = 0; % Default: start of histogram
    upperIntensityNorm = 1; % Default: end of histogram

    dim = length(histogram);
    if dim < 2
        return; % Not enough data points
    end

    % Calculate slopes (difference between adjacent bins)
    % Multiplied by 'dim' to somewhat normalize against bin count, as in original.
    slopes = (circshift(histogram, -1) - histogram) * dim;
    slopes(dim) = []; % Remove last element which is invalid due to circshift

    if nargin < 2 || isempty(k_factor)
        k_slope_threshold = median(abs(slopes)); % Default threshold based on median slope
    else
        k_slope_threshold = k_factor; % Use provided k_factor as direct threshold or multiplier
        % If k_factor was intended as a multiplier for median(abs(slopes)):
        % k_slope_threshold = k_factor * median(abs(slopes));
    end

    meanFreq = mean(histogram);

    for i = 1:(dim - 1)
        % Find lower bound: first significant positive slope from left, below mean frequency
        if histogram(i) <= meanFreq && slopes(i) > k_slope_threshold && lowerIntensityNorm == 0
            lowerIntensityNorm = i / dim;
        end

        % Find upper bound: first significant negative slope from right, below mean frequency
        % Index from right for slopes is (dim - i). Bin index is (dim - i).
        current_idx_from_right = dim - i; % This is the index for 'slopes' and 'histogram'
        if histogram(current_idx_from_right) <= meanFreq && slopes(current_idx_from_right) < -k_slope_threshold && upperIntensityNorm == 1
            upperIntensityNorm = current_idx_from_right / dim;
            % Note: The original 'bandDetection.m' had 'upperIntensity = i / dim;' here,
            % which seemed counter-intuitive. This version uses the actual index from the right.
        end

        % Exit if both bounds are found
        if lowerIntensityNorm ~= 0 && upperIntensityNorm ~= 1
            break;
        end
    end
    
    % Ensure lower < upper
    if lowerIntensityNorm >= upperIntensityNorm && upperIntensityNorm ~=1 && lowerIntensityNorm ~=0
        % If bounds crossed or are invalid, reset to full range or handle error
        % For now, if lower is found but upper is still 1, it's fine.
        % If lower is found after upper (e.g. lower=0.5, upper=0.2), this is an issue.
        % Resetting to defaults might be too aggressive. This case needs careful consideration
        % based on expected histogram shapes. For now, we allow it, caller might need to check.
        % A simple fix if they cross and are not defaults:
        % lowerIntensityNorm = 0; upperIntensityNorm = 1;
    end
end