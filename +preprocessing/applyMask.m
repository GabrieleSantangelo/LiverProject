function maskedVolume = applyMask(inputVolume, mask)
% applyMask Applies a binary mask to an input volume.
%   Sets voxels to zero where the mask is false.
%
% Args:
%   inputVolume (numeric array): The volume to be masked.
%   mask (logical array): The binary mask (same size as inputVolume).
%                         True values indicate regions to keep.
%
% Returns:
%   maskedVolume (numeric array): The masked volume, same type as input.

    % Ensure mask is logical for element-wise multiplication behavior
    maskedVolume = inputVolume .* cast(mask, class(inputVolume));
end