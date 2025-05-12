function stretchedVolume = stretchVolumeContrast(inputVolume, lowerIntensityNorm, upperIntensityNorm, gamma, outputMaxValue, wb_handle)
% stretchVolumeContrast Stretches the contrast of a 3D volume.
%   Each slice is adjusted using imadjust.
%
% Args:
%   inputVolume (numeric array): The input 3D volume (e.g., uint16).
%   lowerIntensityNorm (double): Normalized lower intensity input bound (0 to 1).
%                                This corresponds to the 'low_in' for imadjust.
%   upperIntensityNorm (double): Normalized upper intensity input bound (0 to 1).
%                                This corresponds to the 'high_in' for imadjust.
%   gamma (double): Gamma correction factor.
%   outputMaxValue (double): The maximum value of the inputVolume's data type
%                            (e.g., 65535 for uint16). Used to scale normalized
%                            intensities if inputVolume is integer type.
%   wb_handle (handle, optional): Handle to a waitbar. If provided, the
%                                 waitbar will be updated during processing.
%
% Returns:
%   stretchedVolume (numeric array): The contrast-stretched 3D volume,
%                                    same type as inputVolume.

    stretchedVolume = zeros(size(inputVolume), 'like', inputVolume);
    nSlices = size(inputVolume, 3);

    if nargin < 6
        wb_handle = []; % No waitbar if not provided
    end

    for i = 1:nSlices
        % imadjust handles scaling for integer types if input is [0,1] range
        stretchedVolume(:, :, i) = imadjust(inputVolume(:, :, i), [lowerIntensityNorm, upperIntensityNorm], [0, 1], gamma);
        % Update waitbar if handle is valid
        if ~isempty(wb_handle) && isvalid(wb_handle)
            waitbar(i / nSlices, wb_handle, sprintf('Stretching slice %d of %d...', i, nSlices));
        end
    end
end