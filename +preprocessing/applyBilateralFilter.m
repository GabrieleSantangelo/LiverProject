function filteredVolume = applyBilateralFilter(inputVolume, degreeOfSmoothing, spatialSigma, maxValue, wb_handle)
% applyBilateralFilter Applies a bilateral filter to each slice of a 3D volume.
%
% Args:
%   inputVolume (numeric array): The input 3D volume (e.g., uint16).
%   degreeOfSmoothing (double): Degree of smoothing for imbilatfilt.
%                               Controls intensity similarity.
%   spatialSigma (double): Spatial standard deviation for imbilatfilt.
%                          Controls spatial extent of the filter.
%   maxValue (double): The maximum value of the inputVolume's data type
%                      (e.g., 65535 for uint16). Used to normalize to [0,1]
%                      for imbilatfilt if input is integer.
%   wb_handle (handle, optional): Handle to a waitbar.
%
% Returns:
%   filteredVolume (numeric array): The bilaterally filtered 3D volume,
%                                   same type as inputVolume.

    filteredVolume = zeros(size(inputVolume), 'like', inputVolume);
    nSlices = size(inputVolume, 3);

    if nargin < 5
        wb_handle = []; % No waitbar if not provided
    end

    for i = 1:nSlices
        % Ensure both operands for division are double to avoid type conflicts.
        % maxValue is expected to be a scalar.
        normalizedSlice = double(inputVolume(:,:,i)) ./ double(maxValue); % Normalize to [0,1]
        filteredNormalizedSlice = imbilatfilt(normalizedSlice, degreeOfSmoothing, spatialSigma);
        % For scaling back, ensure maxValue is also double if it wasn't already, for consistency.
        filteredVolume(:,:,i) = cast(filteredNormalizedSlice .* double(maxValue), class(inputVolume)); % Scale back

        % Update waitbar if handle is valid
        if ~isempty(wb_handle) && isvalid(wb_handle)
            waitbar(i / nSlices, wb_handle, sprintf('Applying Bilateral Filter (Slice %d of %d)...', i, nSlices));
        end
    end
end