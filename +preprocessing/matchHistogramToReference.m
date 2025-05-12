function matchedVolume = matchHistogramToReference(inputVolume, referenceSliceIndex, wb_handle)
% matchHistogramToReference Matches histogram of each slice to a reference slice.
%
% Args:
%   inputVolume (numeric array): The 3D input volume.
%   referenceSliceIndex (integer): Index of the slice to be used as reference.
%                                  If empty or not provided, matches slice i to i-1.
%   wb_handle (handle, optional): Handle to a waitbar. If provided, the
%                                 waitbar will be updated during processing.
%
% Returns:
%   matchedVolume (numeric array): Volume with histograms matched.

    matchedVolume = inputVolume;
    nSlices = size(inputVolume, 3);

    if nargin < 3
        wb_handle = []; % No waitbar if not provided
    end

    if nargin < 2 || isempty(referenceSliceIndex)
        % Sequential matching (slice i to slice i-1)
        for i = 2:nSlices
            currentSlice = inputVolume(:, :, i);
            referenceSlice = matchedVolume(:, :, i-1); % Use previously matched slice
            % nBins for imhistmatch can be tricky; default 256 or adaptive
            matchedVolume(:, :, i) = imhistmatch(currentSlice, referenceSlice);

            % Update waitbar if handle is valid
            if ~isempty(wb_handle) && isvalid(wb_handle)
                waitbar(i / nSlices, wb_handle, sprintf('Matching histogram for slice %d of %d...', i, nSlices));
            end
        end
    else
        % Match all slices to a single reference slice
        referenceSlice = inputVolume(:, :, referenceSliceIndex);
        for i = 1:nSlices
            matchedVolume(:, :, i) = imhistmatch(inputVolume(:, :, i), referenceSlice);

            % Update waitbar if handle is valid
            if ~isempty(wb_handle) && isvalid(wb_handle)
                waitbar(i / nSlices, wb_handle, sprintf('Matching histogram for slice %d of %d...', i, nSlices));
            end
        end
    end
end