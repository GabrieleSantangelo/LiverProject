function [imageVolume, labelVolume] = loadNiiVolumes(imagePath, labelPath)
% loadNiiVolumes Loads NIfTI image and label volumes.
%   Handles .gz decompression if files are compressed.
%
% Args:
%   imagePath (char): Path to the image NIfTI file (e.g., 'image.nii.gz').
%   labelPath (char): Path to the label NIfTI file (e.g., 'label.nii.gz').
%
% Returns:
%   imageVolume (numeric array): Loaded image data (e.g., 3D or 4D array).
%   labelVolume (numeric array): Loaded label data.

    % Decompress if necessary
    uncompressedImagePath = attemptDecompression(imagePath);
    uncompressedLabelPath = attemptDecompression(labelPath);

    % Load the .nii files
    imageNii = load_nii(uncompressedImagePath);
    labelNii = load_nii(uncompressedLabelPath);

    % Extract image data
    imageVolume = imageNii.img;
    labelVolume = labelNii.img;

end

function uncompressedPath = attemptDecompression(filePath)
    [~, ~, ext] = fileparts(filePath);
    uncompressedPath = filePath;
    if strcmpi(ext, '.gz')
        uncompressedPath = filePath(1:end-3); % Path without .gz
        if ~isfile(uncompressedPath)
            gunzip(filePath);
        end
    end
end