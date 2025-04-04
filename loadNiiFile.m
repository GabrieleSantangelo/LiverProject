function [trainVolume, labelVolume] = loadNiiFile(imagesPath, labelsPath)

    % Check if the files are already uncompressed; if not, decompress them
    if ~isfile(imagesPath(1:end-3))
        gunzip(imagesPath);
    end
    if ~isfile(labelsPath(1:end-3))
        gunzip(labelsPath);
    end

    % Extract the names of the decompressed files (without the .gz extension)
    [filepathX, nameX, ~] = fileparts(imagesPath);
    [filepathY, nameY, ~] = fileparts(labelsPath);

    niiFileImages = fullfile(filepathX, nameX);
    niiFileLabels = fullfile(filepathY, nameY);
    

    % Load the two .nii files
    trainNII = load_nii(niiFileImages);
    labelNII = load_nii(niiFileLabels);

    % Extract the image data from the loaded NIfTI structures
    trainVolume = trainNII.img;
    labelVolume = labelNII.img;
end