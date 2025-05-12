function openedVolume = openVolume(binaryVolume, structuringElement)
% openVolume Performs a 3D morphological opening on a binary volume.
%
% Args:
%   binaryVolume (logical array): The input 3D binary volume.
%   structuringElement (strel object): The structuring element for the opening.
%                                      Should be a 3D strel (e.g., strel('sphere', r)).
%
% Returns:
%   openedVolume (logical array): The morphologically opened volume.

    openedVolume = imopen(binaryVolume, structuringElement);

end