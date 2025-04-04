function meanVal = computeROIMean(image, roiParams)
    % Crea le matrici di coordinate per i pixel
    [X, Y] = meshgrid(1:size(image,2), 1:size(image,1));
    mask = (X - roiParams.x).^2 + (Y - roiParams.y).^2 <= roiParams.r.^2;
    circleValues = image(mask);
    meanVal = mean(circleValues);
end