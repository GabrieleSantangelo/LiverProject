function [metrics, predictionVolume] = calculateSegmentationMetrics(labelVolume, liverMaskPred, tumorMaskPred)
% calculateSegmentationMetrics Calculates Dice scores and recall for liver and tumor.
%
% Args:
%   labelVolume (numeric array): Ground truth volume (0:bg, 1:liver, 2:tumor).
%   liverMaskPred (logical array): Predicted 3D binary mask for the liver.
%   tumorMaskPred (logical array): Predicted 3D binary mask for tumors.
%
% Returns:
%   metrics (struct): Struct containing 'DiceLiver', 'DiceTumor', 
%                     'RecallLiver', 'RecallTumor'.
%   predictionVolume (numeric array): Combined prediction volume (0,1,2).

    % Prepare Ground Truth Masks
    liver_mask_gt = (labelVolume == 1) | (labelVolume == 2); % Liver GT includes tumor area
    tumor_mask_gt = (labelVolume == 2);

    % Construct Combined Prediction Volume
    predictionVolume = zeros(size(labelVolume), 'like', labelVolume);
    predictionVolume(liverMaskPred) = 1; % Liver prediction
    predictionVolume(tumorMaskPred) = 2; % Tumor prediction (overwrites liver if overlapping)

    % Calculate Dice Scores using MATLAB's dice function
    % Ensure inputs are logical for dice()
    metrics.DiceLiver = dice(liverMaskPred, liver_mask_gt);
    metrics.DiceTumor = dice(tumorMaskPred, tumor_mask_gt);
    if sum(tumor_mask_gt(:)) == 0 && sum(tumorMaskPred(:)) == 0 % Both GT and Pred are empty for tumor
        metrics.DiceTumor = 1.0; % Or NaN, depending on desired behavior for empty sets
    elseif sum(tumor_mask_gt(:)) == 0 && sum(tumorMaskPred(:)) > 0 % GT empty, Pred not
        metrics.DiceTumor = 0.0;
    end

    % Calculate Recall (Sensitivity)
    metrics.RecallLiver = sum(liverMaskPred(:) & liver_mask_gt(:)) / sum(liver_mask_gt(:));
    if sum(tumor_mask_gt(:)) > 0
        metrics.RecallTumor = sum(tumorMaskPred(:) & tumor_mask_gt(:)) / sum(tumor_mask_gt(:));
    else
        metrics.RecallTumor = NaN; % Or 1.0 if tumorMaskPred is also empty
    end
end