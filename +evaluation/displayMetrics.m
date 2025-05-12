function displayMetrics(metrics)
% displayMetrics Displays segmentation metrics.
%
% Args:
%   metrics (struct): Struct containing fields like 'DiceLiver', 'DiceTumor', etc.

    fprintf('\n--- Segmentation Metrics ---\n');
    fprintf('Dice Score Liver: %.4f\n', metrics.DiceLiver);
    fprintf('Dice Score Tumor: %.4f\n', metrics.DiceTumor);
    fprintf('Recall Liver:     %.4f\n', metrics.RecallLiver);
    fprintf('Recall Tumor:     %.4f\n', metrics.RecallTumor);
    fprintf('---------------------------\n');
end