function render3DSegmentation(liverMask, tumorMask, bonesMask, config)
% render3DSegmentation Renders 3D surfaces for liver, tumor, and bones.
%
% Args:
%   liverMask (logical array): 3D binary mask for the liver.
%   tumorMask (logical array): 3D binary mask for tumors.
%   bonesMask (logical array, optional): 3D binary mask for bones.
%   config (struct): Configuration for rendering.
%       .liverColor, .tumorColor, .bonesColor (1x3 double): RGB colors.
%       .liverAlpha, .tumorAlpha, .bonesAlpha (double): Alpha transparency.
%       .smoothing (logical): If true, apply smooth3.
%       .smoothFactor (integer): Gaussian smoothing factor for smooth3.

    figure('Name', '3D Segmentation Render', 'NumberTitle', 'off');
    hold on;

    % Render Liver
    if any(liverMask(:))
        renderSurface(liverMask, config.liverColor, config.liverAlpha, config.smoothing, config.smoothFactor, 'Liver');
    end

    % Render Tumor
    if any(tumorMask(:))
        renderSurface(tumorMask, config.tumorColor, config.tumorAlpha, config.smoothing, config.smoothFactor, 'Tumor');
    end

    % Render Bones (if provided)
    if nargin > 2 && ~isempty(bonesMask) && any(bonesMask(:))
        renderSurface(bonesMask, config.bonesColor, config.bonesAlpha, config.smoothing, config.smoothFactor, 'Bones');
    end

    hold off;
    title('3D Reconstruction');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    grid on; axis equal; view(3);
    camlight; lighting gouraud; rotate3d on;
end

function renderSurface(mask, color, alpha, smoothing, smoothFactor, name)
    fprintf('Rendering %s...\n', name);
    mask_double = double(mask);
    if smoothing
        % smooth3 can be slow on large volumes, consider patch('SmoothMesh') if available
        % or reducing smoothFactor for speed.
        mask_double = smooth3(mask_double, 'gaussian', smoothFactor, smoothFactor/2); % Use smoothFactor also for std dev
    end
    
    % Isosurface level 0.5 is standard for binary masks
    [faces, verts] = isosurface(mask_double, 0.5);
    
    if ~isempty(faces) && ~isempty(verts)
        h_patch = patch('Faces', faces, 'Vertices', verts);
        set(h_patch, 'FaceColor', color, 'EdgeColor', 'none', 'FaceAlpha', alpha);
    else
        fprintf('No surface to render for %s at isosurface level 0.5.\n', name);
    end
    fprintf('%s rendering complete.\n', name);
end