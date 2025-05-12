function visualizeSlices(trainVolume, labelVolume)
    dims = size(trainVolume);
    nSlice = dims(3);
    
    for i = 1:nSlice
        % Extract the i-th slice from both volumes
        trainSlice = trainVolume(:, :, i);
        labelSlice = labelVolume(:, :, i);
        
        % Visualize
        figure(1); clf;
        subplot(1, 2, 1);
        imshow(trainSlice);
        title(['Train Slice', num2str(i)]);
        
        subplot(1, 2, 2);
        imshow(labelSlice,[]);
        title(['Label Slice', num2str(i)]);
        
        pause(0.01);
    end
end