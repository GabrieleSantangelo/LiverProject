function [totalMeanValue, normalizedSlice] = normalizingSlices(TrainVolume, roiParams, maxValue)
    dims = size(TrainVolume);
    nSlice = dims(3);

    log_base = @(x, base) log(x) / log(base);
    
    % Prealloca un array 3D per salvare ogni slice elaborata
    allSlices = zeros(size(TrainVolume,1), size(TrainVolume,2), nSlice, 'uint16');
    partialMeanValues = zeros(nSlice, 1);

    partialMaxValues = zeros(nSlice, 1);

    
    totalMeanValue = 0;
    totalMaxValue = 0;
    
    % Elaborazione di ogni slice e calcolo dei parametri intermedi
    for i = 1:nSlice
        % Estrai la slice originale
        sliceOriginale = TrainVolume(:, :, i);
        
        % Shift e normalizzazione dell'immagine
        TrainSliceShifted = sliceOriginale + 1024;
        TrainSliceNorm = double(TrainSliceShifted) / maxValue;
        in_low = min(TrainSliceNorm(:));
        in_high = max(TrainSliceNorm(:));
        TrainSliceStretched = imadjust(TrainSliceNorm, [in_low in_high], [0 1]);
        TrainSliceStretched_uint16 = uint16(TrainSliceStretched * maxValue);
        
        % Salva la slice elaborata
        allSlices(:, :, i) = TrainSliceStretched_uint16;

        % Calcola la media dei pixel nella ROI circolare per la slice corrente
        [partialMeanValues(i), partialMaxValues(i)] = computeROIMeanMax(TrainSliceStretched_uint16, roiParams);
        
        % Accumula la somma per la media globale
        totalMeanValue = totalMeanValue + partialMeanValues(i);
        totalMaxValue = totalMaxValue + partialMaxValues(i);

    end
    % Calcola la media globale nelle ROI
    totalMeanValue = totalMeanValue / nSlice;
    totalMaxValue = totalMaxValue /nSlice;
    
    % Normalizzazione delle slice: applica il fattore di correzione per ogni slice
    normalizedSlice = zeros(size(allSlices), 'uint16');
    for i = 1:nSlice
        ratio = totalMaxValue / partialMaxValues(i);

        adjPartialMean = partialMeanValues(i)*ratio;

        gamma = log_base(totalMeanValue ,adjPartialMean);

        % factor = totalMeanValue / partialMeanValues(i);
        % Converti la slice in double per eseguire l'operazione
        normalized = (double(allSlices(:, :, i)) .* ratio) .^ gamma;
        % Se necessario, riconverti in uint16
        normalizedSlice(:, :, i) = uint16(normalized);
    end
end
