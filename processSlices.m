function [meanValue, hMeanNorm, hMean_clean, normalizedSlice] = processSlices(TrainVolume, nBins, roiParams)
    dims = size(TrainVolume);
    nSlice = dims(3);
    
    % Prealloca un array 3D per salvare ogni slice elaborata
    allSlices = zeros(size(TrainVolume,1), size(TrainVolume,2), nSlice, 'uint16');
    partialMeanValues = zeros(nSlice, 1);
    
    hTotal = zeros(nBins,1);
    meanValue = 0;
    
    % Elaborazione di ogni slice e calcolo dei parametri intermedi
    for i = 1:nSlice
        % Estrai la slice originale
        sliceOriginale = TrainVolume(:, :, i);
        
        % Shift e normalizzazione dell'immagine
        TrainSliceShifted = sliceOriginale + 1024;
        TrainSliceNorm = double(TrainSliceShifted) / 65535;
        in_low = min(TrainSliceNorm(:));
        in_high = max(TrainSliceNorm(:));
        TrainSliceStretched = imadjust(TrainSliceNorm, [in_low in_high], [0 1]);
        TrainSliceStretched_uint16 = uint16(TrainSliceStretched * 65535);
        
        % Salva la slice elaborata
        allSlices(:, :, i) = TrainSliceStretched_uint16;
        
        % Calcola l'istogramma per la slice (opzionale: istogramma non normalizzato)
        h = imhist(TrainSliceStretched_uint16, nBins);
        hTotal = hTotal + h;

        % Calcola la media dei pixel nella ROI circolare per la slice corrente
        partialMeanValues(i) = computeROIMean(TrainSliceStretched_uint16, roiParams);
        
        % Accumula la somma per la media globale
        meanValue = meanValue + partialMeanValues(i);
    end
    % Calcola la media globale nelle ROI
    meanValue = meanValue / nSlice;
    
    % Istogramma medio (sulle slice elaborate non normalizzate)
    hMean = hTotal / nSlice;
    
    % Rimozione outlier nell'istogramma medio (se necessario)
    hMean_clean = removeHistogramOutliers(hMean);
    
    % Normalizzazione delle slice: applica il fattore di correzione per ogni slice
    normalizedSlice = zeros(size(allSlices), 'uint16');
    for i = 1:nSlice
        factor = meanValue / partialMeanValues(i);
        % Converti la slice in double per eseguire l'operazione
        normalized = double(allSlices(:, :, i)) .* factor;
        % Se necessario, riconverti in uint16
        normalizedSlice(:, :, i) = uint16(normalized);
    end
    
    % Calcola l'istogramma sulle slice normalizzate
    hTotalNorm = zeros(nBins,1);
    for i = 1:nSlice
        h_norm = imhist(normalizedSlice(:, :, i), nBins);
        hTotalNorm = hTotalNorm + h_norm;
    end
    hMeanNorm = hTotalNorm / nSlice;
end
