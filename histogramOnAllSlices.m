function [hMean, hMean_clean] = histogramOnAllSlices(normalizedSlice, maxValue)
    dims = size(normalizedSlice);
    nSlice = dims(3);
    
    hTotal = zeros(maxValue,1);
    
    % Elaborazione di ogni slice e calcolo dei parametri intermedi
    for i = 1:nSlice
        % Calcola l'istogramma per la slice (opzionale: istogramma non normalizzato)
        h = imhist(normalizedSlice, maxValue);
        hTotal = hTotal + h;
    end
    % Istogramma medio (sulle slice elaborate non normalizzate)
    hMean = hTotal / nSlice;
    
    % Rimozione outlier nell'istogramma medio (se necessario)
    hMean_clean = removeHistogramOutliers(hMean);
end
