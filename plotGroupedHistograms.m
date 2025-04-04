function plotGroupedHistograms(binCenters, grouped_hMean, grouped_hMean_clean, nBins, lowerIntensity, upperIntensity)
    figure;
    subplot(2,1,1);
    bar(binCenters, grouped_hMean);
    title('Istogramma medio raggruppato per ampiezza 1000 - Originale');
    xlabel('Intensità');
    ylabel('Frequenza media');
    
    subplot(2,1,2);
    bar(binCenters, grouped_hMean_clean);
    title('Istogramma medio raggruppato per ampiezza 1000 - Con outlier rimossi');
    xlabel('Intensità');
    ylabel('Frequenza media');
    xline(lowerIntensity*nBins, 'r--', 'LineWidth', 2);
    xline(upperIntensity*nBins, 'r--', 'LineWidth', 2);

end