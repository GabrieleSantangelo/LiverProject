function plotHistograms(hMean, hMean_clean, nBins)
    figure;
    subplot(2,1,1);
    bar(hMean);
    title('Istogramma medio originale');
    xlabel('Intensità');
    ylabel('Frequenza media');
    
    subplot(2,1,2);
    bar(hMean_clean);
    title('Istogramma medio con outlier rimossi');
    xlabel('Intensità');
    ylabel('Frequenza media');
    
    % Visualizzazione aggiuntiva dell'istogramma 10 bit dell'ultima slice elaborata
    figure;
    histogram(uint16(hMean), nBins);
    title('Istogramma della slice (10 bit)');
    
end