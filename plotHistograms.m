function plotHistograms(hMean, hMean_clean)
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
    
end