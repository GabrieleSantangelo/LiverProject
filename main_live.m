clear; clc; clf;
disp("Code Run")

% Option to visualize slices
visualizeSlicesFlag = false;


% Load the slices from files
[trainVolume, labelVolume] = loadNiiFile(...
    'data/imagesTr/liver_80.nii.gz', ...
    'data/labelsTr/liver_80.nii.gz'  ...
);


% Visualize the slices
if visualizeSlicesFlag
    visualizeSlices(trainVolume, labelVolume);
end

%%
%[text] Define the region of interest (ROI) parameters for cross normalization
roiParams.x = 35;
roiParams.y = 38;
roiParams.r = 30;


% Maximum value for uint16
maxValue = 65536;

% ----------- START REFACTOR FOLLOWING CODE ----------- %

% Number of bins for histogram
nBins = 65536; % Number of bins for histogram
[meanValue, normalizedSlice] = normalizingSlices(trainVolume, roiParams, maxValue);

% normalizedSlice =  histogramMachingAllSlice(normalizedSlice_temp, 150);
[hMean, hMean_clean] = histogramOnAllSlices(normalizedSlice, nBins);

% visualizeSlices(normalizedSlice, normalizedSlice_temp)
%%
fprintf('Mean value in ROI: %.2f\n', meanValue);

% Visualize mean histogram of all slices
plotHistograms(hMean, hMean_clean);

% Group histogram data
groupSize = 2000;

[grouped_hMean, grouped_hMean_clean, binCenters] = groupHistogramData(hMean, hMean_clean, groupSize, nBins);
[lowerIntensity, upperIntensity] = bandDetection(grouped_hMean_clean); % 30000 con grouped_hMean_clean | 3000000 con grouped_hMean

plotGroupedHistograms(binCenters, grouped_hMean, grouped_hMean_clean, nBins, lowerIntensity, upperIntensity);
%%
% ------------ STOP REFACTOR FOLLOWING CODE ----------- %
% Stretch the slices using the lower and upper intensity values
stretchedSlice = stretchSlices(normalizedSlice, lowerIntensity, upperIntensity, 8);

dims = size(trainVolume);
nSlice = dims(3);

no_bones_slice = zeros(size(trainVolume,1), size(trainVolume,2), nSlice,'uint16');


   
figure;
imshow(stretchedSlice(:,:,195))
title("Porco caneeeeee");
window = 10
for slice_idx=1:nSlice-window
    
    tempSlice = normalizedSlice(:,:,slice_idx);

    % Remove high value (in this case mainly areas representing the bones)
    for j= 1:window
    mask = stretchedSlice(:,:,slice_idx+j-1) > maxValue * 0.95;
    mask = imdilate(double(mask), strel("disk", 4)); % Dilation to remove small holes
    mask = imerode(mask, strel("disk", 2)); % Erosion to remove small holes
    mask = imfill(mask, "holes"); % Fill holes
    mask = imfilter(mask, fspecial("gaussian", 10)); % Gaussian filter to smooth the mask

    tempSlice = uint16(double(tempSlice) .* double(1 - mask));
    end
   
    no_bones_slice(:,:,slice_idx) = tempSlice;
    figure(10);clf;
    imshow(no_bones_slice(:,:,slice_idx))
    pause(0.001)
end


% Visualizzazione della slice modificata
figure;
imshow(no_bones_slice(:,:,slice_idx));
%%

[hMean_no_bones, hMean_clean_no_bones] = histogramOnAllSlices(no_bones_slice, maxValue);

[grouped_hMean_no_bones, grouped_hMean_clean_no_bones, binCenters] = groupHistogramData(hMean_no_bones, hMean_clean_no_bones, groupSize, nBins);


[lowerIntensity_no_bones, upperIntensity_no_bones] = bandDetection(grouped_hMean_clean_no_bones);

plotGroupedHistograms(binCenters, grouped_hMean_no_bones, grouped_hMean_clean_no_bones, nBins, lowerIntensity_no_bones, upperIntensity_no_bones);


%%
flag = true

doubleStretchedSlice_temp = stretchSlices(no_bones_slice, lowerIntensity_no_bones, upperIntensity_no_bones, 6);

%mask_after_no_bones = doubleStretchedSlice_temp >= lowerIntensity_no_bones*maxValue*0.1 & doubleStretchedSlice_temp <= upperIntensity_no_bones*maxValue*1.1;
%doubleStretchedSlice_temp = doubleStretchedSlice_temp.*mask_after_no_bones
for i=1:217
    figure(12);clf;
    imshow(doubleStretchedSlice_temp(:,:,i));
    pause(0.01)
end
% doubleStretchedSlice_temp(1:150,: ,:) = 0;
%%
for slice_idx=1:nSlice


    % Remove high value (in this case mainly areas representing the bones)
    mask = doubleStretchedSlice_temp(:,:,slice_idx) > maxValue * 0.13;
    

    mask = imdilate(double(mask), strel("disk", 4)); % Dilation to remove small holes
    mask = imerode(mask, strel("disk", 4)); % Erosion to remove small holes
    mask = imfill(mask, "holes"); % Fill holes
    mask_array(:,:,slice_idx) = 1 - imfilter(mask, fspecial("gaussian", 10)); % Gaussian filter to smooth the mask

    % tempSlice = uint16(double(doubleStretchedSlice_temp(:,:,slice_idx)) .* double(1 - mask));
    % doubleStretchedSlice_temp2(:,:,slice_idx) = tempSlice;
    
    figure(1);clf;
    imshow(1 - mask_array(:,:,slice_idx))
    title(slice_idx)

    pause(0.001)
end

%%
doubleStretchedSlice_temp2 = histogramMachingAllSlice(doubleStretchedSlice_temp, 150);
% fino a qua è buono
%%
for slice_idx=1:nSlice
    figure(4483); clf;
    doubleStretchedSlice_final(:,:,slice_idx)= uint16(double(doubleStretchedSlice_temp2(:,:,slice_idx)) .* double(1 - mask_array(:,:,slice_idx)));
    imshow(doubleStretchedSlice_final(:,:,slice_idx))
    pause(0.0001)
end
%%
DoS = 0.1; % Regola questo (Degree of Smoothing)
sSigma = 10;

h = waitbar(0);
for slice_idx=1:nSlice
     waitbar(slice_idx/nSlice, h);
    % tempSlice = uint16(double(doubleStretchedSlice_final(:,:,slice_idx)) .* double(1 - mask_array(:,:,slice_idx)));
    % doubleStretchedSlice(:,:,slice_idx) = kuwahara(tempSlice, 19);
    doubleStretchedSlice(:,:,slice_idx) = imbilatfilt(double(double(doubleStretchedSlice_final(:,:,slice_idx))./maxValue), DoS, sSigma);
end

%%
%[text] ## Bilaterale \-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-
% Assumendo che l'input sia 'img_input' (es. tempSlice o l'output di Kuwahara)
% Assicurati che sia double e preferibilmente in scala [0, 1]

DoS = 0.6; % Regola questo (Degree of Smoothing)
sSigma = 5; % Regola questo (Spatial Sigma)

for slice_idx=1:nSlice
smoothed_slice = imbilatfilt(doubleStretchedSlice(:,:,slice_idx), DoS, sSigma);

% 'smoothed_slice' è il risultato da usare per la segmentazione
figure(103); clf;
doubleStretchedSlice(:,:,slice_idx) = smoothed_slice; 
imshow(doubleStretchedSlice(:,:,slice_idx))
title("BILATERALEEEEEEEEEEE")
pause(0.001)
end

%%

% imshow(doubleStretchedSlice(:,:,180))

if flag
    for slice_idx=1:nSlice
    
        figure(1); clf;
        subplot(1, 2, 1);
        imshow(doubleStretchedSlice(:,:,slice_idx))
        title(['Label Slice', num2str(slice_idx)]);
    
        subplot(1, 2, 2);
        imshow(doubleStretchedSlice_temp(:,:,slice_idx))
        title(['Label Slice', num2str(slice_idx)]);
        pause(0.0001);
    end
end


%%
%------------- Segmentazione Morfologica --------------%

slice_idx = 163;
slice = doubleStretchedSlice(:,:,slice_idx);
figure;
imshow(slice, []);
title('Disegna una ROI nel fegato');
roi = drawpolygon; % O drawpolygon, drawellipse...
wait(roi); % Attendi che l'utente finisca di disegnare
mask = createMask(roi);
liver_pixels = slice(mask); % Estrai i pixel della ROI
figure;
imhist(liver_pixels);
title('Istogramma della ROI nel fegato');



%%
% --- Soglie identificate dall'analisi ROI ---
% Assumiamo che doubleStretchedSlice sia effettivamente in scala [0, 1]
% Se così non fosse, dovresti riscalare queste soglie.
lowerThreshold = 0.06;
upperThreshold = 0.3;

% --- Ottieni dimensioni del volume ---
dims = size(doubleStretchedSlice);
nRows = dims(1);
nCols = dims(2);
nSlice = dims(3); % Numero di slice nel volume

% --- Inizializza il volume della maschera binaria ---
% Usiamo 'logical' per efficienza di memoria e calcolo
initialLiverMask = false(nRows, nCols, nSlice);

disp('Applicazione thresholding slice per slice...');
tic; % Inizia a misurare il tempo

% --- Applica il thresholding a ogni slice ---
for slice_idx = 1:nSlice
    % Estrai la slice corrente (è già di tipo double)
    currentSlice = doubleStretchedSlice(:,:,slice_idx);
    
    % Crea la maschera binaria per la slice corrente:
    % Metti a 'true' (1) i pixel con intensità compresa tra le soglie
    mask_slice = (currentSlice >= lowerThreshold) & (currentSlice <= upperThreshold);

    % mask_slice_erode = imerode(mask_slice, strel("disk", 5));
    % mask_slice_dilate = imdilate(mask_slice_erode, strel("disk", 3));
    
    % Salva la maschera binaria nella posizione corretta del volume 3D
    mask_fill = imfill(mask_slice, 18 ,"holes");
    
    mask_opened = imopen(mask_fill, strel('disk', 5));

    % mask_slice_erode = imerode(mask_fill, strel("disk", 6));
    % mask_slice_dilate = imdilate(mask_slice_erode, strel("disk", 3));
    initialLiverMask(:,:,slice_idx) = mask_opened;
    
    % --- (Opzionale) Visualizzazione per controllo ---q
         figure(300); clf;
         subplot(1,2,1); imshow(currentSlice); title(['Slice Filtrata ', num2str(slice_idx)]);
         subplot(1,2,2); imshow(mask_opened); title('Maschera da Thresholding');
         drawnow; % Forza l'aggiornamento della figura
         pause(0.1); % Pausa opzionale
end

elapsedTime = toc; % Ferma il cronometro
fprintf('Thresholding completato in %.2f secondi.\n', elapsedTime);

% Ora 'initialLiverMask' è un volume 3D binario (logical)
% che rappresenta la prima stima della segmentazione del fegato.

% --- Visualizza un esempio di maschera ottenuta ---
slice_to_show = 180; % La stessa slice usata per la ROI o un'altra

figure;
subplot(1,2,1);
imshow(doubleStretchedSlice(:,:,slice_to_show), []);
title(['Slice Filtrata (Kuwahara) ', num2str(slice_to_show)]);
subplot(1,2,2);
imshow(initialLiverMask(:,:,slice_to_show));
title(['Maschera Iniziale (Thresholding) ', num2str(slice_to_show)]);



%%
% --- Raffinamento Morfologico ---
refinedLiverMask = false(size(initialLiverMask)); % Inizializza maschera finale
minLiverArea = 100; % AREA MINIMA in pixel -> REGOLA QUESTO VALORE!
se_close = strel('disk', 7); % Elemento strutturante per 'imclose' (opzionale)
se_open = strel('disk', 5); % Elemento strutturante per 'imopen' (opzionale)

disp('Inizio raffinamento morfologico...');
tic;
for slice_idx = 1:nSlice
    if mod(slice_idx, 20) == 0
        fprintf('Raffinamento slice %d/%d...\n', slice_idx, nSlice);
    end
    
    binarySlice = initialLiverMask(:,:,slice_idx);
    
    % 1. Riempimento buchi interni
    filledSlice = imfill(binarySlice, 'holes');
    
    % 2. (Opzionale) Chiusura morfologica: Chiude piccoli gap
    closedSlice = imclose(filledSlice, se_close);
    
    % 3. (Opzionale) Apertura morfologica: Rimuove piccole connessioni/oggetti
    openedSlice = imopen(closedSlice, se_open); % Se usi la chiusura
    % openedSlice = imopen(filledSlice, se_open); % Se NON usi la chiusura

    % 4. Rimuovi oggetti piccoli (rumore, vasi, ecc.)
    cleanedSlice = bwareaopen(openedSlice, minLiverArea);
    
    % 5. Tieni solo la componente connessa più grande (solitamente il fegato)
    CC = bwconncomp(cleanedSlice);
    if CC.NumObjects > 0
        stats = regionprops(CC, 'Area');
        [~, idx] = max([stats.Area]); % Indice della componente più grande
        finalSlice = false(size(binarySlice));
        finalSlice(CC.PixelIdxList{idx}) = true; % Crea maschera solo con quella
    else
        finalSlice = false(size(binarySlice)); % Nessun oggetto rimasto
    end
    
    refinedLiverMask(:,:,slice_idx) = finalSlice; % Salva la maschera raffinata
end
toc;
disp('Raffinamento morfologico completato.');

% --- Visualizzazione Finale ---
for slice_to_show=1:nSlice
% slice_to_show = 180;
figure(1); clf;
imshow(doubleStretchedSlice(:,:,slice_to_show), []);
hold on;
% Usa visboundaries per mostrare il contorno della segmentazione
visboundaries(refinedLiverMask(:,:,slice_to_show), 'Color', 'r', 'LineWidth', 1);
hold off;
title(['Segmentazione Finale Fegato (Slice ', num2str(slice_to_show), ')']);
pause(0.01)

% Il volume 'refinedLiverMask' contiene la tua segmentazione finale del fegato.
end

%%
% --- Thresholding Applicato all'Intero Volume 3D ---
% lowerThreshold = 0.1; % Dalla tua analisi ROI
%upperThreshold = 0.25;
%disp('Applicazione thresholding 3D...');
%initialLiverMask_3D = (doubleStretchedSlice >= lowerThreshold) & (doubleStretchedSlice <= upperThreshold);

% --- Operazioni Morfologiche 3D ---
disp('Riempimento buchi 3D...');
% imfill funziona anche in 3D
filledVolume_3D = imfill(initialLiverMask, 'holes');

% (Opzionale) Apertura/Chiusura 3D - richiede elementi strutturanti 3D
se_3D = strel('sphere', 4); % Esempio: sfera di raggio 2
openedVolume_3D = imopen(filledVolume_3D, se_3D);
% cleanedVolume_3D = imclose(openedVolume_3D, se_3D); % Esempio

% --- Rimozione Oggetti Piccoli 3D ---
minLiverVolumeVoxels = 300000; % SOGLIA DI VOLUME IN VOXEL -> DA REGOLARE!
disp(['Rimozione oggetti 3D più piccoli di ', num2str(minLiverVolumeVoxels), ' voxels...']);
% Usa connettività 26 per considerare vicini diagonali in 3D
cleanedVolume_3D = bwareaopen(openedVolume_3D, minLiverVolumeVoxels, 18); 
                               % Usa filledVolume_3D o opened/closed se usati

% --- Trova Componenti Connesse 3D ---
disp('Ricerca componenti connesse 3D...');
CC_3D = bwconncomp(cleanedVolume_3D, 18); % Usa connettività 26

refinedLiverMask_3D = false(size(doubleStretchedSlice)); % Inizializza maschera finale

if CC_3D.NumObjects > 0
    disp('Selezione componente 3D più grande...');
    % Calcola il volume (numero di voxel) di ogni componente
    volumes = cellfun(@numel, CC_3D.PixelIdxList);
    % Trova l'indice della componente con il volume massimo
    [~, idxLargest] = max(volumes);

    % Crea la maschera finale solo con la componente più grande
    refinedLiverMask_3D(CC_3D.PixelIdxList{idxLargest}) = true;
    fprintf('Componente più grande selezionata (Volume: %d voxels).\n', volumes(idxLargest));
else
    disp('Nessuna componente trovata dopo la pulizia 3D.');
end

disp('Segmentazione 3D completata.');

for slice_to_show=1:nSlice
% slice_to_show = 180;
figure(121); clf;
subplot(1,2,1); imshow(doubleStretchedSlice(:,:,slice_to_show), []); title(['Slice Filtrata ', num2str(slice_to_show)]);
subplot(1,2,2); imshow(doubleStretchedSlice(:,:,slice_to_show), []); hold on;
visboundaries(refinedLiverMask_3D(:,:,slice_to_show),'Color','g', 'LineWidth', 1); hold off;
title(['Segmentazione Finale 3D (Slice ', num2str(slice_to_show), ')']);
pause(0.01);
end


%%
% --- NUOVO BLOCCO: Post-processing per selezionare la componente più bassa nelle slice ---
disp('Applicazione criterio componente più bassa per slice...');
finalProcessedMask = false(size(refinedLiverMask_3D)); % Maschera per il risultato finale
[nRows, nCols, nSlice] = size(refinedLiverMask_3D); % Ricalcola o usa le dimensioni già note

for slice_idx = 1:nSlice
    currentSliceMask = refinedLiverMask_3D(:,:,slice_idx);

    % Se la slice non contiene alcun pixel della maschera, passa alla successiva
    if ~any(currentSliceMask(:))
        continue;
    end

    % Trova le componenti connesse *in questa slice 2D* (usare conn. 8 per 2D)
    CC_2D = bwconncomp(currentSliceMask, 8);

    % Se c'è solo una componente (o nessuna), la slice va bene così com'è
    if CC_2D.NumObjects <= 1
        finalProcessedMask(:,:,slice_idx) = currentSliceMask;
    else
        % Ci sono più componenti disconnesse in questa slice 2D
        % Trova quella che si estende più in basso (indice di riga massimo)
        maxRowOverall = -1;
        idxLowestComponent = -1;

        for k = 1:CC_2D.NumObjects
            % Trova gli indici di riga dei pixel della componente k
            [rows, ~] = ind2sub([nRows, nCols], CC_2D.PixelIdxList{k});
            currentComponentMaxRow = max(rows); % La riga più bassa di questa componente

            if currentComponentMaxRow > maxRowOverall
                maxRowOverall = currentComponentMaxRow;
                idxLowestComponent = k;
            end
        end

        % Crea una maschera temporanea per la slice contenente solo la componente più bassa
        tempSliceMask = false(nRows, nCols);
        if idxLowestComponent > 0 % Assicurati che sia stato trovato un indice valido
             tempSliceMask(CC_2D.PixelIdxList{idxLowestComponent}) = true;
        end
        finalProcessedMask(:,:,slice_idx) = tempSliceMask;
    end
end
disp('Criterio componente più bassa applicato.');

% --- FINE NUOVO BLOCCO ---


% --- Visualizzazione (usa la maschera finale processata) ---
disp('Visualizzazione maschera finale processata...');
for slice_to_show=1:nSlice
    figure(121); clf; % Usa una nuova figura o sovrascrivi la precedente
    subplot(1,2,1);
    imshow(doubleStretchedSlice(:,:,slice_to_show), []);
    title(['Slice Filtrata ', num2str(slice_to_show)]);

    subplot(1,2,2);
    imshow(doubleStretchedSlice(:,:,slice_to_show), []);
    hold on;
    % Usa finalProcessedMask per la visualizzazione
    visboundaries(finalProcessedMask(:,:,slice_to_show),'Color','b', 'LineWidth', 1); % Cambiato colore in blu per distinguerla
    hold off;
    title(['Segmentazione Finale Processata (Slice ', num2str(slice_to_show), ')']);

    pause(0.01); % Riduci la pausa se troppo lenta
end
%%
% ... dopo aver ottenuto refinedLiverMask_3D con fegato+cuore ...

disp('Calcolo centroide 3D per taglio verticale...');
props = regionprops3(refinedLiverMask_3D, 'Centroid');

if ~isempty(props)
    centroid_3D = props.Centroid(1,:); % Assumendo ci sia solo 1 oggetto nella maschera
    centroid_Z_slice = round(centroid_3D(3)); % Coordinata Z (indice slice) del centroide

    fprintf('Centroide Z trovato alla slice: %d. Rimuovo maschera sopra questa slice.\n', centroid_Z_slice);

    % Rimuovi tutto dalla slice del centroide in su (o regola +1/-1 se necessario)
    refinedLiverMask_3D(:,:, centroid_Z_slice:end) = false; 

else
    disp('Attenzione: impossibile calcolare il centroide (maschera vuota?).');
end


for slice_to_show=1:nSlice
% slice_to_show = 180;
figure(121); clf;
subplot(1,2,1); imshow(doubleStretchedSlice(:,:,slice_to_show), []); title(['Slice Filtrata ', num2str(slice_to_show)]);
subplot(1,2,2); imshow(doubleStretchedSlice(:,:,slice_to_show), []); hold on;
visboundaries(refinedLiverMask_3D(:,:,slice_to_show),'Color','g', 'LineWidth', 1); hold off;
title(['Segmentazione Finale 3D (Slice ', num2str(slice_to_show), ')']);
pause(0.01);
end


%------------- FINE Segmentazione Morfologica --------------%
%%
for slice_idx=1:nSlice
    figure(410); clf;
    imshow(imfill(imclose(edge(imfilter(doubleStretchedSlice(:,:,slice_idx), fspecial("average", 10))), strel("diamond", 10)), "holes"))
    pause(0.0001)
end


%%
for slice_idx=1:nSlice
    figure(415); clf;

    slice = doubleStretchedSlice(:,:,slice_idx) > 0.1;

    % slice = imdilate(slice, strel("diamond", 5));
    slice = imfill(slice, "holes");
    slice = imerode(slice, strel("disk", 8));
    slice = imdilate(slice, strel("diamond", 5));
    slice = imfill(slice, 18 ,"holes");
    slice = imerode(slice, strel("disk", 5));
    slice = imfill(slice, "holes");
    slice = imerode(slice, strel("disk", 5));
    slice = imfill(slice, "holes");
    slice = imerode(slice, strel("disk", 7));
    slice = imfill(slice, "holes");
    % slice = imfilter(double(imfilter(double(slice), fspecial("average", 1)) > 0.5), fspecial("average", 1)) > 0.5
    imshow(slice)

    pause(0.0001)
end

%%

sliceDouble = double(doubleStretchedSlice);
mu = mean(sliceDouble(:));
sigma = std(sliceDouble(:));

lower = mu + 0.7*sigma;
upper = mu + 3*sigma;

% temp_1 = doubleStretchedSlice > maxValue / 14 & doubleStretchedSlice < maxValue / 2; %14 e 2
% temp_1 = sliceDouble > lower & sliceDouble < upper;
for slice_idx=1:nSlice

    level = graythresh(doubleStretchedSlice(:,:,slice_idx));
    temp_1(:,:,slice_idx) = imbinarize(doubleStretchedSlice(:,:,slice_idx), level);

    figure(1); clf;
    imshow(temp_1(:,:,slice_idx))
    title("temp_1")
    pause(0.0001)
end


%%


temp_2 = imfill(temp_1, 18, "holes"); % 18
for slice_idx=1:nSlice
    figure(1); clf;
    imshow(temp_2(:,:,slice_idx))
    title("temp_2")
    pause(0.0001)
end


temp_3 = imerode(temp_2, strel('diamond', 10));
for slice_idx=1:nSlice
    figure(1); clf;
    imshow(temp_3(:,:,slice_idx))
    title("temp_3")
    pause(0.0001)
end

mask = imfill(temp_3, 26, "holes");


for slice_idx=1:nSlice
    figure(1); clf;
    imshow(mask(:,:,slice_idx))
    title("mask")
    pause(0.0001)
end



for slice_idx=1:nSlice
    
    figure(1); clf;
    subplot(1, 2, 1);
    imshow(doubleStretchedSlice(:,:,slice_idx))
    title(['Label Slice', num2str(slice_idx)]);
    
    se = strel('disk', 7); % Prova con diversi raggi
    mask = imopen(mask, se);
    subplot(1, 2, 2);
    imshow(mask(:,:,slice_idx),[])
    title(['Label Slice', num2str(slice_idx)]);
    pause(0.0001);
end
%%
% ------------ START OF TEST CODE ----------- %
%slice_idx = 154
mod = zeros(size(trainVolume,1), size(trainVolume,2), nSlice, 'uint16');

for slice_idx=1:nSlice

    figure(2); clf;
    tempSlice = doubleStretchedSlice(:,:,slice_idx);
    
    [~,threshold] = edge(tempSlice,'sobel');
    fudgeFactor = 0.5;
    BWs = edge(tempSlice,'sobel',threshold * fudgeFactor);

    subplot(1, 2, 1);
    imshow(tempSlice);
    
    se90 = strel('line',3,90);
    se0 = strel('line',3,0);
    
    BWsdil = imdilate(BWs,[se90 se0]);
    %imshow(BWsdil)
    %title('Dilated Gradient Mask')
    
    BWdfill = imfill(BWsdil,'holes');
    %imshow(BWdfill)
    %title('Binary Image with Filled Holes')
    
    
    BWnobord = imclearborder(BWdfill,4);
    %imshow(BWnobord)
    %title('Cleared Border Image')
    
    seD = strel('diamond',1);
    BWfinal = imerode(BWnobord,seD);
    BWfinal = imerode(BWfinal,seD);
    %imshow(BWfinal)
    %title('Segmented Image');
    
    % seD5 = strel('diamond',8);
    % BWfinal = imopen(BWfinal, seD5);
    %imshow(BWfinal)
    %title("AntoMetod")

    subplot(1, 2, 2);
    imshow(BWfinal)
    title(['Label Slice', num2str(slice_idx)]);

    pause(0.01);
    %mod(:,:,slice_idx) = BWfinal;

end

%visualizeSlices(mod, labelVolume)
%%

figure;
% Visualizzazione della slice modificata
imshow(doubleStretchedSlice(:,:,slice_idx));
se_open = strel("disk", 8);

figure;
imshow(imopen(imfill(imfilter(doubleStretchedSlice(:,:,slice_idx), fspecial("disk", 3)) > nBins / 6, "holes"), se_open))

if visualizeSlicesFlag
    visualizeSlices(normalizedSlice, stretchedSlice);
end


disp("Code End")

%[appendix]
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":56.1}
%---
