clear; clc; clf;
%%

disp("Code Run")

% Option to visualize slices
visualizeSlicesFlag = false;


% Load the slices from files
[trainVolume, labelVolume] = loadNiiFile(...
    'data/imagesTr/liver_80.nii.gz', ...
    'data/labelsTr/liver_80.nii.gz'  ...
);

disp(size(trainVolume))
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
    bones_mask(:,:,slice_idx) = double(mask);
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
DoS = 0.04; % Regola questo (Degree of Smoothing)
sSigma = 4;

h = waitbar(0);
for slice_idx=1:nSlice
     waitbar(slice_idx/nSlice, h);
    % tempSlice = uint16(double(doubleStretchedSlice_final(:,:,slice_idx)) .* double(1 - mask_array(:,:,slice_idx)));
    % doubleStretchedSlice(:,:,slice_idx) = kuwahara(tempSlice, 19);
    doubleStretchedSlice(:,:,slice_idx) = imbilatfilt(double(double(doubleStretchedSlice_final(:,:,slice_idx))./maxValue), DoS, sSigma);
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

slice_idx = 156;
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
lowerThreshold = 0.06 %[control:slider:431b]{"position":[18,22]}
upperThreshold = 0.28 %[control:slider:32e1]{"position":[18,22]}

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
% --- INPUT INIZIALE: refinedLiverMask_3D ---
[nRows, nCols, nSlice] = size(refinedLiverMask_3D); % Ottieni dimensioni


% --- FASE 1: Applica criterio componente 2D più bassa e controllo metà inferiore ---
disp('FASE 1: Applicazione criterio componente 2D più bassa e controllo metà inferiore...');
maskAfterStage1 = false(size(refinedLiverMask_3D)); % Maschera INTERMEDIA dopo la fase 1

% Usa il tuo valore per definire la 'metà' inferiore (1.5 non è metà, ma rispetto la tua scelta)
startRowBottomHalf = floor(nRows / 1.5) + 1; 

for slice_idx = 1:nSlice
    currentSliceMask = refinedLiverMask_3D(:,:,slice_idx);

    % Salta se l'intera slice è vuota
    if ~any(currentSliceMask(:))
        continue; 
    end

    % Verifica se c'è almeno un pixel nella parte 'inferiore' definita
    bottomHalfMaskPortion = currentSliceMask(startRowBottomHalf:end, :); 
    if ~any(bottomHalfMaskPortion(:))
        continue; 
    end
    
    % Analisi componenti 2D
    CC_2D = bwconncomp(currentSliceMask, 8);

    % Se 0 o 1 componente, usa la maschera intera (soddisfa già i criteri)
    if CC_2D.NumObjects <= 1 
        maskAfterStage1(:,:,slice_idx) = currentSliceMask;
    else
        % Più componenti: trova quella che si estende più in basso
        maxRowOverall = -1;
        idxLowestComponent = -1;

        for k = 1:CC_2D.NumObjects
            [rows, ~] = ind2sub([nRows, nCols], CC_2D.PixelIdxList{k});
            if isempty(rows), continue; end % Salta componenti vuote (improbabile ma sicuro)
            currentComponentMaxRow = max(rows); 

            if currentComponentMaxRow > maxRowOverall
                maxRowOverall = currentComponentMaxRow;
                idxLowestComponent = k;
            end
        end
        
        % Crea maschera temporanea solo con la componente più bassa trovata
        tempSliceMask = false(nRows, nCols);
        if idxLowestComponent > 0 
             tempSliceMask(CC_2D.PixelIdxList{idxLowestComponent}) = true;
        end
        maskAfterStage1(:,:,slice_idx) = tempSliceMask; % Salva il risultato nella maschera intermedia
    end
end
disp('FASE 1: Completata. Maschera intermedia creata.');
% --- FINE FASE 1 ---


% --- FASE 2: Seleziona la Componente 3D Più Grande dalla Maschera Intermedia ---
disp('FASE 2: Ricerca della componente 3D più grande nella maschera intermedia (risultato Fase 1)...');

finalProcessedMask = false(size(maskAfterStage1)); % Inizializza maschera finale (output)

% Verifica se la maschera intermedia (output della Fase 1) contiene qualcosa
if ~any(maskAfterStage1(:))
    disp('La maschera intermedia (dopo Fase 1) è vuota. La maschera finale sarà vuota.');
else
    % Trova le componenti connesse 3D nella maschera risultante dalla Fase 1
    connectivity = 26; % Connettività 3D
    CC_3D_Stage2 = bwconncomp(maskAfterStage1, connectivity);
    
    % Controlla quante componenti sono rimaste dopo la Fase 1
    if CC_3D_Stage2.NumObjects == 0
        disp('Nessuna componente 3D trovata nella maschera intermedia. Maschera finale vuota.');
    elseif CC_3D_Stage2.NumObjects == 1
         disp('Trovata una sola componente 3D dopo la Fase 1. Verrà usata quella.');
         finalProcessedMask = maskAfterStage1; % È già la più grande
    else
        % Più componenti 3D rimaste: trova la più grande per volume
        disp(['Trovate ', num2str(CC_3D_Stage2.NumObjects), ' componenti 3D (dopo Fase 1). Ricerca della più grande...']);
        numVoxels = cellfun(@numel, CC_3D_Stage2.PixelIdxList);
        [maxSize, idxLargestComponent] = max(numVoxels); 
        fprintf('La componente più grande (dopo Fase 1) ha %d voxel (indice %d).\n', maxSize, idxLargestComponent);
        
        % Crea la maschera finale contenente solo la componente più grande trovata
        finalProcessedMask(CC_3D_Stage2.PixelIdxList{idxLargestComponent}) = true; 
    end
end

disp('FASE 2: Selezione della componente 3D più grande completata.');
% --- FINE FASE 2 ---


% --- Visualizzazione (usa la finalProcessedMask risultante dalla Fase 2) ---
disp('Visualizzazione maschera finale processata (dopo Fase 1 + Fase 2)...');
for slice_to_show=1:nSlice
    figure(121); clf; 
    subplot(1,2,1);
     if slice_to_show <= size(doubleStretchedSlice, 3)
        imshow(doubleStretchedSlice(:,:,slice_to_show), []);
        title(['Slice Filtrata ', num2str(slice_to_show)]);
     else
         title(['Slice Filtrata ', num2str(slice_to_show), ' (Dati non disp.)']);
     end

    subplot(1,2,2);
    if slice_to_show <= size(doubleStretchedSlice, 3)
        imshow(doubleStretchedSlice(:,:,slice_to_show), []);
    else
         imagesc(zeros(nRows, nCols)); 
         axis image; colormap gray;
    end
    hold on;
    % Usa la finalProcessedMask (output della Fase 2) per la visualizzazione
    visboundaries(finalProcessedMask(:,:,slice_to_show),'Color','m', 'LineWidth', 1); % Colore Magenta
    hold off;
    title(['Finale (LowestComp & Largest3D) (Slice ', num2str(slice_to_show), ')']);

    pause(0.01); 
end


%%
image_class = class(no_bones_slice); 
mask_casted = cast(finalProcessedMask, image_class);


sliceFinalMask = no_bones_slice .* mask_casted;

for i=1:nSlice
    figure(404); clf;
    imshow(sliceFinalMask(:,:,i))

    pause(0.001);
end

%%
[hMean_final, hMean_clean_final] = histogramOnAllSlices(sliceFinalMask, maxValue);

[grouped_hMean_final, grouped_hMean_clean_final, binCenters] = groupHistogramData(hMean_final, hMean_clean_final, groupSize, nBins);


[lowerIntensity_final, upperIntensity_final] = bandDetection(grouped_hMean_clean_final, 3000)

plotGroupedHistograms(binCenters, grouped_hMean_final, grouped_hMean_clean_final, nBins, 0.42, 0.55);


stetched_final = stretchSlices(no_bones_slice, 0.42, 0.55, 2);

%%

matchSlices = histogramMachingAllSlice(stetched_final, 150) .* mask_casted;

% gridSlices(stetched_final, 134,142)

%%

for i=1:nSlice
    figure(404); clf;
    subplot(1,2,1);
    imshow(matchSlices(:,:,i),[])
    title(i)
    subplot(1,2,2);
    imshow(labelVolume(:,:,i),[])
    pause(0.001);
end


%%
DoS = 0.04; % Regola questo (Degree of Smoothing)
sSigma = 7;

h = waitbar(0);
for slice_idx=1:nSlice
     waitbar(slice_idx/nSlice, h);
    % tempSlice = uint16(double(doubleStretchedSlice_final(:,:,slice_idx)) .* double(1 - mask_array(:,:,slice_idx)));
    % doubleStretchedSlice(:,:,slice_idx) = kuwahara(tempSlice, 19);
    matchSlices_filter(:,:,slice_idx) = imbilatfilt(double(double(matchSlices(:,:,slice_idx))./maxValue), DoS, sSigma);
end
%%
for i=1:nSlice
    figure(404); clf;
    subplot(1,2,1);
    imshow(matchSlices_filter(:,:,i),[])
    title(i)
    subplot(1,2,2);
    imshow(labelVolume(:,:,i),[])
    pause(0.001);
end


%%
slice_idx = 157;
slice = matchSlices_filter(:,:,slice_idx);
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
disp('--- Inizio Segmentazione Tumori ---');

% --- Input per la segmentazione dei tumori ---
liverIntensityVolume = matchSlices_filter;     % Volume con intensità solo del fegato (già mascherato)
liverMaskVolume = finalProcessedMask; % Maschera binaria 3D del fegato

% --- Parametri per la segmentazione dei tumori (DA REGOLARE!) ---

% 1. Soglia di Intensità:
%    Identifica i pixel candidati ad essere tumore (più scuri).
%    Questo è il parametro PIÙ CRITICO. Ispeziona i valori in matchSlices(:,:,157)
%    nelle aree tumorali (cerchiate) e nel fegato sano circostante.
%    Scegli un valore che sia SOTTO l'intensità del fegato sano ma SOPRA
%    l'intensità dei tumori. Potrebbe richiedere prove.
%    ESEMPIO: Se i tumori hanno intensità ~5000 e il fegato ~15000,
%    una soglia potrebbe essere intorno a 8000-10000.
%    ATTENZIONE: I valori dipendono molto dal preprocessing fatto prima!
%    Guarda l'istogramma di liverIntensityVolume(liverMaskVolume) per aiutarti.
tumorIntensityThreshold_lower = 0.045; % VALORE PURAMENTE INDICATIVO - DA CAMBIARE!
tumorIntensityThreshold_upper = 0.21;

% 2. Dimensione Minima Oggetti (Pulizia):
%    Rimuove piccole regioni (rumore, vasi) dopo il thresholding.
%    Puoi specificare un'area minima per slice (2D) o un volume minimo (3D).
%    3D è generalmente più robusto.
minTumorVolumeVoxels = 5000; % Volume minimo in voxel per considerare una regione come tumore -> DA REGOLARE!

% 3. Elementi Strutturanti Morfologici (Opzionale, per pulizia):
se_open_tumor = strel('disk', 3); % Per rimuovere piccole connessioni/rumore (imopen)
se_close_tumor = strel('disk', 3);% Per chiudere piccoli buchi nei tumori (imclose)
% Per operazioni 3D, potresti usare: strel('sphere', 1) o strel('sphere', 2)

% --- Fase 1: Thresholding Iniziale ---
% Crea una maschera binaria iniziale: pixel nel fegato SOTTO la soglia
% Assicurati che liverIntensityVolume non sia vuoto o contenga solo zeri
initialTumorMask = false(size(liverIntensityVolume));
validLiverPixels = liverMaskVolume & (liverIntensityVolume > 0); % Pixel del fegato con intensità > 0
initialTumorMask(validLiverPixels) = (liverIntensityVolume(validLiverPixels) >= tumorIntensityThreshold_lower & liverIntensityVolume(validLiverPixels) <= tumorIntensityThreshold_upper);

disp('Thresholding iniziale per tumori completato.');

% --- Fase 2: Raffinamento Morfologico e Selezione Componenti (3D) ---
disp('Inizio raffinamento morfologico 3D e rimozione piccoli oggetti...');
tic;

% 1. (Opzionale) Apertura Morfologica 3D: Rimuove piccole protuberanze/rumore
%    Potrebbe rimuovere anche tumori molto piccoli se se_open_tumor è grande
openedMask_3D = imopen(initialTumorMask, strel('sphere', 5)); % Usa strel 3D
% Potresti preferire applicare l'apertura 2D slice-by-slice se la 3D è troppo aggressiva
% openedMask_3D = false(size(initialTumorMask));

currentMask = openedMask_3D; % Maschera da usare per i passi successivi

% 2. Rimuovi Oggetti Piccoli (Basato sul Volume 3D)
%    Questa è spesso la pulizia più efficace. Usa bwareaopen in modalità 3D.
%    La connettività 18 o 26 è comune per il 3D.
connectivity_3d = 18; % o 26
cleanedTumorMask_3D = bwareaopen(currentMask, minTumorVolumeVoxels, connectivity_3d);

% 3. (Opzionale) Chiusura Morfologica 3D: Riempie piccoli buchi interni
% closedMask_3D = imclose(cleanedTumorMask_3D, se_3D_tumor); % Usa strel 3D
% Oppure slice-by-slice:
closedMask_3D = false(size(cleanedTumorMask_3D));
for k=1:size(cleanedTumorMask_3D, 3)
     closedMask_3D(:,:,k) = imclose(cleanedTumorMask_3D(:,:,k), se_close_tumor);
     closedMask_3D(:,:,k) = imfill( closedMask_3D(:,:,k),"holes");
end
finalTumorMask = closedMask_3D; % Maschera finale dopo la pulizia

% 4. (Opzionale Extra) Riempimento Buchi Completo (se necessario)
finalTumorMask = imfill(finalTumorMask, 'holes'); % Funziona in 3D

elapsedTimeTumor = toc;
fprintf('Raffinamento morfologico 3D completato in %.2f secondi.\n', elapsedTimeTumor);

% --- Visualizzazione Finale Tumori ---
disp('Visualizzazione segmentazione finale tumori...');
for slice_idx = 120:size(liverIntensityVolume, 3)
    figure(500); clf; % Usa un nuovo numero di figura

    subplot(1,2,1);
    imshow(labelVolume(:,:,slice_idx), []);
    title(['Fegato Preprocessato (Slice ', num2str(slice_idx), ')']);

    subplot(1,2,2);
    imshow(liverIntensityVolume(:,:,slice_idx), []);
    hold on;
    % Mostra contorno fegato (verde) e tumori (rosso)
    visboundaries(liverMaskVolume(:,:,slice_idx), 'Color', 'g', 'LineWidth', 0.5, 'LineStyle', ':'); % Contorno fegato
    visboundaries(finalTumorMask(:,:,slice_idx), 'Color', 'r', 'LineWidth', 1); % Contorno tumori
    hold off;
    title(['Fegato (Verde) e Tumori Rilevati (Rosso) - Slice ', num2str(slice_idx)]);

    drawnow;
    pause(0.01); % Pausa per visualizzazione
end

disp('--- Fine Segmentazione Tumori ---');

% La variabile 'finalTumorMask' ora contiene la maschera 3D binaria
% delle masse tumorali identificate.

%%
% --- INIZIO CODICE VALUTAZIONE ---
% Qui finisce il tuo codice di segmentazione. Le variabili importanti sono:
% labelVolume: Ground Truth (0, 1, 2)
% liverMaskVolume: Predizione maschera Fegato (logica 0/1)
% finalTumorMask: Predizione maschera Tumore (logica 0/1)

disp('--- Inizio Calcolo Metriche di Valutazione ---');

% --- Preparazione Maschere Predette (assicurati siano logiche) ---
liver_mask_pred = logical(liverMaskVolume);
tumor_mask_pred = logical(finalTumorMask);

% --- Preparazione Maschere Ground Truth ---
tumor_mask_gt = (labelVolume == 2);
liver_mask_gt = (labelVolume == 1) | tumor_mask_gt; % SOLO fegato GT

% --- 1. Calcolo Recall (Sensitivity) ---
%    (Questa parte usa le maschere predette separate ed è corretta)

fprintf('\n--- Recall (Sensitivity) ---\n');
% Valutazione Fegato (solo label 1)
tp_liver = sum(liver_mask_pred(:) & liver_mask_gt(:));
total_gt_liver = sum(liver_mask_gt(:));
if total_gt_liver > 0
    liver_recall = tp_liver / total_gt_liver;
else
    liver_recall = NaN; % O 1.0 se tp_liver è 0? Dipende da interpretazione
end
fprintf('Recall Fegato (TP / GT Fegato): %.4f%%  (TP=%d, Totale GT=%d)\n', liver_recall *100, tp_liver, total_gt_liver);

% Valutazione Tumore (label 2)
tp_tumor = sum(tumor_mask_pred(:) & tumor_mask_gt(:));
total_gt_tumor = sum(tumor_mask_gt(:));
if total_gt_tumor > 0
    tumor_recall = tp_tumor / total_gt_tumor;
else
    tumor_recall = NaN; % O 1.0 se tp_tumor è 0?
end
fprintf('Recall Tumore (TP / GT Tumore): %.4f%%  (TP=%d, Totale GT=%d)\n', tumor_recall*100, tp_tumor, total_gt_tumor);


% --- 2. Calcolo Dice Score (Metrica Standard) ---
%    (Questa parte richiede la costruzione di predictionVolume con etichette 0, 1, 2)

fprintf('\n--- Dice Score ---\n');
% COSTRUZIONE di predictionVolume dalle maschere predette:
predictionVolume = zeros(size(labelVolume), 'like', labelVolume); % Inizializza a 0 (background)
predictionVolume(liver_mask_pred) = 1; % Imposta fegato predetto a 1
predictionVolume(tumor_mask_pred) = 2; % Imposta tumore predetto a 2 (sovrascrive fegato se necessario)
disp('Volume di Predizione (0,1,2) costruito.');

% Ora puoi chiamare la funzione per il Dice Score
% Assicurati che la funzione sia nel path di MATLAB o nella stessa cartella
try
    scores = calculateDiceScoresMATLAB(labelVolume, predictionVolume, [1, 2]); % Passa GT e Predizione (0,1,2)
    disp('Dice Scores calcolati:');
    disp(scores);

    % Calcolo opzionale dello score medio (ignorando NaN)
    all_scores = struct2array(scores);
    valid_scores = all_scores(~isnan(all_scores)); % Rimuove NaN
    if ~isempty(valid_scores)
        mean_dice = mean(valid_scores);
        fprintf('Mean Dice Score (su classi valide): %.4f\n', mean_dice);
    else
        disp('Nessuno score Dice valido calcolato per la media.');
    end

catch ME
    warning('Errore durante il calcolo del Dice Score. Assicurati che la funzione "calculateDiceScoresMATLAB" sia disponibile.');
    disp(ME.message);
    scores = struct(); % Crea una struct vuota per evitare errori dopo
end

disp('--- Fine Calcolo Metriche di Valutazione ---');


% --- Parte Grafica (già presente nel tuo codice) ---
% La visualizzazione che hai implementato alla fine della segmentazione
% mostra già l'overlay della maschera predetta (finalTumorMask in rosso)
% sull'immagine del fegato preprocessato, che è un ottimo modo per
% vedere graficamente il risultato della *tua segmentazione*.

% Se vuoi confrontare graficamente la TUA predizione (rossa) con il
% GROUND TRUTH (blu), puoi modificare quella visualizzazione:

disp('Visualizzazione comparativa GT vs Predizione (slice per slice)...');
[row, col, nSlice] = size(labelVolume); % Prendi nSlice da labelVolume
for slice_idx = 1:nSlice % Usa nSlice corretto
    figure(600); clf; % Usa un altro numero di figura per non sovrascrivere

    imshow(liverIntensityVolume(:,:,slice_idx), []); % Mostra immagine base
    hold on;

    % Contorni Ground Truth (es. blu e ciano)
    visboundaries(tumor_mask_gt(:,:,slice_idx), 'Color', 'b', 'LineWidth', 1, 'LineStyle', '-'); % GT Tumore
    visboundaries(liver_mask_gt(:,:,slice_idx), 'Color', 'c', 'LineWidth', 0.5, 'LineStyle', '-'); % GT Fegato (solo label 1)

    % Contorni Predizione (es. rosso e giallo)
    visboundaries(tumor_mask_pred(:,:,slice_idx), 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--'); % Predizione Tumore
    visboundaries(liver_mask_pred(:,:,slice_idx), 'Color', 'y', 'LineWidth', 0.5, 'LineStyle', '--'); % Predizione Fegato

    hold off;
    title(sprintf('Slice %d: GT (blu/ciano) vs Pred (rosso/giallo)', slice_idx));

    drawnow;
    pause(0.05); % Pausa leggermente più lunga per vedere
end


%%
disp('--- Inizio Ricostruzione 3D (isosurface) ---');

bones_mask_2 = imopen(bones_mask, strel("sphere", 2));

% --- Parametri ---
liverColor = [0 1 0]; % Verde per il fegato
tumorColor = [1 0 0]; % Rosso per il tumore
bonesColor = [1 1 0.94];


liverAlpha = 0.4;     % Trasparenza per il fegato (per vedere il tumore dentro)
tumorAlpha = 1.0;     % Opaco per il tumore
bonesAlpha = 1.0;
smoothing = true;     % Applica smoothing per superfici meno "a blocchi"?
smoothFactor = 5;     % Fattore di smoothing (se attivo)

% Assicurati che le maschere siano logiche
liver_mask_pred = logical(liverMaskVolume);
tumor_mask_pred = logical(finalTumorMask);

% --- Crea la Figura ---
figure;
hold on; % Mantiene gli oggetti nella stessa figura

% --- Ricostruzione Fegato ---
fprintf('Ricostruzione superficie fegato...\n');
% Converte in double per isosurface e smoothing
mask_double_liver = double(liver_mask_pred);
if smoothing
    mask_double_liver = smooth3(mask_double_liver, 'gaussian', smoothFactor);
end
% Calcola la superficie (al livello 0.5 per maschere binarie)
[faces_liver, verts_liver] = isosurface(mask_double_liver, 0.5);
% Disegna la patch (superficie)
h_liver = patch('Faces', faces_liver, 'Vertices', verts_liver);
% Imposta proprietà grafiche
set(h_liver, 'FaceColor', liverColor, 'EdgeColor', 'none', 'FaceAlpha', liverAlpha);
fprintf('Fegato disegnato.\n');

% --- Ricostruzione Tumore ---
% (Solo se ci sono voxel tumorali)
if any(tumor_mask_pred(:))
    fprintf('Ricostruzione superficie tumore...\n');
    mask_double_tumor = double(tumor_mask_pred);
     if smoothing
        mask_double_tumor = smooth3(mask_double_tumor, 'gaussian', smoothFactor);
    end
    [faces_tumor, verts_tumor] = isosurface(mask_double_tumor, 0.5);
    h_tumor = patch('Faces', faces_tumor, 'Vertices', verts_tumor);
    set(h_tumor, 'FaceColor', tumorColor, 'EdgeColor', 'none', 'FaceAlpha', tumorAlpha);
    fprintf('Tumore disegnato.\n');
else
    fprintf('Nessun tumore trovato nella maschera predetta, skip ricostruzione tumore.\n');
end


% --- Ricostruzione Ossa ---
% (Solo se ci sono voxel tumorali)
if any(bones_mask(:))
    fprintf('Ricostruzione superficie tumore...\n');
    mask_double_bones = double(bones_mask_2);
     if smoothing
        mask_double_bones = smooth3(mask_double_bones, 'gaussian', smoothFactor);
    end
    [faces_bones, verts_bones] = isosurface(mask_double_bones, 0.5);
    h_bones = patch('Faces', faces_bones, 'Vertices', verts_bones);
    set(h_bones, 'FaceColor', bonesColor, 'EdgeColor', 'none', 'FaceAlpha', bonesAlpha);
    fprintf('Ossa disegnate.\n');
else
    fprintf('Nessun tumore trovato nella maschera predetta, skip ricostruzione tumore.\n');
end

% --- Impostazioni Finali Vista 3D ---
hold off;
title('Ricostruzione 3D - Fegato (Verde Trasparente) e Tumore (Rosso Opaco)');
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;        % Mostra griglia
axis equal;     % Assicura che le proporzioni siano corrette
view(3);        % Imposta la vista 3D standard
camlight;       % Aggiunge una luce per migliorare la visualizzazione
lighting gouraud; % Tipo di illuminazione
rotate3d on;    % Permette la rotazione interattiva con il mouse

disp('--- Fine Ricostruzione 3D (isosurface) ---');

%%
for i = 1:207
    figure(123);clf;
    imshow(bones_mask(:,:,i),[]);
    pause(0.01);
end

%[appendix]
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":43.7}
%---
%[control:slider:431b]
%   data: {"defaultValue":0,"label":"lowerThreshold","max":0.2,"min":0,"run":"Section","runOn":"ValueChanging","step":0.01}
%---
%[control:slider:32e1]
%   data: {"defaultValue":0,"label":"upperThreshold","max":0.3,"min":0,"run":"Section","runOn":"ValueChanging","step":0.01}
%---
