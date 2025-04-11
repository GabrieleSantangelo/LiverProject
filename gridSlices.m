function [] = gridSlices(slices, startSlice, endSlice)

num_slices_to_display = endSlice - startSlice; 


grid_cols = ceil(sqrt(num_slices_to_display));
grid_rows = ceil(num_slices_to_display / grid_cols);

figure; 

for i = 1:num_slices_to_display
    % Indice della slice corrente nella matrice 3D originale
    slice_index = startSlice + i - 1;

    % Seleziona la slice 2D
    currentSlice = slices(:, :, slice_index);

    % Attiva l'asse corretto nella griglia (subplot)
    % subplot(numero_righe, numero_colonne, indice_corrente)
    subplot(grid_rows, grid_cols, i);

    imshow(currentSlice,[]);

    % Aggiungi un titolo a ciascuna subplot (opzionale)
    title(['Slice ', num2str(slice_index)]);

    % Rimuovi i segni di graduazione dagli assi per un look più pulito (opzionale)
     set(gca, 'XTick', [], 'YTick', []);
    % oppure usa axis off; (ma nasconde anche il riquadro)

end

% Aggiungi un titolo generale all'intera figura (richiede MATLAB più recente)
if exist('sgtitle', 'file')
    sgtitle(['Visualizzazione Griglia Slices da ', num2str(startSlice), ' a ', num2str(startSlice + num_slices_to_display - 1)]);
end
end