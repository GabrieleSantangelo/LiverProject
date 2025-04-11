function outputSlice = gammaCorrection(inputSlice, gammaValue)
% Applica la correzione gamma a una slice di immagine (matrice 2D).
%
% Args:
%   inputSlice: La matrice 2D di input (es. uint8, uint16, double).
%   gammaValue: Il valore gamma da applicare (es. 0.5, 1.0, 2.2).
%
% Returns:
%   outputSlice: La matrice 2D con la correzione gamma applicata,
%                nello stesso tipo di dato dell'input.

    % 1. Memorizza la classe (tipo di dato) originale per la conversione finale
    originalClass = class(inputSlice);

    % 2. Converti la slice in double precision nell'intervallo [0, 1]
    %    im2double gestisce uint8, uint16, double, single, logical ecc.
    sliceDouble = im2double(inputSlice);

    % 3. Applica la correzione gamma
    %    Usa l'operatore di potenza elemento per elemento '.^'
    %    È buona norma aggiungere un piccolo epsilon per evitare log(0) o 0^gamma
    %    se gamma è negativo, anche se per la correzione gamma standard non serve.
    epsilon = 1e-10; % Valore piccolo per stabilità numerica (spesso non necessario)
    gammaCorrectedDouble = (sliceDouble + epsilon) .^ gammaValue;

    % 4. Riconverti al tipo di dato originale
    %    Le funzioni im2uint8/im2uint16 scalano automaticamente da [0, 1]
    %    all'intervallo intero corretto.
    switch originalClass
        case 'uint8'
            outputSlice = im2uint8(gammaCorrectedDouble);
        case 'uint16'
            outputSlice = im2uint16(gammaCorrectedDouble);
        case {'double', 'single'}
            % Se l'originale era già float, di solito si mantiene nell'intervallo [0, 1]
            % o si scala se necessario. Qui lo manteniamo [0,1].
             if strcmp(originalClass, 'single')
                 outputSlice = single(gammaCorrectedDouble);
             else
                 outputSlice = gammaCorrectedDouble; % Già double
             end
        otherwise
            warning('Tipo di dato di input (%s) non gestito esplicitamente per la riconversione. Restituisco double [0, 1].', originalClass);
            outputSlice = gammaCorrectedDouble;
    end
end