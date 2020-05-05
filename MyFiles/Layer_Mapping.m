function symb_PDSCH = Layer_Mapping(mod_symb,N_laers)
    if ~iscell(mod_symb)
        codeword = {mod_symb};
    else
        codeword = mod_symb;
    end
    % Получим количество кодовых слов
    N_codeword = numel(codeword);
    % Получим количество слоев для каждого кодового слова
    N_laers = floor((N_laers + (0:N_codeword-1))/N_codeword);
    % Получим количество символов модуляции на слой
    N_Symb_Layer = uint32(length(codeword{1})/N_laers(1));
    symb_PDSCH = zeros(N_Symb_Layer,sum(N_laers),'like',(codeword{1}));
    for idx = 1:N_laers(1)
        symb_PDSCH(:,idx) = codeword{1}(idx:N_laers(1):end);
    end
    if N_codeword == 2
        for idx2 = 1:N_laers(2)
            symb_PDSCH(:,N_laers(1) + idx2) = codeword{2}(idx2:N_laers(2):end);
        end
    end
end