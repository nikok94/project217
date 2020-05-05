function symb_PDSCH = PDSCH(codeword,modulation,N_layers,N_ID,RNTI)
    %% PDSCH второй блок генерации сигнала
    % ?анна§ функци§ включает в себ§ процессы скремблинга, модул§ции и записи
    % символов модул§ции на уровни передачи(слои). 
    % На вход вунуции идет одно или два кодовых слова PDSCH (codeword), как указано в
    % стд.38.212 разд.7.2.6., тип модул§ции (modulation) дл§ одного или
    % нескольких кодовых слов, количество уровней передачи (N_laers),
    % идентификатор скремблировани§ (N_ID), и временной идентификатор радиосети (RNTI).
    %   - codeword - может быть вектор столбцом (представл§ющим одно кодовое слово)
    % или массив §чеек из одного или двух вектор столбцов (представл§ющих одно или два кодовых слова).
    %   - modulation - может быть указана как АQPSKњ, А16QAMњ, А64QAMњ, А256QAMњ. ?сли codeword содержит два кодовых слова, этот 
    % пор§док модул§ции будет примен§тьс§ к обоим кодовым словам. ђ качестве альтернативы, строка
    % массив или §чейка массива символьных векторов могут быть использованы дл§ указани§
    % различных схем модул§ции дл§ каждого кодового слова.
    %   - N_laers - количество уровней передачи (1 ... 4 дл§ одного кодового слова,
    % 5 ... 8 дл§ двух кодовых слов).
    %   - N_ID - представл§ет собой либо идентификатор §чейки
    % NCellID (0 ... 1007) или параметр более высокого уровн§ dataScrabingIdentityPDSCH (0 ... 1023).
    %   - RNTI - принимает значени§ в диапазоне (0 ... 65535).
    
    N_codeword = 1 + (N_layers > 4);
    
        if ~iscell(codeword)
            cell_codeword = {codeword};
        else
            cell_codeword = codeword;
        end
        if ~iscell(modulation)
            modulation_cell = {modulation};
        else
            modulation_cell = modulation;
        end
    for q = 1:N_codeword
        E = length(cell_codeword{q}); 
        c_init = (double(RNTI) * 2^15) + (double(q-1) * 2^14) + double(N_ID);
        % получим скремблирующую последовательность
        scramb_seq = scrambling(c_init,E);
        % скремблинг в соответствии с стд. 38.211 разд. 7.3.1.1
        scrambl_out{q} = xor(cell_codeword{q},scramb_seq');
        % Юодул§ци§ в соответствии с стд.38.211 разд. 7.3.1.2
        %mod_symb{q} = bfModulator(double(scrambl_out{q}),modulation_cell{q}); 
        mod_symb{q} = nrSymbolModulate(double(scrambl_out{q}),modulation_cell{q}); 
    end
    
    % запись символов модул§ции на уровни передачи(слои) , стд. 38.211 раздел 7.3.1.3
    symb_PDSCH = Layer_Mapping(mod_symb,N_layers);
end