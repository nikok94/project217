function ResourceElementGrids = bwpToResourceElementGrids(waveconfig)
  %  waveconfig включает в себя следующие структуры :
  %  
  %  ssburst = waveconfig.SSBurst;
    carriers = waveconfig.Carriers;
    bwp = waveconfig.BWP;
    coreset = waveconfig.CORESET;
    pdcch = waveconfig.PDCCH;
    pdsch = waveconfig.PDSCH;
    
    % Находим максимальное количесво layers в PSCH
    maxlayers = -1;
    for i = 1 : length(pdsch)
        if pdsch(i).Enable == 1 && ~isempty(pdsch(i).NLayers)
            if pdsch(i).NLayers > maxlayers
                maxlayers = pdsch(i).NLayers;
            end
        end
    end
    
    % Создаем массив из length(bwp) элементов, каждый элемент это трехмерный 
    % массив M*S*L - сетка ресурс элементов для одного из указанных BWP
    % M = NRB*12, NRB текущего BWP
    % S - количество символов во фрейме для текущего BWP
    % L - максимальное количество layer-ов из всех BWP
    carrierscs = [carriers.SubcarrierSpacing]; 
    ResourceElementGrids = cell(1, length(bwp));
    for i = 1 : length(bwp)
        t = bwp(i).SubcarrierSpacing == carrierscs;
        cidx = find(t,1);
        if strcmpi(bwp(i).CyclicPrefix,'Normal')
            symbolsperslot = 14;
        else
            symbolsperslot = 12;
        end
        ResourceElementGrids{i} = zeros(bwp(i).NRB*12,waveconfig.NumSubframes*1*symbolsperslot*fix(bwp(i).SubcarrierSpacing/15),maxlayers);
        bwp(i).CarrierIdx = cidx;
    end 
    
    for nch = 1 : length(pdcch)
        ch = pdcch(nch);
        if ~ch.Enable
            continue;
        end
        
        if ch.BWP < 1 || ch.BWP > length(bwp)
            error('bwpToResourceElementGrids: cch BWP index error');
        end
        
        if ch.CORESET < 1 || ch.CORESET > length(coreset)
            error('bwpToResourceElementGrids: cch CORESET index error');
        end
        
        if strcmpi(bwp(ch.BWP).CyclicPrefix,'Normal')
            symbperslot = 14;
        else
            symbperslot = 12;
        end
        % 
        %%
        % копируем нужный CORESET
        cset = coreset(pdcch(nch).CORESET);
        % проверяем какие выделенные для CORESET символы выходят за пределы
        % одного слота
        slotsymbs = cset.AllocatedSymbols(cset.AllocatedSymbols+cset.Duration < symbperslot);
        % если выделеные символы выходят за пределы одного слота то эти
        % символы не используются при дальнейшей конфигурации сетки
        % ресурсов
        if length(slotsymbs) ~= length(cset.AllocatedSymbols)
            warning('CORESET %d (%d symbol duration) in BWP %d includes positions which fall outside the slot (0...%d). Using only CORESET locations within a slot.',ch.CORESET,cset.Duration,ch.BWP,symbperslot-1);
            cset.AllocatedSymbols = slotsymbs;
        end
        
        alocatedSlots = expandbyperiod(cset.AllocatedSlots,cset.AllocatedPeriod,waveconfig.NumSubframes, bwp(nch).SubcarrierSpacing);
%        alocatedSlots1 = RepeatByPeriod(cset.AllocatedSlots, cset.AllocatedPeriod, waveconfig.NumSubframes*1*fix(bwp(nch).SubcarrierSpacing/15));
        potentialsymbols = reshape(symbperslot*alocatedSlots + cset.AllocatedSymbols',1,[]);
        
      %  allocslotindices = RepeatByPeriod(ch.AllocatedSearchSpaces,ch.AllocatedPeriod,numel(potentialsymbols));
        allocslotindices = expandbyperiod(ch.AllocatedSearchSpaces,ch.AllocatedPeriod,numel(potentialsymbols));
        
        
        % Находим все начальные интервалы символов задействованных CORESET
        % для значений PDCCH
        allocatedsymbols = potentialsymbols(1+allocslotindices);
        
        if max(cset.AllocatedPRB) >= bwp(ch.BWP).NRB
            error('The PRB allocation part of CORESET %d (set of 6 PRB/REG CCEs, max 0-based PRB index = %d) exceeds the number of available RB (%d) in BWP %d.',ch.CORESET,max(cs.AllocatedPRB),bwp(ch.BWP).NRB,ch.BWP);
        end
        % выравнивание индексов CCE так, чтобы было кратно 6, проверка
        % если нет необходимости, то проверку можно пропустить  
        ncces = unique(fix(cset.AllocatedPRB/6)); 
        ncces = ncces(ncces < fix(bwp(ch.BWP).NRB/6));
        % получаем индексы PRB для всех ССE в символе 
        potentialPRB = expander(ncces*6,6);
        
        % все индексы REG для всех возможных CCE
        if ch.CCEInterleaver
            if ch.CCEREGbundleSize < cset.Duration
                error('CCEREGbundleSize invalid : ch.CCEREGbundleSize < cset.Duration');
            end
        end
        [REGindexForCCE, f] = bfPDCCHinterleaver(length(cset.AllocatedPRB)*6*cset.Duration, ch.CCEREGbundleSize, ch.CCEInterleaverSize, ch.CCEnshift, ch.CCEInterleaver);

        CCH_CCE_index = zeros(ch.NumCCE*6,1);
        % Выбираем индексы для нужного количества CCE
        m = 0;
        for k = ch.StartCCE : ch.StartCCE + ch.NumCCE - 1
            CCH_CCE_index(m + 1 : m + 6, 1) = REGindexForCCE(:, k + 1);
            m = m + 6;
        end
        CCH_CCE_index = CCH_CCE_index';
        
        prblocations = [mod(CCH_CCE_index,cset.Duration)',potentialPRB(1+fix(CCH_CCE_index/cset.Duration))'];
        celloc1d = cell(1,cset.Duration);
        for c = 1:cset.Duration
          celloc1d{c} = reshape(prblocations(prblocations(:,1)==(c-1),2),1,[]);  % Store each PRB set as a row
        end  
        
        indpdcch = reSub2Ind(bwp(ch.BWP).NRB,celloc1d,0:cset.Duration-1);
        repdcch = expander(12*indpdcch,12,4,1,true)';  % Step by 4 RE with an offset of 1     
        redmrs = expander(12*indpdcch,12,4,1,false)';  % RE carrying DM-RS
        % Заполняем символы значениями PDCCH
        
        datasource = hVectorDataSource(ch.DataSource); 
        
        nID = ch.NID;
        Gd = length(repdcch);
        G = Gd*2;
        
        for s = allocatedsymbols
            
            if ~isfield(ch,'EnableCoding') || ch.EnableCoding
                % Генерируем случайную последовательность бит, вместо DCI
                dcibits = datasource.getPacket(ch.DataBlkSize);

                % Кодируем 
                codeword = bfDCIEncode(dcibits,ch.RNTI,G);
            else
                % получем пакет бит 
                codeword = datasource.getPacket(G);
            end
            
            % Скремблинг 
            % На вход данной функции поступает битовая последовательность длинной E, которая является
            % выходом процедуры "Согласование скоростей", значение N_ID (лежит в диапазоне значений {0,1,2,...,65535}, если 
            % задан параметр верхнего уровня pdcch-DMRS-ScramblingID и мы работаем в специальном пространстве поиска, 
            % либо равен N_ID_cell) и n_RNTI(лежит в диапазоне значений {0,1,2,...,65535}, если 
            % задан параметр верхнего уровня pdcch-DMRS-ScramblingID и мы работаем в специальном пространстве поиска, 
            % либо равен 0). Выходом функции является последовательность длинной N,
            % которая получена в соответствии с стд.38.211 разд.7.3.2.3
    
            % Сгенерируем скремблирующую последовательность scramb_seq в
            % соответствии с стд.38.211 разд.5.2.1
            c_init = mod(double(ch.RNTI)*2^16 + double(nID),2^31);
            scramb_seq = scrambling(c_init, length(codeword));
            cSeq = scramb_seq';
    
            % Скремблинг, на выходе получаем сигнал диной E.
            codeword_scrambl = xor(codeword,cSeq);
    
            % QPSK модуляция
            % На вход поступает последовательность длиной E, выходомибудет набор комплексных
            % символов QPSK модуляции.
            symbols = bfModulator(codeword_scrambl, 'QPSK');
           % symbols1 = nrPDCCH(codeword,nID,ch.RNTI);
           %  if symbols~= symbols1
           %     error('');
           % end 
            offset = 1 + s*12*bwp(ch.BWP).NRB;
            ResourceElementGrids{ch.BWP}(offset+repdcch) = ResourceElementGrids{ch.BWP}(offset+repdcch) + symbols*db2mag(ch.Power);
            nslot = mod(fix(s/symbperslot), bwp(ch.BWP).SubcarrierSpacing/15 * 10);  % Номер слота длительностью 10мс
            nsym = mod(s,symbperslot); % Номер символа в слоте
            % Сформируем DMRS для PDCCH для каждого символа OFDM
            dmrssym = [];   
            crboffset = carriers(bwp(ch.BWP).CarrierIdx).RBStart + bwp(ch.BWP).RBOffset; % первый CRB, связанный с началом BWP;
            
            for i = 1:cset.Duration

                    prb = celloc1d{i};     % PRB индексы для текущего OFDM символа
                    slen = max(prb)+1;     % Количество PRB необходимое для формирования непрерывной последовательности PRBS 

                    % Сгенерируем инициализирующую последовательность c_init для псевдо-случайной
                    % последовательности (c) в соответствии с стд. 38.211. разд. 7.4.1.3.1
                    cinit = mod(2^17*(symbperslot*nslot+nsym+1)*(2*nID+1)+2*nID,2^31);
                    nsc = 6;   % количество SC необходимое для записи DMRS на один RB,3 символа DM-RS на RB и 2 бита PRBS, необходимых на символ QPSK 
                    % Сгенерируем псевдо-случайную последовательность c в соответствии с
                    % стд.38.211 разд. 5.2.1
                    n = (crboffset+slen) * nsc; 
                    seq = scrambling(cinit,n);
                    c_seq = seq(end - slen*nsc + 1:end);
                    % Зададим соответствие DMRS символов и PRB выделенных под CORESET в данном
                    % символе OFDM
                    k = 0;
                    for l = 0:slen - 1
                        j = 0;
                        while j <= (nsc - 1)           
                            DMRS_seq(j+1,l+1) = c_seq(k+1);
                            j = j + 1; 
                            k = k + 1;
                        end
                        if l == slen-1
                            break
                        end
                    end

                    DMRS_seq = DMRS_seq(:,prb+1); 
                    dmrssym = [dmrssym; bfModulator(DMRS_seq(:),'QPSK')];
                    nsym = nsym + 1;  % переходим на следующий символ 
            end
            ResourceElementGrids{ch.BWP}(offset+redmrs) = ResourceElementGrids{ch.BWP}(offset+redmrs) + dmrssym*db2mag(ch.Power + ch.PowerDMRS);          
        end
    end 
    
    dlsch = nrDLSCH('MultipleHARQProcesses',false);
    codeRate = dlsch.TargetCodeRate;
    for nch = 1 : length(pdsch)
        
        ch = pdsch(nch);
        
        if ~ch.Enable 
            continue;
        end

       % allocatedSlots1 = RepeatByPeriod(ch.AllocatedSlots, ch.AllocatedPeriod, waveconfig.NumSubframes*1*fix(bwp(nch).SubcarrierSpacing/15));
        allocatedSlots = expandbyperiod(ch.AllocatedSlots,ch.AllocatedPeriod, waveconfig.NumSubframes, bwp(nch).SubcarrierSpacing);
        
        if strcmpi(bwp(ch.BWP).CyclicPrefix,'Normal')
            symbperslot = 14;
        else
            symbperslot = 12;
        end
        
        nrb = bwp(ch.BWP).NRB;
        nID = ch.NID;
        ch.PRBSet = ch.AllocatedPRB;
        ch.SymbolSet = ch.AllocatedSymbols;
        crboffset = carriers(bwp(ch.BWP).CarrierIdx).RBStart + bwp(ch.BWP).RBOffset; % First CRB associated with start of BWP
        ch.PRBRefPoint = crboffset; % First CRB associated with start of BWP
        
        dmrsSymbol = [];
        dmrsIndexe = [];
        
        for i = 1 : length(allocatedSlots)

            % Get current slot number 
            s = allocatedSlots(i);
            ch.NSlot = s;

            % Create an empty slot grid to contain a single PDSCH instance
            slotgrid = ones(12*nrb,symbperslot,ch.NLayers).*-1;
            
            %slotgrid(dmrsIndexe) = dmrsSymbol;
            
            % Get the slot-oriented PDSCH indices, DM-RS indices and DM-RS symbol values  
            [pdschREindices,dmrsREindices,dmrsSymbols,modinfo] = ...
                hPDSCHResources(struct('NRB',bwp(ch.BWP).NRB,'CyclicPrefix',bwp(ch.BWP).CyclicPrefix,...
                                       'SubcarrierSpacing',bwp(ch.BWP).SubcarrierSpacing),ch);
                       
            for k = 0 : ch.NLayers - 1
            ch.NLayer = k + 1;
            [dmrsSymb, dmrsIndex, ptrsSymb, ptrsIndex] = getPDSCH_REF_values(ch, nrb, symbperslot, 0);
            dmrsSymbol(:,ch.NLayer)  =  dmrsSymb;
            dmrsIndexe(:,ch.NLayer) = dmrsIndex + 12*nrb*symbperslot*k;
            end
                                   
            if ~isfield(ch,'EnableCoding') || ch.EnableCoding
                % Get the RV value for this transmission instance
                rvidx = mod(i-1,length(ch.RVSequence))+1;
                rv = ch.RVSequence(rvidx); 

                % For the first RV in a sequence, get a new transport block from 
                % the data source and pass it to the DL-SCH processing
                if rvidx == 1
                   trblksize = hPDSCHTBS(ch,modinfo.NREPerPRB-ch.Xoh_PDSCH);
                   trblk = datasource.getPacket(trblksize);
                   setTransportBlock(dlsch,trblk);
                end

                % DL-SCH processing to create a codeword
                codeword = dlsch(ch.Modulation,ch.NLayers,modinfo.G,rv);
                codeword1 = bfDLSCH({trblk}, {ch.Modulation}, ch.NLayers, {rv}, {modinfo.G}, {codeRate});
                if length(codeword1) == 1
                   codeword1 = codeword1{1};
                end
                if codeword ~= codeword1
                    error('codeword ~= codeword1');
                end
                
                
            else
                % If transport coding is not enabled then get the codeword
                % directly from the data source
                codeword = datasource.getPacket(modinfo.G);
            end
            
            % Формируем QAM для PDSCH
            symbols = nrPDSCH(codeword,ch.Modulation,ch.NLayers,nID,ch.RNTI);
            symb_PDSCH = PDSCH(codeword,ch.Modulation,ch.NLayers,nID,ch.RNTI);
            
            if symbols ~= symb_PDSCH
               error('symbols ~= symb_PDSCH');
            end

            % Write the PDSCH and DM-RS symbols in the slot grid
            slotgrid(pdschREindices) = symb_PDSCH*db2mag(ch.Power);
            slotgrid(dmrsREindices) = dmrsSymbols*db2mag(ch.Power+ch.PowerDMRS);

            % Combine PDSCH instance with the rest of the BWP grid
            ResourceElementGrids{ch.BWP}(:,s*symbperslot+(1:symbperslot),1:ch.NLayers) = ResourceElementGrids{ch.BWP}(:,s*symbperslot+(1:symbperslot),1:ch.NLayers) + slotgrid; 

        end
    end
end

function out = RepeatByPeriod(allocatedArrayInCommon, period, CommonArrayLength) 
% allocatedArrayInCommon набор значений вектор столбец
% period задает значение через которое будут повторяться значения указанные
% в allocatedArrayInCommon
% CommonArrayLength длина конечной последовательности 
% функция RepeatByPeriod выбирает из allocatedArrayInCommon все значения 
% попадающие в границы period и последовательно уладывает их N раз, где N =
% CommonArrayLength/period;

        NPeriods = ceil(CommonArrayLength/period);
        nAllocated = 0;
        
        for k = 1 : length(allocatedArrayInCommon)
            if allocatedArrayInCommon(1, k) < period
                nAllocated = nAllocated + 1;
            end
        end
        
        out = zeros(1, NPeriods * nAllocated);
        nAllocated = 0;
        for nPer = 0 : NPeriods - 1
            for k = 1 : length(allocatedArrayInCommon)
                if allocatedArrayInCommon(1, k) < period
                    nAllocated = nAllocated + 1;
                    out(1, nAllocated) = allocatedArrayInCommon(1, k) + nPer*period;
                end
            end
        end
end

function expanded = expander(d,e,s,o,excl)
    if nargin < 5
        excl = 0;
    end
    if nargin < 4
        o = 0;
    end
    if nargin < 3
        s = 1;
    end
    eseq = (o:s:e-1)';
    if excl
        eseq = setdiff((0:e-1)',eseq);
    end
    expanded = reshape(reshape(d,1,[]) + eseq,1,[]);  % Use column expansion
end

function ind = reSub2Ind(nrb,prb,symbols)

    % Use column expansion
    if iscell(prb)
        slen = min(length(prb),length(symbols));
        ind = cell2mat(cellfun(@(x,y)reshape(reshape(x,[],1)+nrb*reshape(y,1,[]),1,[]),prb(1:slen),num2cell(symbols(1:slen)),'UniformOutput',false));
    else
        ind = reshape(reshape(prb,[],1) + nrb*reshape(symbols,1,[]),1,[]);
    end
    
end

function sp = expandbyperiod(s,p,nsf,scs)

    if nargin > 3
        % Expand s by period p for ts length
        ts = nsf*1*fix(scs/15);
    else
        ts = nsf;
    end
    % Is the period is empty then the pattern doesn't repeat, so doesn't need extending
    if isempty(p)
        sp = s;
    else
        g = (0:ceil(ts/max(p,1))-1);
        
        sp = reshape(s(s<p),[],1) + p*g;
    end
    
    if ~isempty(ts)
        sp = reshape(sp(sp < ts),1,[]);            % Trim any excess
    else
        sp = ones(1,0);
    end
end

function scramb_seq = scrambling(c_init,E)
scramb_seq=zeros(1,E);
c_init=uint32(c_init);
Nc=1600;
x1=zeros(1,Nc+E);
x1(1)=1;
x2=zeros(1,Nc+E);
for i=0:30
    x2(i+1)=bitand(c_init,uint32(1));
    c_init=bitshift(c_init,-1);
end
for n=0:Nc+E-31-1
    x1(n+31+1)=mod(x1(n+3+1)+x1(n+1),2);
    x2(n+31+1)=mod(x2(n+3+1)+x2(n+2+1)+x2(n+1+1)+x2(n+1),2);
end

for n=0:E-1
    scramb_seq(n+1)=mod(x1(n+Nc+1)+x2(n+Nc+1),2);
end
end 

