% config структура которая состоит из следующих полей:
% NCellID               - идентификатор ячейки физического уровня(0...1007)

% Carriers - Структурный массив конфигураций несущих, специфичных для SCS, 
%   с полями:
%   SubcarrierSpacing - конфигурация разнесения поднесущих в кГц
%       Стандартные конфигурации 15,30,60,120,240 для нормального CP 
%       и 60 для расширенного CP
%   NRB - Количество блоков ресурсов
%   RBStart - начальный индекс CRB несущей SCS относительно «точки A»

% BWP - Структурный массив конфигураций части полосы пропускания с полями:
%   SubcarrierSpacing - конфигурация разнесения поднесущих в кГц
%       Стандартные конфигурации 15,30,60,120,240 для нормального CP 
%       и 60 для расширенного CP

%   CyclicPrefix - циклический префикс («обычный», «расширенный»)
%   NRB - количество блоков ресурсов в части полосы пропускания
%   RBOffset - начальный индекс BWP в несущей SCS


function [res] = bwpToGrids(config)

ssburst = config.SSBurst;
carriers = config.Carriers;
bwp = config.BWP;
coreset = config.CORESET;
pdcch = config.PDCCH;
pdsch = config.PDSCH;

bwCP = bwp.CyclicPrefix;

maxLayers = -1;
for n = 1 : length(pdsch)
    if pdsch(n).Enable == 1
        if pdsch(n).NLayers > maxLayers
            maxLayers = pdsch(n).NLayers;
        end
    end
end

if maxLayers < 1 
    maxLayers = 1;
end 

for bp = 1 : length(bwp)
    cidx = NaN;
    carrierscs = [config.Carriers.SubcarrierSpacing];
    for i = 1 : length(carrierscs)
        if (carrierscs(i) == bwp(bp).SubcarrierSpacing) 
            cidx = i;
            continue; 
        end
    end 
    
    if isnan(cidx)
        error('A SCS specific carrier configuration for SCS = %d kHz has not been defined. This carrier definition is required for BWP %d.',bwp(bp).SubcarrierSpacing,bp); 
    end
    
   bwp(bp).CarrierIdx = cidx;
    if strcmpi(bwp(bp).CyclicPrefix,'Normal')
        symbperslot = 14;
    else
        symbperslot = 12;
    end
    ResourceElementGrids{bp} = zeros(bwp(bp).NRB*12,config.NumSubframes*1*symbperslot*fix(bwp(bp).SubcarrierSpacing/15), maxLayers);
end

% процесс размещения PDCCH на сетке ресурс элементов
for Ncch = 1 : length(pdcch)
    
    if ~pdcch(Ncch).Enable
        continue;
    end
        
    % Находим количество символов в слоте для соответствующего BWP (зависит от CP)
    if strcmpi(bwp(pdcch(Ncch).BWP).CyclicPrefix,'Normal')
        symbperslot = 14;
    else
        symbperslot = 12;
    end
    
    % копируем необходимый CORESET
    cset = coreset(pdcch(Ncch).CORESET);
    
    slotsymbs = [];
    % исключаем все CORESET выходящие за пределы слота
    for i = 1 : length(cset.AllocatedSymbols)
        if cset.AllocatedSymbols(i) + cset.Duration < symbperslot
            slotsymbs = [slotsymbs, cset.AllocatedSymbols(i)];
        end 
    end
    
    if length(slotsymbs) ~= length(cset.AllocatedSymbols)
        warning('CORESET %d (%d symbol duration) in BWP %d includes positions which fall outside the slot (0...%d). Using only CORESET locations within a slot.',ch.CORESET,cset.Duration,ch.BWP,symbperslot-1);
        cset.AllocatedSymbols = slotsymbs;
    end
    
   
    
    pSlots = getSlotIndex(cset.AllocatedSlots, cset.AllocatedPeriod, config.NumSubframes, fix(bwp(pdcch(Ncch).BWP).SubcarrierSpacing/15)*config.NumSubframes);
    
    pSymb = getSymbIndex(pSlots, symbperslot, cset.AllocatedSymbols);
   
    allocslotindices = getSlotIndex(pdcch(Ncch).AllocatedSearchSpaces, pdcch(Ncch).AllocatedPeriod, numel(pSymb), 1);
    
    allocatedsymbols = pSymb(allocslotindices + 1);
    
    nbwp = 1;
    nrb = bwp(nbwp).NRB;
            % Establish 0-based CCE indices defined by the PRB allocation
    ncces = unique(fix(cset.AllocatedPRB/6)); 
    % Remove any partial CCE from the allocation 
    ncces = ncces(ncces < fix(nrb/6));     % Use only complete CCEs that fall within NRB
    potentialPRB = expander(6*ncces,6);    % Expand CCE to associated PRB set 
    
    expcceindices = 6*pdcch(Ncch).StartCCE + (0:6*pdcch(Ncch).NumCCE-1);     
    prblocations = [mod(expcceindices,cset.Duration)',potentialPRB(1+fix(expcceindices/cset.Duration))'];
    
     % Turn this into a cell array per symbol for reservation      
        celloc1d = cell(1,cset.Duration);
        for c = 1:cset.Duration
          celloc1d{c} = reshape(prblocations(prblocations(:,1)==(c-1),2),1,[]);  % Store each PRB set as a row
        end  
        % Repeat the 1 slot allocation cell array across all transmission occasion
    
    % Turn subscripts into 0-based RE indices for a single PDCCH instance
    indpdcch = reSub2Ind(bwp(pdcch(Ncch).BWP).NRB,celloc1d,0:cset.Duration-1);
    repdcch = expander(12*indpdcch,12,4,1,true)';  % Step by 4 RE with an offset of 1     
    redmrs = expander(12*indpdcch,12,4,1,false)';  % RE carrying DM-RS

    % PDCCH instance symbol/bit capacity
    Gd = length(repdcch);
    G = Gd*2;
    
    datasource = hVectorDataSource(pdcch(Ncch).DataSource); 
    
    for s = allocatedsymbols
        if pdcch(Ncch).EnableCoding == 1 || ~isfield(pdcch(Ncch),'EnableCoding')
           dcibits = datasource.getPacket(pdcch(Ncch).DataBlkSize);
           codeword = bfDCIEncode(dcibits,pdcch(Ncch).RNTI,G);
        else
           codeword = datasource.getPacket(G);
        end
        
        symbols = nrPDCCH(codeword,pdcch(Ncch).NID,pdcch(Ncch).RNTI);
        offset = 1+s*12*bwp(pdcch(Ncch).BWP).NRB;
        ResourceElementGrids{pdcch(Ncch).BWP}(offset+repdcch) = ResourceElementGrids{pdcch(Ncch).BWP}(offset+repdcch) + symbols*db2mag(pdcch(Ncch).Power); 
        % Construct and map the PDCCH DM-RS     
        nslot = mod(fix(s/symbperslot), bwp(pdcch(Ncch).BWP).SubcarrierSpacing/15 * 10);  % Slot number, in a 10ms frame
        nsym = mod(s,symbperslot); % Symbol number in slot

        % Construct a single symbol column vector of DM-RS for the PDCCH
        dmrssym = [];   
        crboffset = carriers(bwp(pdcch(Ncch).BWP).CarrierIdx).RBStart + bwp(pdcch(Ncch).BWP).RBOffset; % First CRB associated with start of BWP
        for i = 1:cset.Duration

             prb = celloc1d{i};     % PRB indices associated with current symbol
             slen = max(prb)+1;     % Number of PRB that the continuous PRBS sequence needs to cover 

             % Construct PRBS for the transmission DM-RS, offseting the 
             % sequence to account for the origin of the BWP
             cinit = mod(2^17*(symbperslot*nslot+nsym+1)*(2*pdcch(Ncch).NID+1)+2*pdcch(Ncch).NID,2^31);
             nsc = 6;   % 3 DM-RS symbols per RB and 2 PRBS bits needed per QPSK symbol
             cSeq = nrPRBS(cinit,nsc*[crboffset slen]);
             % Select onto the subset required to match the PRB
             cSeq = reshape(cSeq,nsc,[]);
             cSeq = cSeq(:,prb+1);
             % Create associated complex DM-RS symbols
             dmrssym = [dmrssym; nrSymbolModulate(cSeq(:),'QPSK')]; %#ok<AGROW>
             nsym = nsym + 1;  % Increment symbol number 
        end
            % Combine PDCCH with the grid
       ResourceElementGrids{pdcch(Ncch).BWP}(offset+redmrs) = ResourceElementGrids{pdcch(Ncch).BWP}(offset+redmrs) + dmrssym*db2mag(pdcch(Ncch).Power + pdcch(Ncch).PowerDMRS);    
    end 
end

for nch = 1:length(pdsch)
           % Get a copy of the current PDSCH channel parameters
        ch = pdsch(nch);

        % Only process configuration if enabled
        if ~ch.Enable 
            continue;
        end

        % Establish whether transport coding is enabled
        trcoding = ~isfield(ch,'EnableCoding') || ch.EnableCoding;
        
        % Check the referenced BWP index
        checkIndex('PDSCH',pdsch,nch,'BWP',bwp);

        % Expand the allocated slot sequence by the repetition period, across
        % the length of the waveform
        allocatedSlots = expandbyperiod(ch.AllocatedSlots,ch.AllocatedPeriod,waveconfig.NumSubframes,bwp(ch.BWP).SubcarrierSpacing);
        if any(ch.AllocatedPRB >= bwp(ch.BWP).NRB)
            error('The allocated PRB indices (0-based, largest value = %d) for PDSCH %d exceed the NRB (%d) for BWP %d.',max(ch.AllocatedPRB),nch,bwp(ch.BWP).NRB,ch.BWP);
        end        
        % Ensure that the allocated symbols for the slot are within a slot for the BWP CP
        symbperslot = symbolsPerSlot(bwp(ch.BWP));
        slotsymbs = ch.AllocatedSymbols(ch.AllocatedSymbols < symbperslot);
        if length(slotsymbs) ~= length(ch.AllocatedSymbols)
            warning('The slot-wise symbol allocation for PDSCH %d in BWP %d includes 0-based symbol indices which fall outside a slot (0...%d). Using only symbols within a slot.',nch,ch.BWP,symbperslot-1);
            ch.AllocatedSymbols = slotsymbs;
        end
    
        % Reserved PRB-level resources associated with SS burst
        rs = ssbreserved(bwp(ch.BWP).CarrierIdx);
        rs.PRB = rs.PRB - bwp(ch.BWP).RBOffset;
        ch.Reserved(end+1) = rs;         % Configure the channel with the pattern
        
        % Convert rate match pattern configurations into a format suitable 
        % for the hPDSCHResource function
        % 
        % Turn reserved CORESET indices into reserved patterns
        % with format, reserved = struct('Name',{},'PRB',{},'Symbols',{},'Period',{});
        for rmp = ch.RateMatch
            % Process CORESET ratematch pattern part
            if isfield(rmp,'CORESET')         
                for cidx = rmp.CORESET                
                    % Expand and project CORESET into the BWP 
                    % Pattern representation is single vector of PRB across all symbols  

                    % Check the CORESET index
                    if cidx < 1 || cidx > length(coreset)
                        error('For PDSCH %d, the ratematch CORESET index (%d) must be between 1 and the number of CORESET defined (%d)',...
                                         nch,cidx,length(coreset));
                    end

                    % Get a copy of the CORESET configuration
                    cs = coreset(cidx);

                    % Expand the allocated slots across the repetition period
                    rmallocatedSlots = expandbyperiod(cs.AllocatedSlots,cs.AllocatedPeriod,waveconfig.NumSubframes,bwp(ch.BWP).SubcarrierSpacing);
                    % All CORESET symbols in each allocated slot

                    % Expand to identify all symbols included in this CORESET sequence
                    slotsymbs = cs.AllocatedSymbols(cs.AllocatedSymbols+cs.Duration < symbperslot);
                    csetsymbols = expander(slotsymbs,cs.Duration);
                    rmallocatedsymbols = reshape(symbperslot*rmallocatedSlots+csetsymbols',1,[]);  % Use column expansion

                    % Turn the allocated PRB parameter into blocks of 6 RB/REG CCEs
                    allocatedPRB = expander(unique(6*fix(cs.AllocatedPRB/6)),6);

                    % Check that the associated PRB set fits within the associated BWP NRB
                    if max(allocatedPRB) >= bwp(ch.BWP).NRB
                        error('For PDSCH %d, the effective PRB allocation part of RateMatch CORESET %d (set of 6 PRB/REG CCEs, max PRB index = %d) exceeds the number of RB (%d) in BWP %d.',nch,cidx,max(allocatedPRB),bwp(ch.BWP).NRB,ch.BWP);
                    end

                    % Create reserved configuration structure and push it onto the copy of the PDSCH parameters
                    rs.Name = sprintf('Reserved for CORESET %d',cidx);
                    rs.PRB = allocatedPRB;           % Reserved PRB (0-based indices, defined as a vector or cell array)
                    rs.Symbols = rmallocatedsymbols; % OFDM symbols associated with reserved PRB (0-based indices, spanning one or more slots)
                    rs.Period = [];                  % Total number of slots in the pattern period (empty means don't cyclically repeat)
                    ch.Reserved(end+1) = rs;         % Configure the channel with the pattern
                end
            end
            % Process bitmap derived ratematch pattern part
            if isfield(rmp,'Pattern')
                for rmpat = rmp.Pattern
                    % Name this pattern for identification purposes
                    rs.Name = sprintf('Reserved for rate-matching pattern');

                    % Check that the associated PRB set fits within the associated BWP NRB
                    if max(rmpat.AllocatedPRB) >= bwp(ch.BWP).NRB
                        error('For PDSCH %d, the PRB allocation part of the RateMatch pattern exceeds the number of RB (%d) in BWP %d.',nch,bwp(ch.BWP).NRB,ch.BWP);
                    end
                    rs.PRB = rmpat.AllocatedPRB;
                    % Need to combine allocated symbols and allocated slots into a single list             
                    rmallocslots = expandbyperiod(rmpat.AllocatedSlots,rmpat.AllocatedPeriod,waveconfig.NumSubframes,bwp(ch.BWP).SubcarrierSpacing);
                    rmallocsymbols = reshape(symbperslot*rmallocslots + rmpat.AllocatedSymbols',1,[]);    
                    rs.Symbols = rmallocsymbols;    % OFDM symbols associated with reserved PRB (0-based indices, spanning one or more slots)

                    rs.Period = [];                % Total number of slots in the pattern period (empty means don't repeat)              
                    ch.Reserved(end+1) = rs;       % Configure the PDSCH channel with the pattern                 
                end       
            end
        end       

        % Waveform generation RE level processing 
        %
        % The hPDSCHResources uses a slot-level set of parameters so map the
        % relevant parameter from the waveform level down to the slot level
        nrb = bwp(ch.BWP).NRB;
        ch.PRBSet = ch.AllocatedPRB;
        ch.SymbolSet = ch.AllocatedSymbols;     
        crboffset = carriers(bwp(ch.BWP).CarrierIdx).RBStart + bwp(ch.BWP).RBOffset; % First CRB associated with start of BWP
        ch.PRBRefPoint = crboffset; % First CRB associated with start of BWP

        % Create a data source for this PDSCH sequence
        datasource = hVectorDataSource(ch.DataSource);

        % Configure the DL-SCH processing object for this PDSCH sequence
        if trcoding
            dlsch.TargetCodeRate = ch.TargetCodeRate;
        end
        
        % Loop over all the allocated slots
        for i = 1:length(allocatedSlots)

            % Get current slot number 
            s = allocatedSlots(i);
            ch.NSlot = s;

            % Create an empty slot grid to contain a single PDSCH instance
            slotgrid = zeros(12*nrb,symbperslot,ch.NLayers);

            % Get the slot-oriented PDSCH indices, DM-RS indices and DM-RS symbol values  
            [pdschREindices,dmrsREindices,dmrsSymbols,modinfo] = ...
                hPDSCHResources(struct('NRB',bwp(ch.BWP).NRB,'CyclicPrefix',bwp(ch.BWP).CyclicPrefix,...
                                       'SubcarrierSpacing',bwp(ch.BWP).SubcarrierSpacing),ch);
            if trcoding
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
            else
                % If transport coding is not enabled then get the codeword
                % directly from the data source
                codeword = datasource.getPacket(modinfo.G);
            end
            
            % PDSCH processing to create the PDSCH QAM symbols
            nID = ch.NID;
            symbols = nrPDSCH(codeword,ch.Modulation,ch.NLayers,nID,ch.RNTI);

            % Write the PDSCH and DM-RS symbols in the slot grid
            slotgrid(pdschREindices) = symbols*db2mag(ch.Power);
            slotgrid(dmrsREindices) = dmrsSymbols*db2mag(ch.Power+ch.PowerDMRS);

            % Combine PDSCH instance with the rest of the BWP grid
            ResourceElementGrids{ch.BWP}(:,s*symbperslot+(1:symbperslot),1:ch.NLayers) = ResourceElementGrids{ch.BWP}(:,s*symbperslot+(1:symbperslot),1:ch.NLayers) + slotgrid; 

        end
end



end

function s = getSlotIndex(AllocatedRegion, Period,  SpaceSize, SpaceNumber)
Space = SpaceSize*SpaceNumber;
PeriodNumber = ceil(Space/Period);
nAllocated = 0;

for k = 1 : length(AllocatedRegion)
    if AllocatedRegion(1, k) < Period
        nAllocated = nAllocated + 1;
    end
end

s = zeros(1, nAllocated*PeriodNumber);
nAllocated = 0;

for nPer = 0 : PeriodNumber - 1
    for k = 1 : length(AllocatedRegion)
      if AllocatedRegion(1, k) < Period
        nAllocated = nAllocated + 1;
        s(1, nAllocated) = AllocatedRegion(1, k) + nPer*Period;
      end
    end
end
    
end

function symb = getSymbIndex(SlotsIndexs, SymbInSlot, AllocatedSymb)
symb = zeros(1, length(SlotsIndexs)*length(AllocatedSymb));
    for i = 0 : length(SlotsIndexs) - 1
        for j = 1 : length(AllocatedSymb(1,:))
            symb(1, i * length(AllocatedSymb(1,:))  + j) = SlotsIndexs(i + 1) * SymbInSlot + AllocatedSymb(j);
        end
    end
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
