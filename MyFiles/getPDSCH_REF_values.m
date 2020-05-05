function [dmrsRESymbols, dmrsREIndex, ptrsRESymbols, ptrsREIndex] = getPDSCH_REF_values(PDSCHx, NRB, symbperslot, DMRSdownlink_r16)
% PDSCHMappingType  тип отображения на сетку ресурсов, может принимать
%                   следуюющие значения 'A' или 'B' TS 38.211 п. 7.4.1.1.2
% dmrsType          тип конфигурации dmrs, может принимать значения 1 или 2

   if PDSCHx.PDSCHMappingType ~= 'A' && PDSCHx.PDSCHMappingType ~= 'B'
       error('GetDMRSPositionTable: PDSCHMappingType error');
   end
      
   if PDSCHx.PDSCHMappingType == 'A'
       typeB = 0;
       l_0 = PDSCHx.DMRSTypeAPosition;
   else
       typeB = 1;
       l_0 = 0;
   end

l = PDSCHx.AllocatedSymbols(1);
m = l;
for i = 2 : length(PDSCHx.AllocatedSymbols)
    if PDSCHx.AllocatedSymbols(i) > m
        m = PDSCHx.AllocatedSymbols(i);
    end
    if PDSCHx.AllocatedSymbols(i) < l
        l = PDSCHx.AllocatedSymbols(i);
    end
end
if typeB == 0
    l = 0;
end

ld = m - l + 1;

DMRSPositions = GetDMRSPositionTable(typeB, ld, PDSCHx);
if isempty(DMRSPositions)
   dmrsREIndex = [];
   dmrsRESymbols = [];
   return;
end

if PDSCHx.PDSCHMappingType == 'A'
    DMRSPositions(1) = l_0;
end

switch PDSCHx.DMRSConfigurationType
    case 1 
        ConfTypeFactor = 12/4;
        ndmrsre = 12/4*2;
    case 2 
        ConfTypeFactor = 12/6;
        ndmrsre = 12/6*2;
    otherwise 
        error('DMRSConfigurationType invalid');  
end
% вычисляем количество dmrs символов в слоте
dmrsREIndex = zeros(length(PDSCHx.AllocatedPRB)*ConfTypeFactor*PDSCHx.DMRSLength*length(DMRSPositions), 1);
dmrsRESymbols = dmrsREIndex;

kREref = get_kREref(PDSCHx.DMRSConfigurationType,PDSCHx.NLayer,PDSCHx.DMRSAdditionalPosition);

if mod(length(PDSCHx.AllocatedPRB),PDSCHx.kPTRS) == 0
    kRBref = mod(PDSCHx.RNTI, PDSCHx.kPTRS);
else
    kRBref = mod(PDSCHx.RNTI, mod(length(PDSCHx.AllocatedPRB),PDSCHx.kPTRS));
end

PTRS_kIndex = zeros(ceil((length(PDSCHx.AllocatedPRB) - kREref)/PDSCHx.kPTRS), 1);
PTRS_kValue = PTRS_kIndex;

params = getDMRSParameters(PDSCHx.NLayer, PDSCHx.DMRSConfigurationType);
lambda = params{1};
delta = params{2}; 

wf = params{3};
wt = params{4};

prbMin = min(PDSCHx.PRBSet);
nn = 1;
iPTRS = 0;
t = 0;
dmrsLPosition = zeros(1, PDSCHx.DMRSLength);
for ndmrs = DMRSPositions
    for lpos = 0 : PDSCHx.DMRSLength - 1
        l = lpos + ndmrs;
        dmrsLPosition(1, lpos + 1) = l;
        symbols_DMRS = prbs_DMRS_Seq(PDSCHx.PRBSet, PDSCHx.PRBRefPoint, ndmrsre, symbperslot, PDSCHx.NSlot, l, PDSCHx.NSCID, DMRSdownlink_r16, PDSCHx.NIDNSCID,  lambda);
        for p = PDSCHx.AllocatedPRB
            for pp = 0 : ConfTypeFactor - 1
                for kk = 0 : 1
                    n = p*ConfTypeFactor + pp;
                    dmrsi = (p-prbMin)*ConfTypeFactor + pp;
                    if PDSCHx.DMRSConfigurationType == 1
                        k = 4*n + 2*kk + delta;
                    else
                        k = 6*n + kk + delta;
                    end
                    if l == l_0
                        % если текущая позиция символа равна l0 то
                        % находим индексы PTRS_kIndex для PTRS согласно п. 7.4.1.2.2 TS 38.211
                        % поочереди начиная с самого меньшего значения 
                        % находим символы DMRS соответсвующие индексам PTRS_kIndex
                        while kREref + (iPTRS*PDSCHx.kPTRS + kRBref)*12 < k
                            iPTRS = iPTRS + 1;
                        end
                        PRTSCurInd = kREref + (iPTRS*PDSCHx.kPTRS + kRBref)*12;
                        if k == PRTSCurInd
                           PTRS_kIndex(t + 1) = PRTSCurInd + 1;
                           PTRS_kValue(t + 1) = symbols_DMRS(2*dmrsi + kk + 1);
                           iPTRS = iPTRS + 1;
                           t = t + 1;
                        end
                    end
                    dmrsREIndex(nn) = l * NRB * 12 + k + 1;
                    dmrsRESymbols(nn) = symbols_DMRS(2*dmrsi + kk + 1)*wf(kk + 1)*wt(lpos+1);
                    nn = nn + 1;
                end
            end
        end
    end
end

% находим в каких символах должны находиться PTRS  
i = 0;
lref = 0;
l = min(PDSCHx.AllocatedSymbols);
lIndex = [];
while lref + i*PDSCHx.PTRSLength + l<= max(PDSCHx.AllocatedSymbols)
ind = (max(lref + (i-1)*PDSCHx.PTRSLength + 1, lref) : lref + i*PDSCHx.PTRSLength) + l;
   if overlap(dmrsLPosition, ind)
       i = 1;
       lref = max(dmrsLPosition);
       ind = (max(lref + (i-1)*PDSCHx.PTRSLength + 1, lref) : lref + i*PDSCHx.PTRSLength) + l;
   end
   lIndex = [lIndex, ind(end)];
   i = i + 1;
end

PTRS_lIndex = getIntersection(lIndex, PDSCHx.AllocatedSymbols);
ptrsREIndex = zeros(length(PTRS_lIndex) * length(PTRS_kIndex), 1);
ptrsRESymbols = ptrsREIndex;

k = 1;
for i = PTRS_lIndex
    t = 1;
    for j = PTRS_kIndex'
        ptrsREIndex(k) = i * NRB * 12 + j;
        ptrsRESymbols(k) = PTRS_kValue(t);
        t = t + 1;
        k = k + 1;
    end
end


end

function out = getIntersection(a, b)
% функция getIntersection возвращает только те элементы, которые 
% имеются в сразу в двух массивах a и b

out = [];
    for i = 1 : length(a)
        for j = 1 : length(b)
           if b(j) == a(i)
            out = [out, a(i)];
           end
        end
    end
end
function symbols = prbsDMRSSequence(pxsch,ndmrssc,prbset,prbrefpoint,nslot,nsymbol,ldash,symbperslot) %#ok<INUSL>
    
    % Cache the scrambling IDs
    nidnscid = pxsch.NIDNSCID;
    nscid = pxsch.NSCID;
    
    if ~isempty(prbset)
        % Generate PRBS for DM-RS sequence which covers the PRB allocation set range 
        [minprb,maxprb] = bounds(prbset);
        cinit = mod(2^17*(symbperslot*nslot + nsymbol + 1)*(2*nidnscid + 1) + 2*nidnscid + nscid,2^31);
        prbs = reshape(nrPRBS(cinit,2*ndmrssc*[prbrefpoint+minprb maxprb-minprb+1]),2*ndmrssc,[]);
    
        % Extract PRBS values associated with PRB and turn into complex DM-RS symbols
        bpsk = 1/sqrt(2)*(1-2*reshape(prbs(:,prbset-minprb+1),2,[])');
        symbols = complex(bpsk(:,1),bpsk(:,2));
    else
        symbols = complex([]);
    end
end  

function out = overlap(array, value)
out = 0;
    for i = 1 : length(value)
        for j = 1 : length(array)
           if array(j) == value(i)
            out = 1;
            return;
           end
        end
    end
end

function DMRSPositions = GetDMRSPositionTable(typeB, symbDuration, PDSCHx)
% TS 38.211 п. 7.4.1.1.2
% PDSCHMappingType = ['A', 'B']
% DMRSLength = [1 , 2]
% DMRSAdditionalPosition = [0:3]
% количество доступных символов
   
        % Type A, single-symbol, 0,1,2,3 *additional* symbols
        dmrs_singleA = {
            [],[],  [],  [];                %  1 symbol duration
            [],[],  [],  [];                %  2 symbol duration
            0,  0,  0,    0;                %  3 symbol duration
            0,  0,  0,    0;                %  4 symbol duration
            0,  0,  0,    0;                %  5 symbol duration
            0,  0,  0,    0;                %  6 symbol duration
            0,  0,  0,    0;                %  7 symbol duration
            0,  [0,7],  [0,7],  [0,7];      %  8 symbol duration
            0,  [0,7],  [0,7],  [0,7];      %  9 symbol duration
            0,  [0,9], [0,6,9], [0,6,9];    % 10 symbol duration
            0,  [0,9], [0,6,9], [0,6,9];    % 11 symbol duration
            0,  [0,9], [0,6,9], [0,5,8,11]; % 12 symbol duration
            0, [0,-1], [0,7,11],[0,5,8,11]; % 13 symbol duration,  -1 represents l_1 = 11/12
            0, [0,-1], [0,7,11],[0,5,8,11]; % 14 symbol duration,  -1 represents l_1 = 11/12
        };
        % Type B, single-symbol, 0,1,2,3 *additional* symbols
        dmrs_singleB = {
            [],[],[],[];         %  1 symbol duration
             0, 0, 0, 0;         %  2 symbol duration
            [],[],[],[];         %  3 symbol duration
             0, 0, 0, 0;         %  4 symbol duration
            [],[],[],[];         %  5 symbol duration 
             0,[0,4],[0,4],[0,4];%  6 symbol duration (extended CP, half slot)
             0,[0,4],[0,4],[0,4];%  7 symbol duration (normal CP, half slot)
            [],[],[],[];         %  8 symbol duration
            [],[],[],[];         %  9 symbol duration
            [],[],[],[];         % 10 symbol duration
            [],[],[],[];         % 11 symbol duration
            [],[],[],[];         % 12 symbol duration
            [],[],[],[];         % 13 symbol duration
            [],[],[],[];         % 14 symbol duration
        }; 
        % Type A, double-symbol, 0,1(,2) *additional* symbol *pairs*
        dmrs_doubleA = {
            [],[];     %  1 symbol duration
            [],[];     %  2 symbol duration
            [],[];     %  3 symbol duration
             0, 0;     %  4 symbol duration
             0, 0;     %  5 symbol duration
             0, 0;     %  6 symbol duration
             0, 0;     %  7 symbol duration
             0, 0;     %  8 symbol duration
             0, 0;     %  9 symbol duration
             0,[0,8];  % 10 symbol duration
             0,[0,8];  % 11 symbol duration
             0,[0,8];  % 12 symbol duration
             0,[0,10]; % 13 symbol duration
             0,[0,10]; % 14 symbol duration
        };
        % Type B, double-symbol, 0,1(,2) *additional* symbol *pairs*
        dmrs_doubleB = {
            [],[];    %  1 symbol duration    
            [],[];    %  2 symbol duration
            [],[];    %  3 symbol duration
            [],[];    %  4 symbol duration
            [],[];    %  5 symbol duration       
             0, 0;    %  6 symbol duration
             0, 0;    %  7 symbol duration
            [],[];    %  8 symbol duration
            [],[];    %  9 symbol duration
            [],[];    % 10 symbol duration
            [],[];    % 11 symbol duration
            [],[];    % 12 symbol duration
            [],[];    % 13 symbol duration
            [],[];    % 14 symbol duration
        };   
   dmrs_add_pos = { dmrs_singleA, dmrs_doubleA;
                    dmrs_singleB, dmrs_doubleB }; 
   if isempty(PDSCHx.DMRSLength)
       PDSCHx.DMRSLength = 1;
   end
                
   positionstable = dmrs_add_pos{typeB+1, PDSCHx.DMRSLength};
    
   if PDSCHx.DMRSAdditionalPosition < size(positionstable,2)
       DMRSPositions = positionstable{symbDuration,PDSCHx.DMRSAdditionalPosition + 1};
   else
       DMRSPositions = [];
   end
   
   if ~isempty(DMRSPositions) && DMRSPositions(end)==-1
       DMRSPositions(end) = 9 + PDSCHx.DMRSTypeAPosition;
   end
end
function out = get_kREref(DMRSConfigurationType,DMRSAntennaPort,DMRSAdditionalPosition)
cfg1 = {    [0, 2, 6, 8];
            [2, 4, 8, 10];
            [1, 3, 7, 9];
            [3, 5, 9, 11] 
    };
cfg2 = {    [0, 1, 6, 7];
            [1, 6, 7, 0];
            [2, 3, 8, 9];
            [3, 8, 9, 2];
            [4, 5, 10, 11];
            [5, 10, 11, 4]
            
    };

switch DMRSConfigurationType
    case 1 
        if DMRSAntennaPort > 3
            error('DMRSAntennaPort error');
        end
        out = cfg1{DMRSAntennaPort + 1}(DMRSAdditionalPosition + 1);
    case 2
        if DMRSAntennaPort > 5
            error('DMRSAntennaPort error');
        end
        out = cfg2{DMRSAntennaPort + 1}(DMRSAdditionalPosition + 1);
    otherwise
        error('DMRSConfigurationType error');
        
end
end
function params = getDMRSParameters(NLayer, DMRSConfigurationType)
 config1params = {
                   {0, 0, [1 , 1], [1, 1]};   % p=1000
                   {0, 0, [1 , -1], [1, 1]};  % p=1001
                   {1, 1, [1 , 1], [1, 1]};   % p=1002
                   {1, 1, [1 , -1], [1, 1]};  % p=1003
                   {0, 0, [1 , 1], [1, -1]};  % p=1004
                   {0, 0, [1 , -1], [1, -1]}; % p=1005
                   {1, 1, [1 , 1], [1, -1]};  % p=1006
                   {1, 1, [1 , -1], [1, -1]}  % p=1007
            };
 config2params = {
                   {0, 0, [1 , 1], [1, 1]};   % p=1000
                   {0, 0, [1 , -1], [1, 1]};  % p=1001
                   {1, 2, [1 , 1], [1, 1]};   % p=1002
                   {1, 2, [1 , -1], [1, 1]};  % p=1003
                   {2, 4, [1 , 1], [1, 1]};   % p=1004
                   {2, 4, [1 , -1], [1, 1]};  % p=1005
                   {0, 0, [1 , 1], [1, -1]};  % p=1006
                   {0, 0, [1 , -1], [1, -1]}; % p=1007
                   {1, 2, [1 , 1], [1, -1]};  % p=1008
                   {1, 2, [1 , -1], [1, -1]}; % p=1009
                   {2, 4, [1 , 1], [1, -1]};  % p=1010
                   {2, 4, [1 , -1], [1, -1]}  % p=1011
            };
 switch DMRSConfigurationType
    case 1 
        params = config1params{NLayer};
    case 2 
        params = config2params{NLayer};
    otherwise 
        error('DMRSConfigurationType invalid');  
 end
  
end
