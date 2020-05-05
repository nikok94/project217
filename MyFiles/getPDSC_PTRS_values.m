function  [kIndex, lIndex]= getPDSC_PTRS_values(DMRSLPosition, PDSCHx, DMRSAntennaPort, l0_dmrs_index, l0_dmrs_value)
i = 0;
lref = 0;
l = min(PDSCHx.AllocatedSymbols);
PTRSInd = [];
while lref + i*PDSCHx.PTRSLength + l<= max(PDSCHx.AllocatedSymbols)
ind = (max(lref + (i-1)*PDSCHx.PTRSLength + 1, lref) : lref + i*PDSCHx.PTRSLength) + l;
   if overlap(DMRSLPosition, ind)
       i = 1;
       lref = max(DMRSLPosition);
       ind = (max(lref + (i-1)*PDSCHx.PTRSLength + 1, lref) : lref + i*PDSCHx.PTRSLength) + l;
   end
   PTRSInd = [PTRSInd, ind(end)];
   i = i + 1;
end

lIndex = getIntersection(PTRSInd, PDSCHx.AllocatedSymbols);

kREref = get_kREref(PDSCHx.DMRSConfigurationType,DMRSAntennaPort,PDSCHx.DMRSAdditionalPosition);
NRB = length(PDSCHx.AllocatedPRB);
if mod(NRB,PDSCHx.kPTRS) == 0
    kRBref = mod(PDSCHx.RNTI, PDSCHx.kPTRS);
else
    kRBref = mod(PDSCHx.RNTI, mod(NRB,PDSCHx.kPTRS));
end
% вычисляем количество индексов задейсвованных для PT-RS
kIndex = zeros(ceil((length(PDSCHx.AllocatedPRB) - kREref)/PDSCHx.kPTRS), 1);

% находим индексы PT-RS в пределах между минимальным и максимальным
% значениями ресурс блоков выделенных для PDSCH
for p = 0 : length(kIndex) - 1
    kIndex(p+1) = kREref + (p*PDSCHx.kPTRS + kRBref)*12;
end

end

function out = getIntersection(a, b)
out = [];
    for i = 1 : length(a)
        for j = 1 : length(b)
           if b(j) == a(i)
            out = [out, a(i)];
           end
        end
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