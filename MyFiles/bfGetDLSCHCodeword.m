function codeword = bfGetDLSCHCodeword(dBits, codeRate, rv, modulation, nlayers, outlen)

A = length(dBits(:,1));
% добавляем к входным битам CRC п.7.2.1 TS 38_212
if  A > 3824
    dBits_CRC = bfCRCEncode(dBits(:,1), '24A', 0);
else
    dBits_CRC = bfCRCEncode(dBits(:,1), '16', 0);
end

% Выбор базового графа LDPC п.7.2.2 TS 38_212
if (A <= 292) || ((A <= 3824) && (codeRate <= 0.67)) || codeRate <= 0.25
    ldpcGraph = 2;
else
    ldpcGraph = 1;
end

% Разбиваем на блоки кодирования и присоединеняем CRC п.7.2.3 TS 38_212
[Zc, iLS, c] = bfCBSegmentANDCRCAttach(dBits_CRC, ldpcGraph);

% LDPC кодирование
enc = nrLDPCEncode(c,ldpcGraph);
%enc1 = bfLDPCEncode(c, ldpcGraph, iLS, Zc);

% Согласование скорости и объединение кодовых блоков
codeword = bfRateMatchLDPC(enc,outlen,rv,modulation,nlayers);

end

