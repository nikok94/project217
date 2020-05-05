function out_CRC = bfCRCEncode(in, polynom, mask)
%% out_CRC = bfGetCRC(in, polynom) возвращает входные значения in (вектор столбец) с
% добавлеными битами CRC out вектор столбец длинною length(in) + lengthPoly
%% Для CRC определено 6 возможных полиномов п. 5.1 TS 38.212
% ниже описан процесс выбора соответствующего полинома и его длины
    switch polynom
        case '6'
            Poly = [1 1 0 0 0 0 1]';
            lengthPoly = 6;
        case '11'
            Poly = [1 1 1 0 0 0 1 0 0 0 0 1]';
            lengthPoly = 11;
        case '16'
            Poly = [1 0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 1]';
            lengthPoly = 16;
        case {'24a','24A'}
            Poly = [1 1 0 0 0 0 1 1 0 0 1 0 0 1 1 0 0 1 1 1 1 1 0 1 1]';
            lengthPoly = 24;
        case {'24b','24B'}
            Poly = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 1 1]';
            lengthPoly = 24;
        otherwise % {'24c','24C'}
            Poly = [1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1]';
            lengthPoly = 24;
    end
    %% Заведем расширенную до длины length(in) + lengthPoly переменную 
    % и скопируем входные значения в начало
    inWid = zeros(length(in) + lengthPoly, 1);
    inWid(1 : length(in) ,1)= in(:,1);
    
    %% Процедура нахождения CRC сводится к нахождению остатка от деления
    % двоичного значения заданного вектором in на двоичное значение 
    % двоичного вектора представленного значением Poly
    % деление осуществляется в двоичном коде и необходимо оперировать
    % знаковыми числами, поэтому заведем знаковую переменную длинной 
    % lengthPoly + 1, младший индекс отвечает за знак
    
    crc = zeros(lengthPoly + 1,1); % знаковая переменная числителя 
    crc(2 : end) = inWid(1 : lengthPoly,1);
    
    for i = 1:length(in)
        numerator = [crc(2:end); inWid(i+lengthPoly)];
        if numerator(1) == 1
            crc = rem(Poly+numerator,2);
        else
            crc = numerator;
        end
    end

    out_CRC = [in; crc(2:end)];
    if mask
        maskBits = rem(floor(double(mask)*pow2(1-lengthPoly:0)),2)';
        out_CRC(end-lengthPoly+1:end,:) = xor(out_CRC(end-lengthPoly+1:end,:), ...
            repmat(maskBits>0,[1 1]));
    end
end