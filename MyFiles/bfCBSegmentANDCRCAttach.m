function [Zc, iLS, c_out] = bfCBSegmentANDCRCAttach(dataIn, LDPCgraph)
% данна€ функци€ определена в п.5.2.2 TS 38.212
% входна€ битова€ последовательность dataIn это одномерный массив в виде
% вектора столбца размер length(dataIn(:,1)) > 0. ≈сли размер входной
% последовательности больше чем максимальный размер кодового блока Kcb 
% (Kcb = 8448 дл€ графа 1 и  Kcb = 3840 дл€ графа 2)дл€ 
% данного графа, то выполн€етс€ сегментаци€ входной последовательности и к
% каждому сегменту добавл€етс€ CRC длины 24
% Zc - минимальное значение Zc из таблицыв п.5.3.2-1 TS 38.212 Zlist 
% при котором Kb*Zc >= Kd
% iLS - индекс значени€ Zc
B = length(dataIn);

if LDPCgraph == 1
    % ƒл€ базового графа 1 LDPC максимальный размер кодового блока составл€ет:
    Kcb = 8448;
else
    % ƒл€ базового графа 2 LDPC максимальный размер кодового блока составл€ет:
    Kcb = 3840;
end

% определ€ем количество кодовых блоков C 

if B <= Kcb
    % если размер входной последовательности бит меньше чем максимальный размер 
    % кодового блока Kcb, тогда   
    L = 0;
    C = 1;
    Bd = B;
else
    L = 24;
    C = ceil(B/(Kcb-L));
    Bd = B + C*L;
end

% результат c_out это многомерный массив [Kr, r], где 1 <= r < C + 1 (т.к.
% в матлабе отсчет от 1) Kr = K, K определ€етс€ следующим образом:

Kd = ceil(Bd/C);

if LDPCgraph == 1
  Kb = 22;
else
  if B > 640
    Kb = 10;
  elseif B > 560
    Kb = 9;
  elseif B > 192
    Kb = 8;
  else
    Kb = 6;
  end
end

% получаем таблицу Zlist размеров LDPC определенную в п.5.3.2-1
% TS 38.212
Zlist = getZsetsOfLDPC;

% находим минимальное значение Zc из таблицы Zlist при котором Kb*Zc >= Kd
Zc = -1;
for i = 1 : length(Zlist)
    for j = 1 : length(Zlist{i})
        zl = Zlist{i}(j);
        if  zl*Kb >= Kd
            if Zc == -1
                Zc = zl;
                iLS = i;
            else
                if Zc > zl
                    Zc = zl;
                    iLS = i;
                end 
            end
        end 
    end
end 

if LDPCgraph == 1
  K = 22*Zc;
else
  K = 10*Zc;
end

s = 1;
c_out = zeros(K, C);
  for r = 1 : C
    for k = 1 : Kd - L
        c_out(k, r) = dataIn(s);
        s = s + 1;
    end
    if C > 1
        % дл€ каждого блока вычисл€ем CRC, использу€ полином "24B"
        cCRC = bfCRCEncode(c_out(1 : Kd - L, r), '24B', 0);
        c_out(1 : Kd, r) = cCRC;
    end
    
    c_out(Kd + 1 : K, r) = -1;
  end
end