function out = bfRateMatchLDPC(in, outlen, rv, modulation, Nl,varargin)
% кодируемые биты in - двумерный массив [Nr, r], r - количесво блоков кодирования 
% Nr количесво бит в блоке r. Общее количество блоков обозначается буквой C, каждый 
% блок кодирования индивидуально согласовывается по скорости согласно п. 5.4.2 TS 38.212, Значение 
% LBRM = 1
% Выход out представляет собой соласованные по скорости биты - двумерный массив [Er,r], где 
% Er кодированные биты r-го блока
% outlen - длина блока после кодирования
% rv - версия резервирования выхода, RV =(0,1,2,3)
% modulation вид модуляции, возможные значения  {'pi/2-BPSK','QPSK','16QAM','64QAM','256QAM'}
% Nl - общее количество уровней передачи, связанных с транспортным блоком (1 ... 4). 

    narginchk(5,6);
    if nargin==5
        Nref = [];
    else
        Nref = varargin{1};
    end

[N,C] = size(in);
ZSets = getZsetsOfLDPC;

Nisvalid = -1;
for i = 1 : length(ZSets)
  for j = 1 : length(ZSets{i})
    if N == ZSets{i}(j)*50
      bgn = 2;
      ncwnodes = 50;
      Nisvalid = 1;
      break;
    end
    if N == ZSets{i}(j)*66
      bgn = 1;
      ncwnodes = 66;
      Nisvalid = 1;
      break;
    end
  end
  if Nisvalid == 1
    break;
  end
end 

if Nisvalid == -1
    error('bfRateMatchLDPC: invalid N');
end

    Zc = N/ncwnodes;

% определяем порядок модуляции
    switch modulation
    case {'pi/2-BPSK', 'BPSK'}
        Qm = 1;
    case 'QPSK'
        Qm = 2;
    case '16QAM'
        Qm = 4;
    case '64QAM'
        Qm = 6;
    otherwise   % '256QAM'
        Qm = 8;
    end


if ~isempty(Nref)
   Ncb = min(N,Nref);
else   
   Ncb = N;
end

% находим k0 из таблицы 5.4.2.1-2 TS 38.212
    if bgn == 1
        if rv == 0
            k0 = 0;
        elseif rv == 1
            k0 = floor(17*Ncb/N)*Zc;
        elseif rv == 2
            k0 = floor(33*Ncb/N)*Zc;
        else % rv is equal to 3
            k0 = floor(56*Ncb/N)*Zc;
        end
    else
        if rv == 0
            k0 = 0;
        elseif rv == 1
            k0 = floor(13*Ncb/N)*Zc;
        elseif rv == 2
            k0 = floor(25*Ncb/N)*Zc;
        else % rv is equal to 3
            k0 = floor(43*Ncb/N)*Zc;
        end
    end

% 
out = [];
  for r = 0 : C - 1 
    if r <= C - mod(outlen/(Nl*Qm), C) - 1
        E = Nl*Qm*floor(outlen/(Nl*Qm*C));
    else
        E = Nl*Qm*ceil(outlen/(Nl*Qm*C));
    end
    
    % 5.4.2.1 TS 38.212
    k = 0;
    j = 0;
    e = zeros(E,1,class(in(:, r+1)));
    while k < E
        if in(mod(k0+j,Ncb)+1, r+1) ~= -1    
            e(k+1) = in(mod(k0+j,Ncb)+1, r+1);
            k = k+1;
        end
        j = j+1;
    end
    
    % Чередование битов п. 5.4.2.2 TS 38.212
    f = zeros(E, 1);
    for j = 0 : E/Qm - 1
      for i = 0 : Qm - 1
        f(i+j*Qm + 1) = e (i*E/Qm + j + 1);
      end
    end
    % объединение кодовых блоков п.5.5 TS 38.212
    out = [out; f];
  end
  
end


