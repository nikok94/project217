function di = bfModulator(bi,mode)
% di = NRSymbolMode(mode, bi) отображает битовую последовательность,
% записанную в кодовом слове 'bi', в сложные символы модуляции в
% соответствии с разделом 5.1 стандарта TS 38.211.
% Схема модуляции, 'mode' должна быть одной из 'pi / 2-BPSK', 'BPSK',
% 'QPSK', '16QAM','64QAM', '256QAM'. 'bi' вектор столбец.
% Пример:
% Генерация случайной последовательности двоичных значений длиной 256. 
% Генерация модулированных символов с использованием модуляции QPSK.
% bi=randi([0 1],256,1);
% y=NRSymbolMode('QPSK',bi);

tmp=complex(0+1i*0);
switch mode
    % pi/2-BPSK
    case 'pi/2-BPSK'
        di=zeros(length(bi),1,'like',tmp);
        for i= 0:length(bi)-1
            di(i+1,1)=exp(1i*pi/2*mod(i,2))/sqrt(2)*((1-2*bi(i+1,1))+1i*(1-2*bi(i+1,1)));
        end
%         'BPSK'
    case 'BPSK' 
        di=zeros(length(bi),1,'like',tmp);
        for i=0:length(bi)-1
            di(i+1,1)=1/sqrt(2)*((1-2*bi(i+1,1))+1i*(1-2*bi(i+1,1)));
        end
%         'QPSK'
    case 'QPSK' 
        if mod(length(bi),2) == 0
            di=zeros(length(bi)/2,1,'like',tmp);
            for i= 0 : length(bi)/2-1
                di(i+1,1)=1/sqrt(2)*((1-2*bi(2*i + 1,1))+1i*(1-2*bi(2*i+2,1)));
            end
        else 
            disp('неверный тип данных');
        end
%         '16QAM'
    case '16QAM'
        if mod(length(bi),4) == 0
            di=zeros(length(bi)/4,1,'like',tmp);
            for i= 0 : length(bi)/4-1
                di(i+1,1)=1/sqrt(10)*((1-2*bi(4*i+1,1))*(2-(1-2*bi(4*i+2+1,1)))+1i*(1-2*bi(4*i+1+1,1))*(2-(1-2*bi(4*i+3+1,1))));
            end
        else
            disp('неверный тип данных');
        end
%         '64QAM'
    case '64QAM'
        if mod(length(bi),6) == 0
            di=zeros(length(bi)/6,1,'like',tmp);
            for i= 0 : length(bi)/6-1
                di(i+1,1)=1/sqrt(42)*((1-2*bi(6*i+1,1))*(4-(1-2*bi(6*i+3,1))*(2-(1-2*bi(6*i+5,1))))+1i*(1-2*bi(6*i+2,1))*(4-(1-2*bi(6*i+4,1))*(2-(1-2*bi(6*i+6,1)))));
            end
        else
            disp('неверный тип данных');
        end 
%         '256QAM'
    case '256QAM'
     if mod(length(bi),8) == 0 
        di=zeros(length(bi)/8,1,'like',tmp);
        for i= 0 : length(bi)/8-1
            di(i+1,1)=1/sqrt(170)*((1-2*bi(8*i+1,1))*(8-(1-2*bi(8*i+3,1))*(4-(1-2*bi(8*i+5,1))*(2-(1-2*bi(8*i+7,1)))))+1i*(1-2*bi(8*i+2,1))*(8-(1-2*bi(8*i+4,1))*(4-(1-2*bi(8*i+6,1))*(2-(1-2*bi(8*i+8,1))))));  
        end
     else
         disp('неверный тип данных');
     end   
    otherwise  disp('неверный тип модуляции, введите из следующих типов(pi/2-BPSK, BPSK, QPSK, 16QAM, 64QAM, 256QAM)');
end








