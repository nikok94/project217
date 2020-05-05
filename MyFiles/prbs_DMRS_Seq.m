function symbols_DMRS = prbs_DMRS_Seq(PRBSet, PRBRefPoint, ndmrsre, symbperslot, Nslot, nsymbol, ncsid, DMRSdownlink_r16, nidnscid,  lambda)
%% Генерация сигнала DMRS для PDSCH
% На вход функции поступает: 
%   k_dash - k' в стд.38.211 разд.7.4.1.1.2
%   PRBSet - PRB, выделенные каналу;
%   PRBRefPoint - Индекс PRB, к которому относится PRBSet.
%   ndmrssc - RE индексы в PRB. Зависит от типа конфигурации DMRS, если
% DMRSConfigurationType1, то ndmrssc = 6, если DMRSConfigurationType2, то
% ndmrssc = 4;
%   symbperslot - количество символов OFDM в слоте (12 или 14);
%   Nslot - номер слота передачи;
%   scs - разнесение поднесущих;
%   nsymbol - номер текущего символа DMRS;
%   ncsid - задается полем инициализации последовательности DM-RS 1 бит в DCI, связанном с передачей PDSCH, 
% если используется формат DCI 1_1 или 1_2 в (стд.38.212 разд.7.3.1.2.2 и разд. 7.3.1.2.3), 
% в противном случае ncsid = 0;
%   DMRSdownlink_r16 - параметр более высокого уровня, задается в
% DMRS-DownlinkConfig, если задан то задаем значение 1, если нет 0;
%   nidnscid - принимает значение {0,1,...65535}, задаются параметрами более высокого уровня scrabingID0 и scrabingID1 соответственно 
% в IE DMRS-DownlinkConfig, если он предоставляется, и PDSCH планируется с помощью PDCCH с использованием 
% формата DCI 1_1 или 1_2 с CRC, скремблированным C-RNTI, MCS-C-RNTI, или
% CS-RNTI. Или задается параметром scrabingID0 верхнего уровня в IE DMRS-DownlinkConfig, если он предоставляется, 
% и PDSCH планируется с помощью PDCCH с использованием формата DCI 1_0 с
% CRC, скремблированным C-RNTI, MCS-C-RNTI или CS-RNTI. Иначе nidnscid = N_Cell_ID.
%   lambda - является группой CDM, определенной в стд.38.211 разд.
% 7.4.1.1.2.Принимает значения (0,1,2).
% На выходе функции имеем значения DMRS для одного символа OFDM.


	if DMRSdownlink_r16 == 1
        if (lambda == 1)
            nlambdscid = 1 - ncsid;
        else
            nlambdscid = ncsid;
        end
        lamdaX = lambda;
    else
        nlambdscid = ncsid;
        lamdaX = 0;
	end   

    if ~isempty(PRBSet)
        % Сгенерируем инициальзирующую последовательность c_init в соответствии со
        % стд.38.211 разд. 7.4.1.1.1
        [minprb,maxprb] = bounds(PRBSet);
        c_init = mod(2^17*(symbperslot*Nslot + nsymbol + 1)*(2*nidnscid + 1) + 2^17*floor(lamdaX/2) + 2*nidnscid + nlambdscid,2^31);
        E = 2*ndmrsre*((maxprb-minprb+1) + (PRBRefPoint + minprb)); 
        c = scrambling(c_init,E);
        DMRS_Seq = zeros(length(c)/2, 1);
        for i = 0 : length(c)/2 - 1
            DMRS_Seq(i + 1) = 1/sqrt(2)*(1 - 2*c(2*i+1))+1i*1/sqrt(2)*(1 - 2*c(2*i+2));
        end
        symbols_DMRS = DMRS_Seq(end - (maxprb-minprb+1)*ndmrsre + 1 : end);
        
        
    else
        symbols_DMRS = complex([]);
    end
end