function symbols_DMRS = prbs_DMRS_Seq(PRBSet, PRBRefPoint, ndmrsre, symbperslot, Nslot, nsymbol, ncsid, DMRSdownlink_r16, nidnscid,  lambda)
%% ��������� ������� DMRS ��� PDSCH
% �� ���� ������� ���������: 
%   k_dash - k' � ���.38.211 ����.7.4.1.1.2
%   PRBSet - PRB, ���������� ������;
%   PRBRefPoint - ������ PRB, � �������� ��������� PRBSet.
%   ndmrssc - RE ������� � PRB. ������� �� ���� ������������ DMRS, ����
% DMRSConfigurationType1, �� ndmrssc = 6, ���� DMRSConfigurationType2, ��
% ndmrssc = 4;
%   symbperslot - ���������� �������� OFDM � ����� (12 ��� 14);
%   Nslot - ����� ����� ��������;
%   scs - ���������� ����������;
%   nsymbol - ����� �������� ������� DMRS;
%   ncsid - �������� ����� ������������� ������������������ DM-RS 1 ��� � DCI, ��������� � ��������� PDSCH, 
% ���� ������������ ������ DCI 1_1 ��� 1_2 � (���.38.212 ����.7.3.1.2.2 � ����. 7.3.1.2.3), 
% � ��������� ������ ncsid = 0;
%   DMRSdownlink_r16 - �������� ����� �������� ������, �������� �
% DMRS-DownlinkConfig, ���� ����� �� ������ �������� 1, ���� ��� 0;
%   nidnscid - ��������� �������� {0,1,...65535}, �������� ����������� ����� �������� ������ scrabingID0 � scrabingID1 �������������� 
% � IE DMRS-DownlinkConfig, ���� �� ���������������, � PDSCH ����������� � ������� PDCCH � �������������� 
% ������� DCI 1_1 ��� 1_2 � CRC, ���������������� C-RNTI, MCS-C-RNTI, ���
% CS-RNTI. ��� �������� ���������� scrabingID0 �������� ������ � IE DMRS-DownlinkConfig, ���� �� ���������������, 
% � PDSCH ����������� � ������� PDCCH � �������������� ������� DCI 1_0 �
% CRC, ���������������� C-RNTI, MCS-C-RNTI ��� CS-RNTI. ����� nidnscid = N_Cell_ID.
%   lambda - �������� ������� CDM, ������������ � ���.38.211 ����.
% 7.4.1.1.2.��������� �������� (0,1,2).
% �� ������ ������� ����� �������� DMRS ��� ������ ������� OFDM.


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
        % ����������� ���������������� ������������������ c_init � ������������ ��
        % ���.38.211 ����. 7.4.1.1.1
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