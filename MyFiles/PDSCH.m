function symb_PDSCH = PDSCH(codeword,modulation,N_layers,N_ID,RNTI)
    %% PDSCH ������ ���� ��������� �������
    % ?���� ������ �������� � ��� �������� �����������, �������� � ������
    % �������� �������� �� ������ ��������(����). 
    % �� ���� ������� ���� ���� ��� ��� ������� ����� PDSCH (codeword), ��� ������� �
    % ���.38.212 ����.7.2.6., ��� �������� (modulation) �� ������ ���
    % ���������� ������� ����, ���������� ������� �������� (N_laers),
    % ������������� �������������� (N_ID), � ��������� ������������� ��������� (RNTI).
    %   - codeword - ����� ���� ������ �������� (������������� ���� ������� �����)
    % ��� ������ ����� �� ������ ��� ���� ������ �������� (������������� ���� ��� ��� ������� �����).
    %   - modulation - ����� ���� ������� ��� �QPSK�, �16QAM�, �64QAM�, �256QAM�. ?��� codeword �������� ��� ������� �����, ���� 
    % ������ �������� ����� ��������� � ����� ������� ������. � �������� ������������, ������
    % ������ ��� ������ ������� ���������� �������� ����� ���� ������������ �� �������
    % ��������� ���� �������� �� ������� �������� �����.
    %   - N_laers - ���������� ������� �������� (1 ... 4 �� ������ �������� �����,
    % 5 ... 8 �� ���� ������� ����).
    %   - N_ID - ����������� ����� ���� ������������� ������
    % NCellID (0 ... 1007) ��� �������� ����� �������� ����� dataScrabingIdentityPDSCH (0 ... 1023).
    %   - RNTI - ��������� ������� � ��������� (0 ... 65535).
    
    N_codeword = 1 + (N_layers > 4);
    
        if ~iscell(codeword)
            cell_codeword = {codeword};
        else
            cell_codeword = codeword;
        end
        if ~iscell(modulation)
            modulation_cell = {modulation};
        else
            modulation_cell = modulation;
        end
    for q = 1:N_codeword
        E = length(cell_codeword{q}); 
        c_init = (double(RNTI) * 2^15) + (double(q-1) * 2^14) + double(N_ID);
        % ������� �������������� ������������������
        scramb_seq = scrambling(c_init,E);
        % ���������� � ������������ � ���. 38.211 ����. 7.3.1.1
        scrambl_out{q} = xor(cell_codeword{q},scramb_seq');
        % ������� � ������������ � ���.38.211 ����. 7.3.1.2
        %mod_symb{q} = bfModulator(double(scrambl_out{q}),modulation_cell{q}); 
        mod_symb{q} = nrSymbolModulate(double(scrambl_out{q}),modulation_cell{q}); 
    end
    
    % ������ �������� �������� �� ������ ��������(����) , ���. 38.211 ������ 7.3.1.3
    symb_PDSCH = Layer_Mapping(mod_symb,N_layers);
end