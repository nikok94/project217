function cwout = bfDLSCH(trblk, modulation, nlayers, rv, outlen, codeRate)
% out ������������ � ������ ������� �������� ���� ��� ��� ������� �����
% �������������� � ���� ������� �������, ������� ����� ����� �����
% ��������� �����
% nlayers ��������� ����� [1 : 8] TS 38.214 �. 5.1.1.1
% indata �������� 1 ��� 2 ������������ �����, ���� nlayers > 4, ��  trblk
% ������� �� 2 ������������ ������, ����� �� 1-��
% modulation ������ ���������� ���� ��� ��� ���� ��������� ��� ������������
% ������� �������� �����
% rv ������ �� 1:2 ������ �������������� ������, ������ �� ������� ����� 
% ���� �����  RV [0:3]
if nlayers > 8
    error('bfDLSCH: invalid nlayers');
end

if nlayers > 4
    if length(trblk) < 2 || length(modulation) < 2
        error('bfDLSCH: invalid input data size');
    end
    is2CW = true;
    nl1 = floor(nlayers/2);
    nl2 = ceil(nlayers/2);
else
    is2CW = false;
    nl1 = nlayers;
end

  if is2CW
    cw1 = bfGetDLSCHCodeword(trblk{1}, codeRate{1}, rv{1}, modulation{1}, nl1, outlen{1}); 
    cw2 = bfGetDLSCHCodeword(trblk{2}, codeRate{2}, rv{2}, modulation{2}, nl2, outlen{2});
    cwout = {cw1, cw2};
  else    % ���� ������� �����
    cwout = {bfGetDLSCHCodeword(trblk{1}, codeRate{1}, rv{1}, modulation{1}, nl1, outlen{1})};
  end

end