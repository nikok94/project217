function [REG, f] = bfPDCCHinterleaver(N_REG_CORESET, L, R, nshift, varargin)
% �.7.3.2.2 TS 38.211
% ������� bfPDCCHinterleaver ���������� ������� f ��� CCE 
% � ������� REG(resource-element groups) 
%(CCE ������� �� 6
% REG(resource-element groups) ������ �� ������� ������� �� 12 ������ 
% ��������� (������ ������� - ���� ������ �� ����� ��������)
% N_REG_CORESET ���������� ������ ������ ��� �������� CORESET
% ����������� CCE � REG ����� ����������� ��� ��������������, ��� � ���
% interleaver = 1 ���������� ������� �������������, ������������� ��������������
% ������� REG bundle L - ������ ���������, �������� �������� ������, 
% ������ ����� REG bundle ��� ����� ���������� ������� ������ ������ REG
% ����� �������������� ��� ���������� CCE 
% R ������ ����������� R = {0,3,6}, �������� ��������� �������� ������
% nshift = N_ID_cell ��� ������������ CORESET ��� PBCH ��� SIB1
% nshift?{0,1,...,274} � ���� ������

% REG ������������ ����� ��������� ������ ����������� [6, N_REG_CORESET/6]
% � ������ ������� REG ��� ������� �������� CCE, ����� N_REG_CORESET/6 C�E 
REG = zeros(6, N_REG_CORESET/6);

narginchk(4,5);
if nargin==4
   interleaver = 0;
else
   interleaver = varargin{1};
end

if interleaver 
  % ����������� �������������� � ������������  
  C = N_REG_CORESET/(L*R);
  f = zeros(C*R, 1);
  for r = 0 : R - 1
      for c = 0 : C - 1
         j = c*R + r + 1;
         f(j,1) = mod((r*C + c + nshift),(N_REG_CORESET/L));
      end
  end
  
  for j = 0 : N_REG_CORESET/6 - 1
    for k = 0 : 6/L - 1
        for l = 0 : L - 1
            REG(k*L + l + 1, j + 1) = f(6*j/L + k + 1)*6 + l;
        end
    end
  end
else
  % ����������� �������������� ��� �����������    
  f = zeros(N_REG_CORESET/6, 1);
  for j = 0 : N_REG_CORESET/6 - 1
     f(j + 1, 1) = j;
     for l = 0 : 5
         REG(l + 1, j + 1) = j * 6 + l;
     end
  end 
end

end