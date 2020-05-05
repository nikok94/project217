function [Zc, iLS, c_out] = bfCBSegmentANDCRCAttach(dataIn, LDPCgraph)
% ������ ������� ���������� � �.5.2.2 TS 38.212
% ������� ������� ������������������ dataIn ��� ���������� ������ � ����
% ������� ������� ������ length(dataIn(:,1)) > 0. ���� ������ �������
% ������������������ ������ ��� ������������ ������ �������� ����� Kcb 
% (Kcb = 8448 ��� ����� 1 �  Kcb = 3840 ��� ����� 2)��� 
% ������� �����, �� ����������� ����������� ������� ������������������ � �
% ������� �������� ����������� CRC ����� 24
% Zc - ����������� �������� Zc �� �������� �.5.3.2-1 TS 38.212 Zlist 
% ��� ������� Kb*Zc >= Kd
% iLS - ������ �������� Zc
B = length(dataIn);

if LDPCgraph == 1
    % ��� �������� ����� 1 LDPC ������������ ������ �������� ����� ����������:
    Kcb = 8448;
else
    % ��� �������� ����� 2 LDPC ������������ ������ �������� ����� ����������:
    Kcb = 3840;
end

% ���������� ���������� ������� ������ C 

if B <= Kcb
    % ���� ������ ������� ������������������ ��� ������ ��� ������������ ������ 
    % �������� ����� Kcb, �����   
    L = 0;
    C = 1;
    Bd = B;
else
    L = 24;
    C = ceil(B/(Kcb-L));
    Bd = B + C*L;
end

% ��������� c_out ��� ����������� ������ [Kr, r], ��� 1 <= r < C + 1 (�.�.
% � ������� ������ �� 1) Kr = K, K ������������ ��������� �������:

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

% �������� ������� Zlist �������� LDPC ������������ � �.5.3.2-1
% TS 38.212
Zlist = getZsetsOfLDPC;

% ������� ����������� �������� Zc �� ������� Zlist ��� ������� Kb*Zc >= Kd
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
        % ��� ������� ����� ��������� CRC, ��������� ������� "24B"
        cCRC = bfCRCEncode(c_out(1 : Kd - L, r), '24B', 0);
        c_out(1 : Kd, r) = cCRC;
    end
    
    c_out(Kd + 1 : K, r) = -1;
  end
end