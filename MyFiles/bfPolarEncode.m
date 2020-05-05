function out = bfPolarEncode(in,E,varargin)
%bfPolarEncode �������� ����� - ������� ������� 5GNR 
  % out = bfPolarEncode(in,E) ���������� ������� ������������ ��������
  % ������� ������ in
  % in ������ ������� ������ K, ������� ���� CRC, ���� ������� ������������
  % out ������ ������� ����� N 
  % varargin - ������� �������� ���������� �����, �������� � ����
  % �������������� ��������� nMax � iIL
  % nMax - ������������ �������� n, ������ ����� �������� �������� N = 2^n
  % ���������� � �. 7.3.3 TS 38.212
  % iIL - ����������� �����, ��������� �������� true ��� false
  % ���������� � �. 5.3.1.1 TS 38.212
  % E - ����� �������� ������������������ ��� ������������ ��������
  % �.5.4.1. TS 38.212
  
  narginchk(2,4); % ������ ���������� ������� ���������� ���� 2, ���� 4
  %% �������� �������� ��������
    if nargin==2
        % �������� ���������� ������ (Downlink) 
        % bfPolarEncode(in,E)
        nMax = 9;       % maximum n value for N
        iIL = true;     % input interleaving
    elseif nargin==3
        error('bfPolarEncode_InvalidNumInputs');
    else
        % bfPolarEncode(in,E,nMax,iIL)
        nMax = varargin{1};
        iIL = varargin{2};
    end
  checkValidIn(in,E,nMax,iIL);
  K = length(in(:,1));
  %% ������������ ������� �������� �.5.3.1.1 TS 38.212
  if iIL
    IPIndex = bfGetInterlivIndex(K); % ������� ����������� 
    inInterleaving = in(IPIndex);
  else
    inInterleaving = in;
  end
  %% ������� �������� N �.5.3.1 TS 38.212
  if ((E <= (9/8)*2^(ceil(log2(E)) - 1)) && (K/E < 9/16))
    n1 = ceil(log2(E)) - 1;
  else
    n1 = ceil(log2(E));
  end
  rMin = 1/8;
  n2 = ceil(log2(K/rMin));
  nMin = 5;
  n = max(min([n1 n2 nMax]),nMin);
  N = 2^n; 
  
  %% UCI ��� PC-Polar, � 6.3.1.3.1 TS 38.212
    if (K>=18 && K<=25) % �. 5.3.1 TS 38.212
        nPC = 3;
        if ((E - K + 3) > 192)
            nPCwm = 1;
        else
            nPCwm = 0;
        end
    else                 % CA-Polar
        nPC = 0;
        nPCwm = 0;
    end
  %% 
    qTable = polarQTab; % �. 5.3.1.2-1 TS 38.212
    idx = (qTable < N);
    qSeq = qTable(idx);                            % 0-based

  %% ������� 5.4.1.1-1
     SbIP = [ 0 1 2 4 ...
              3 5 6 7 ...
              8 16 9 17 ...
              10 18 11 19 ...
              12 20 13 21 ...
              14 22 15 23 ...
              24 25 26 28 ...
              27 29 30 31]';
      
 %% 5.4.1 Rate matching for Polar code 
 % ������� Subblock interleaving pattern ��� Polar rate-matching �. 5.4.1.1.TS 38.212
  jn = zeros(N,1);
    
  for n = 0:N-1
     i = floor(32*n/N);
     jn(n+1) = SbIP(i+1)*(N/32)+ mod(n,N/32);
  end
    
qNFtmp = [];

    if E < N
        if K/E <= 7/16  % puncturing
            for i = 0:(N-E-1)
                qNFtmp = [qNFtmp; jn(i+1)];   
            end
            if E >= 3*N/4
                uLim = ceil(3*N/4-E/2);
                qNFtmp = [qNFtmp; (0:uLim-1).'];
            else
                uLim = ceil(9*N/16-E/4);
                qNFtmp = [qNFtmp; (0:uLim-1).'];
            end
            qNFtmp = unique(qNFtmp);
        else            % shortening
            for i = E:N-1
                qNFtmp = [qNFtmp; jn(i+1)];   
            end
        end
    end
%% ������� qI �� qNFtmp � qSeq
    qI = zeros(K+nPC,1);
    j = 0;
    for i = 1:N
        ind = qSeq(N-i+1);      
        if any(ind==qNFtmp)
            continue;
        end
        j = j+1;
        qI(j) = ind;
        if j==(K+nPC)
            break;
        end
    end
%% ������� qF 
 qF = zeros(length(qSeq)-length(qI),1);
    c = 0;
    for i = 1:length(qSeq)
        tmp = qSeq(i);
        k = 0;
        for j = 1:length(qI)
            if qI(j)==tmp
                break;
            else
                k = k+1;
            end
        end
        if k == length(qI)
            c = c+1;
            qF(c) = tmp;
        end
    end
%%
  F = zeros(N,1);
  F(qF+1) = ones(length(qF),1);
%%
 qPC = zeros(nPC,1);
    if nPC > 0
        qPC(1:(nPC-nPCwm),1) = qI(end-(nPC-nPCwm)+1:end); % least reliable

        if nPCwm>0  % assumes ==1, if >0.
            G = bfGetKron(N);
            wg = sum(G,2);              % row weight

            qtildeI = qI(1:end-nPC,1);
            wt_qtildeI = wg(qtildeI+1);
            minwt = min(wt_qtildeI);    % minimum weight
            allminwtIdx = find(wt_qtildeI==minwt);

            % most reliable, minimum row weight is first value
            qPC(nPC,1) = qtildeI(allminwtIdx(1));
        end
    end
%%
 N = length(F);
    nPC = length(qPC);

    % Generate u
    u = zeros(N,1);     % doubles only
    if nPC > 0
        % Parity-Check Polar (PC-Polar)
        y0 = 0; y1 = 0; y2 = 0; y3 = 0; y4 = 0;
        k = 1;
        for idx = 1:N
            yt = y0; y0 = y1; y1 = y2; y2 = y3; y3 = y4; y4 = yt;
            if F(idx)   % frozen bits
                u(idx) = 0;
            else        % info bits
                if any(idx==(qPC+1))
                    u(idx) = y0;
                else
                    u(idx) = inInterleaving(k); % Set information bits (interleaved)
                    k = k+1;
                    y0 = double(xor(y0,u(idx)));
                end
            end
        end
    else
        % CRC-Aided Polar (CA-Polar)
        u(F==0) = inInterleaving;   % Set information bits (interleaved)
    end
%% 
    G = bfGetKron(N);

    % Encode using matrix multiplication
    outd = mod(u'*G,2)';
    out = cast(outd,class(in));
end

function checkValidIn(in,E,nMax,iIL)
  if nMax ~= 9 && nMax ~= 10
    error('bfPolarEncode_InvalidnMax');
  end
  K = length(in);
  
  if nMax==9 && iIL && (K < 36 || K > 164)
    error('bfPolarEncode_InvalidInLength');
  end
  
  if nMax==10 && ~iIL && (K<18 || (K>25 && K<31)|| K>1023)
    error('bfPolarEncode_InvalidInLength');
  end

  if (K>=18 && K<=25) % ���� ������������ PC-Polar
    nPC = 3;
  else
    nPC = 0;
  end
  
  if (E > 8192 || E < K+nPC)
    error('bfPolarEncode_InvalidELength');
  end
end

function out = polarQTab()
% ������� �. 5.3.1.2-1 
  table = [  0    518     94    214    364    414    819     966
             1     54    204    309    654    223    814     755
             2     83    298    188    659    663    439     859
             4     57    400    449    335    692    929     940
             8    521    608    217    480    835    490     830
            16    112    352    408    315    619    623     911
            32    135    325    609    221    472    671     871
             3     78    533    596    370    455    739     639
             5    289    155    551    613    796    916     888
            64    194    210    650    422    809    463     479
             9     85    305    229    425    714    843     946
             6    276    547    159    451    721    381     750
            17    522    300    420    614    837    497     969
            10     58    109    310    543    716    930     508
            18    168    184    541    235    864    821     861
           128    139    534    773    412    810    726     757
            12     99    537    610    343    606    961     970
            33     86    115    657    372    912    872     919
            65     60    167    333    775    722    492     875
            20    280    225    119    317    696    631     862
           256     89    326    600    222    377    729     758
            34    290    306    339    426    435    700     948
            24    529    772    218    453    817    443     977
            36    524    157    368    237    319    741     923
             7    196    656    652    559    621    845     972
           129    141    329    230    833    812    920     761
            66    101    110    391    804    484    382     877
           512    147    117    313    712    430    822     952
            11    176    212    450    834    838    851     495
            40    142    171    542    661    667    730     703
            68    530    776    334    808    488    498     935
           130    321    330    233    779    239    880     978
            19     31    226    555    617    378    742     883
            13    200    549    774    604    459    445     762
            48     90    538    175    433    622    471     503
            14    545    387    123    720    627    635     925
            72    292    308    658    816    437    932     878
           257    322    216    612    836    380    687     735
            21    532    416    341    347    818    903     993
           132    263    271    777    897    461    825     885
            35    149    279    220    243    496    500     939
           258    102    158    314    662    669    846     994
            26    105    337    424    454    679    745     980
           513    304    550    395    318    724    826     926
            80    296    672    673    675    841    732     764
            37    163    118    583    618    629    446     941
            25     92    332    355    898    351    962     967
            22     47    579    287    781    467    936     886
           136    267    540    183    376    438    475     831
           260    385    389    234    428    737    853     947
           264    546    173    125    665    251    867     507
            38    324    121    557    736    462    637     889
           514    208    553    660    567    442    907     984
            96    386    199    616    840    441    487     751
            67    150    784    342    625    469    695     942
            41    153    179    316    238    247    746     996
           144    165    228    241    359    683    828     971
            28    106    338    778    457    842    753     890
            69     55    312    563    399    738    854     509
            42    328    704    345    787    899    857     949
           516    536    390    452    591    670    504     973
            49    577    174    397    678    783    799    1000
            74    548    554    403    434    849    255     892
           272    113    581    207    677    820    964     950
           160    154    393    674    349    728    909     863
           520     79    283    558    245    928    719     759
           288    269    122    785    458    791    477    1008
           528    108    448    432    666    367    915     510
           192    578    353    357    620    901    638     979
           544    224    561    187    363    630    748     953
            70    166    203    236    127    685    944     763
            44    519     63    664    191    844    869     974
           131    552    340    624    782    633    491     954
            81    195    394    587    407    711    699     879
            50    270    527    780    436    253    754     981
            73    641    582    705    626    691    858     982
            15    523    556    126    571    824    478     927
           320    275    181    242    465    902    968     995
           133    580    295    565    681    686    383     765
            52    291    285    398    246    740    910     956
            23     59    232    346    707    850    815     887
           134    169    124    456    350    375    976     985
           384    560    205    358    599    444    870     997
            76    114    182    405    668    470    917     986
           137    277    643    303    790    483    727     943
            82    156    562    569    460    415    493     891
            56     87    286    244    249    485    873     998
            27    197    585    595    682    905    701     766
            97    116    299    189    573    795    931     511
            39    170    354    566    411    473    756     988
           259     61    211    676    803    634    860    1001
            84    531    401    361    789    744    499     951
           138    525    185    706    709    852    731    1002
           145    642    396    589    365    960    823     893
           261    281    344    215    440    865    922     975
            29    278    586    786    628    693    874     894
            43    526    645    647    689    797    918    1009
            98    177    593    348    374    906    502     955
           515    293    535    419    423    715    933    1004
            88    388    240    406    466    807    743    1010
           140     91    206    464    793    474    760     957
            30    584     95    680    250    636    881     983
           146    769    327    801    371    694    494     958
            71    198    564    362    481    254    702     987
           262    172    800    590    574    717    921    1012
           265    120    402    409    413    575    501     999
           161    201    356    570    603    913    876    1016
           576    336    307    788    366    798    847     767
            45     62    301    597    468    811    992     989
           100    282    417    572    655    379    447    1003
           640    143    213    219    900    697    733     990
            51    103    568    311    805    431    827    1005
           148    178    832    708    615    607    934     959
            46    294    588    598    684    489    882    1011
            75     93    186    601    710    866    937    1013
           266    644    646    651    429    723    963     895
           273    202    404    421    794    486    747    1006
           517    592    227    792    252    908    505    1014
           104    323    896    802    373    718    855    1017
           162    392    594    611    605    813    924    1018
            53    297    418    602    848    476    734     991
           193    770    302    410    690    856    829    1020
           152    107    649    231    713    839    965    1007
            77    180    771    688    632    725    938    1015
           164    151    360    653    482    698    884    1019
           768    209    539    248    806    914    506    1021
           268    284    111    369    427    752    749    1022
           274    648    331    190    904    868    945    1023 ];
       
 out = table(:);

end