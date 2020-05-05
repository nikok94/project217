%% 5G NR Downlink Carrier Waveform Generation
% This example implements a 5G NR downlink carrier waveform generator using
% 5G Toolbox(TM).

% Copyright 2018-2019 The MathWorks, Inc.

%% Introduction
% This example shows how to parameterize and generate a 5G New Radio (NR)
% downlink waveform. The following channels and signals are generated;
% 
% * PDSCH and its associated DM-RS
% * PDCCH and its associated DM-RS
% * PBCH and its associated DM-RS
% * PSS and SSS 
% 
% This example supports the parameterization and generation of multiple 
% SCS specific carriers and multiple bandwidth parts (BWP). Multiple 
% instances of the PDSCH and PDCCH channels can be generated over the
% different BWPs. Sets of CORESETs and search space monitoring opportunities
% can be configured for mapping the PDCCHs. Note that no precoding is
% applied to the physical channels and signals in this example.

%% Waveform and Carrier Configuration
% This section sets the SCS specific carrier bandwidths in resource blocks,
% the cell ID, and the length of the generated waveform in subframes. You
% can visualize the generated resource grids by setting the |DisplayGrids|
% field to 1. The channel bandwidth and frequency range parameters are used
% to display the associated minimum guardbands on a schematic diagram of
% the SCS carrier alignment.

waveconfig = [];
waveconfig.NCellID = 0;            % Cell identity
waveconfig.ChannelBandwidth = 40;  % Channel bandwidth (MHz)
waveconfig.FrequencyRange = 'FR1'; % 'FR1' or 'FR2'
waveconfig.NumSubframes = 10;      % Number of 1ms subframes in generated waveform (1,2,4,8 slots per 1ms subframe, depending on SCS)
waveconfig.DisplayGrids = 1;       % Display the resource grids after signal generation

% Define a set of SCS specific carriers, using the maximum sizes for a 
% 40 MHz NR channel. See TS 38.101-1 for more information on defined
% bandwidths and guardband requirements
carriers(1).SubcarrierSpacing = 15;
carriers(1).NRB = 216;
carriers(1).RBStart = 0;

carriers(2).SubcarrierSpacing = 30;
carriers(2).NRB = 106;
carriers(2).RBStart = 1;

%% SS Burst 
% In this section you can set the parameters for the SS burst. The
% numerology of the SS burst can be different from other parts of the
% waveform. This is specified via the block pattern parameter as specified
% in TS 38.213 Section 4.1. A bitmap is used to specify which blocks are
% transmitted in a 5ms half-frame burst. The periodicity in milliseconds
% and the power of the burst can also be set here. Other SS burst
% parameters not shown here can also be set. For the full list see the help
% for |hSSBurst|.

% SS burst configuration
ssburst = [];
ssburst.Enable = 1;                     % Enable SS Burst
ssburst.BlockPattern = 'Case B';        % Case B (30kHz) subcarrier spacing
ssburst.SSBTransmitted = [1 1 1 1];     % Bitmap indicating blocks transmitted in a 5ms half-frame burst
ssburst.SSBPeriodicity = 20;            % SS burst set periodicity in ms (5, 10, 20, 40, 80, 160)
ssburst.FrequencySSB = 0*5000;          % Frequency offset of SS burst (Hz), relative to waveform center (multiples of 5kHz)
ssburst.Power = 0;                      % Power scaling in dB

%% Bandwidth Parts
% A BWP is formed by a set of contiguous resources sharing a numerology on
% a given carrier. This example supports the use of multiple BWPs using a
% structure array. Each entry in the array represents a BWP. For each BWP
% you can specify the subcarrier spacing (SCS), the cyclic prefix (CP)
% length and the bandwidth. The |SubcarrierSpacing| parameter maps the BWP
% to one of the SCS specific carriers defined earlier. The |RBOffset|
% parameter controls the location of the BWP in the carrier. This is
% expressed in terms of the BWP numerology. Different BWPs can overlap with
% each other.
% 
% <<../bwp.png>>

% Bandwidth parts configurations
bwp = [];
 
bwp(1).SubcarrierSpacing = 15;          % BWP Subcarrier Spacing
bwp(1).CyclicPrefix = 'Normal';         % BWP Cyclic prefix for 15 kHz
bwp(1).NRB = 25;                        % Size of BWP
bwp(1).RBOffset = 10;                   % Position of BWP in SCS carrier

bwp(2).SubcarrierSpacing = 30;          % BWP Subcarrier Spacing
bwp(2).CyclicPrefix = 'Normal';         % BWP Cyclic prefix for 30 kHz
bwp(2).NRB = 50;                        % Size of BWP
bwp(2).RBOffset = 50;                   % Position of BWP in SCS carrier

%% CORESET and Search Space Configuration
% The parameters in this section specify the control resource set (CORESET)
% and the PDCCH search space configuration. The CORESET specifies the
% possible locations (in time and frequency) of the control channel for a
% given numerology. This example supports multiple CORESETs. The following
% parameters can be specified:
% 
% * Allocated OFDM symbols: specifies the first symbol of each CORESET
% monitoring opportunity in a slot
% * The allocated slots within a period
% * Periodicity of the allocation. If this is set to empty it indicates
% no repetition
% * CORESET duration in symbols, either 1, 2 or 3.
% * The first PRB of the allocation. Note that the allocation is in blocks
% of 6 PRBs.
% 
% Note that this example only supports non-interleaved CCE-to-REG mapped
% CORESETs.
% 
% The figure below shows the meaning of the CORESET parameters.
% 
% <<../coresetAlloc.png>>

% CORESET/search configurations
coreset = [];
coreset(1).AllocatedSymbols = [0,7];    % First symbol of each CORESET monitoring opportunity in a slot
coreset(1).AllocatedSlots = [0,1];      % Allocated slots within a period
coreset(1).AllocatedPeriod = 5;         % Allocated slot period (empty implies no repetition)
coreset(1).Duration = 3;                % CORESET symbol duration (1,2,3)
coreset(1).AllocatedPRB = 6*[0,1,3];    % 6 REG sized indices, relative to BWP (RRC - frequencyDomainResources)

%% PDCCH Instances Configuration
% This section specifies the parameters for the set of PDCCH instances in
% the waveform. Each element in the structure array defines a PDCCH
% sequence instance. The following parameters can be set:
% 
% * Enable/disable the PDCCH sequence
% * Specify the BWP carrying the PDCCH
% * PDCCH instance power in dB
% * Enable/disable DCI channel coding
% * Allocated search spaces: indices within the corresponding CORESET where
% to map the PDCCH
% * CORESET which carries the PDCCH instance
% * Periodicity of the allocation. If this is set to empty it indicates no 
% repetition
% * Number of control channel elements (CCEs) in this PDCCH
% * |NumCCE| and |StartCCE| specify the elements used for the transmission of 
% this PDCCH
% * RNTI
% * Scrambling NID for this PDCCH and its associated DM-RS
% * DM-RS power boosting
% * DCI message payload size
% * DCI message data source. You can use one of the following standard PN
% sequences: 'PN9-ITU', 'PN9', 'PN11', 'PN15', 'PN23'. The seed for the
% generator can be specified using a cell array in the form |{'PN9',seed}|.
% If no seed is specified, the generator is initialized with all ones

pdcch = [];
pdcch(1).Enable = 1;                    % Enable PDCCH sequence
pdcch(1).BWP = 1;                       % Bandwidth part
pdcch(1).Power = 1.1;                   % Power scaling in dB
pdcch(1).EnableCoding = 1;              % Enable DCI coding
pdcch(1).AllocatedSearchSpaces = [0,3]; % Index within the CORESET 
pdcch(1).CORESET = 1;                   % Control resource set ID which carries this PDCCH
pdcch(1).AllocatedPeriod = 4;           % Allocation slot period (empty implies no repetition)
pdcch(1).NumCCE = 8;                    % Number of CCE used by PDCCH {1,2,4,8,16}
pdcch(1).StartCCE = 0;                  % Starting CCE of PDCCH
pdcch(1).CCEInterleaver = 1;            % Enable interleaver 
pdcch(1).CCEREGbundleSize = 6;          % L in TS 32.211 7.3.2.2
pdcch(1).CCEInterleaverSize = 3;        % R in TS 32.211 7.3.2.2
pdcch(1).CCEnshift = 0;
pdcch(1).RNTI = 0;                      % RNTI
pdcch(1).NID = 1;                       % PDCCH and DM-RS scrambling NID 
pdcch(1).PowerDMRS = 0;                 % Additional power boosting in dB
pdcch(1).DataBlkSize = 20;              % DCI payload size
pdcch(1).DataSource = 'PN9';            % DCI data source

%% PDSCH Instances Configuration
% This section specifies the set of PDSCH instances in the waveform. Each
% element in the structure array defines a PDSCH sequence instance. This
% example defines two PDSCH sequence instances.
%
% *General Parameters* 
%
% The following parameters are set for each instance:
% 
% * Enable/disable this PDSCH sequence
% * Specify the BWP carrying the PDSCH. The PDSCH will use the SCS
% specified for this BWP
% * Power scaling in dB
% * Enable/disable DL-SCH transport channel coding
% * Transport block data source. You can use one of the following standard
% PN sequences: 'PN9-ITU', 'PN9', 'PN11', 'PN15', 'PN23'. The seed for the
% generator can be specified using a cell array in the form |{'PN9', seed}|.
% If no seed is specified, the generator is initialized with all ones
% * Target code rate used to calculate the transport block sizes
% * Overhead parameter
% * Symbol modulation
% * Number of layers
% * Redundancy version (RV) sequence

pdsch = [];
pdsch(1).Enable = 1;                    % Enable PDSCH sequence
pdsch(1).BWP = 1;                       % Bandwidth part
pdsch(1).Power = 0;                     % Power scaling in dB
pdsch(1).EnableCoding = 1;              % Enable DL-SCH transport channel coding
pdsch(1).DataSource = 'PN9';            % Channel data source 
pdsch(1).TargetCodeRate = 0.4785;       % Code rate used to calculate transport block sizes
pdsch(1).Xoh_PDSCH = 0;                 % Rate matching overhead
pdsch(1).Modulation = 'QPSK';           % 'QPSK', '16QAM', '64QAM', '256QAM'
pdsch(1).NLayers = 2;                   % Number of PDSCH layers
pdsch(1).RVSequence = [0,2,3,1];        % RV sequence to be applied cyclically across the PDSCH allocation sequence

%%
% *Allocation*
%
% The following diagram represents some of the parameters used in the PDSCH 
% allocation.
% 
% <<../pdschAlloc.png>>
% 
% You can set the following parameters to control the PDSCH allocation.
% Note that these parameters are relative to the BWP. The specified PDSCH
% allocation will avoid the locations used for the SS burst.
% 
% * Symbols in a slot allocated to each PDSCH instance
% * Slots in a frame used for the sequence of PDSCH
% * Period of the allocation in slots. If this is empty it indicates no
% repetition
% * The allocated PRBs are relative to the BWP
% * RNTI. This value is used to links the PDSCH to an instance of the PDCCH
% * NID for scrambling the PDSCH bits

pdsch(1).AllocatedSymbols = 2:10;      % Range of symbols in a slot
pdsch(1).AllocatedSlots = 0:9;         % Allocated slot indices for PDSCH sequence
pdsch(1).AllocatedPeriod = 15;         % Allocation period in slots (empty implies no repetition)
pdsch(1).AllocatedPRB = [0:5, 10:20];  % PRB allocation
pdsch(1).RNTI = 0;                     % RNTI
pdsch(1).NID = 1;                      % Scrambling for data part

%% 
% Note that the generator in this example does not check for inter-channel
% conflict. However, additional parameters can be specified for rate 
% matching around other resources, if required
% 
% * The PDSCH can be rate matched around one or more CORESETs
% * The PDSCH can be rate matched around other resource allocations 

pdsch(1).RateMatch(1).CORESET = [1];                  % Rate matching pattern, defined by CORESET IDs
pdsch(1).RateMatch(1).Pattern.AllocatedPRB = [];      % Rate matching pattern, defined by set of 'bitmaps'
pdsch(1).RateMatch(1).Pattern.AllocatedSymbols = [];
pdsch(1).RateMatch(1).Pattern.AllocatedSlots = [];
pdsch(1).RateMatch(1).Pattern.AllocatedPeriod = [];

%%
% *PDSCH DM-RS Configuration*
%
% The following DM-RS parameters can be set

% Antenna port and DM-RS configuration (TS 38.211 section 7.4.1.1)
pdsch(1).PortSet = 0:pdsch(1).NLayers-1; % DM-RS antenna ports used
pdsch(1).PDSCHMappingType = 'A';         % PDSCH mapping type ('A'(slot-wise),'B'(non slot-wise))
pdsch(1).DMRSTypeAPosition = 2;          % Mapping type A only. First DM-RS symbol position (2,3)
pdsch(1).DMRSLength = 2;                 % Number of front-loaded DM-RS symbols (1(single symbol),2(double symbol))   
pdsch(1).DMRSAdditionalPosition = 0;     % Additional DM-RS symbol positions (max range 0...3)
pdsch(1).DMRSConfigurationType = 2;      % DM-RS configuration type (1,2)
pdsch(1).NumCDMGroupsWithoutData = 0;    % CDM groups without data (max range 0...3)
pdsch(1).NIDNSCID = 1;                   % Scrambling identity (0...65535)
pdsch(1).NSCID = 0;                      % Scrambling initialization (0,1)
pdsch(1).PowerDMRS = 0;                  % Additional power boosting in dB

pdsch(1).PTRSLength = 2;                 % может принимать следующие значения (1,2,4)
pdsch(1).kPTRS = 2;                      % может принимать следующие значения (2,4)

%%
% *Specifying Multiple PDSCH Instances*
%
% A second PDSCH sequence instance is specified next using the second BWP.

pdsch(2) = pdsch(1);
pdsch(2).Enable = 1;
pdsch(2).BWP = 2;                           % PDSCH mapped to 2nd BWP
pdsch(2).AllocatedSymbols = 0:11;
pdsch(2).AllocatedSlots = [2:4,6:20];
pdsch(2).AllocatedPRB = [25:30, 35:38];     % PRB allocation, relative to BWP

%% Waveform Generation
% This section collects all the parameters into the carrier configuration
% and generates the waveform.

% Collect together channel oriented parameter sets into a single
% configuration
waveconfig.SSBurst = ssburst;
waveconfig.Carriers = carriers;
waveconfig.BWP = bwp;
waveconfig.CORESET = coreset;
waveconfig.PDCCH = pdcch;
waveconfig.PDSCH = pdsch;

% Generate complex baseband waveform
%[waveform,bwpset] = hNRDownlinkWaveformGenerator(waveconfig);

waveform = bwpToResourceElementGrids(waveconfig);

