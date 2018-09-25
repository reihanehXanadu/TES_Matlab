function defs = AlazarDefs(card)
% -------------------------------------------------------------------------
% modified by NT 110503 to add card specific entries (and defs. structure)
% copied from:
% Title:   AlazarDefs.m
% Version: 6.0.0
% Date:    2011/01/24
% --------------------------------------------------------------------------

switch card
	case 'ATS460'
		defs.card = 'ATS460';
		defs.minsample = 128;
			% worked using S.RecordLength  = 2^4, S.RecordCountPerChannel = 2^3;
		defs.maxbuffersize = 2^22;	% 4 Mb
			% from software manual v5_6_0, pg 133
			% is this same as memsize?
		defs.SampleRates = [125e6,100e6,50e6,20e6,10e6,5e6,...
			2e6,1e6,500e3,200e3,100e3,50e3,20e3,10e3];
		defs.inputranges_1000000 = [0.02,0.04,0.05,0.08,0.1,...
				0.2,0.4,0.5,0.8,1,2,4,5,8,10];
		defs.inputranges_50 = [0.02,0.04,0.05,0.08,0.1,...
				0.2,0.4,0.5,0.8,1,2,4];
		defs.maxPreTrigger = -64;	% code for (RecordLength - 64)
		defs.BWLimit = 20e6;
		defs.FIFO = 0;
		defs.numchannels = 2;
	case 'ATS660'
		defs.card = 'ATS660';
		defs.minsample = 128;
		defs.maxbuffersize = 2^22;	% 4 Mb
			% from software manual v5_6_0, pg 133
		defs.SampleRates = [125e6,100e6,50e6,20e6,10e6,5e6,...
			2e6,1e6,500e3,200e3,100e3,50e3,20e3,10e3,5e3,2e3,1e3];
		defs.inputranges_1000000 = [0.2,0.4,0.8,2,4,8,16];
		defs.inputranges_50 = [0.2,0.4,0.8,2,4];
		defs.maxPreTrigger = -64;	% code for (RecordLength - 64)
		defs.BWLimit = 20e6;
		defs.FIFO = 0;
		defs.numchannels = 2;
	case 'ATS9462'
		defs.card = 'ATS9462';
		defs.minsample = 256;
		defs.maxbuffersize = 2^24;	% 16 Mb
			% from software manual v5_6_0, pg 133
		defs.SampleRates = [180e6,160e6,125e6,100e6,50e6,20e6,10e6,5e6,...
			2e6,1e6,500e3,200e3,100e3,50e3,20e3,10e3,5e3,2e3,1e3];
		defs.inputranges_1000000 = [0.2,0.4,0.8,2,4,8,16];
		defs.inputranges_50 = [0.2,0.4,0.8,2,4];
		defs.maxPreTrigger = 2048;
		defs.BWLimit = 20e6;
		defs.FIFO = 1;
		defs.numchannels = 2;
     case 'ATS9440'
		defs.card = 'ATS9440';
		defs.minsample = 128;
			% worked using S.RecordLength  = 2^4, S.RecordCountPerChannel = 2^3;
		defs.maxbuffersize = 2^24;	% 16 Mb
			% from software manual v5_6_0, pg 133
			% is this same as memsize?
		defs.SampleRates = [125e6,100e6,50e6,20e6,10e6,5e6,...
			2e6,1e6,500e3,200e3,100e3,50e3,20e3,10e3,5e3,2e3,1e3];
		defs.inputranges_1000000 = [];
		defs.inputranges_50 = [0.1,0.2,0.4,1,2,4];
		defs.maxPreTrigger = -64;	% code for (RecordLength - 64)
		defs.BWLimit = [];
		defs.FIFO = 1;
		defs.numchannels = 4;
end
%--------------------------------------------------------------------------
% Return codes
%--------------------------------------------------------------------------

defs.ApiSuccess                              = int32(512);
defs.ApiFailed								= int32(513);
defs.ApiAccessDenied                         = int32(514);
defs.ApiDmaChannelUnavailable				= int32(515);
defs.ApiDmaChannelInvalid					= int32(516);
defs.ApiDmaChannelTypeError					= int32(517);
defs.ApiDmaInProgress						= int32(518);
defs.ApiDmaDone								= int32(519);
defs.ApiDmaPaused							= int32(520);
defs.ApiDmaNotPaused                         = int32(521);
defs.ApiDmaCommandInvalid					= int32(522);
defs.ApiDmaManReady							= int32(523);
defs.ApiDmaManNotReady						= int32(524);
defs.ApiDmaInvalidChannelPriority			= int32(525);
defs.ApiDmaManCorrupted						= int32(526);
defs.ApiDmaInvalidElementIndex				= int32(527);
defs.ApiDmaNoMoreElements					= int32(528);
defs.ApiDmaSglInvalid						= int32(529);
defs.ApiDmaSglQueueFull						= int32(530);
defs.ApiNullParam							= int32(531);
defs.ApiInvalidBusIndex						= int32(532);
defs.ApiUnsupportedFunction					= int32(533);
defs.ApiInvalidPciSpace						= int32(534);
defs.ApiInvalidIopSpace						= int32(535);
defs.ApiInvalidSize							= int32(536);
defs.ApiInvalidAddress						= int32(537);
defs.ApiInvalidAccessType					= int32(538);
defs.ApiInvalidIndex                         = int32(539);
defs.ApiMuNotReady							= int32(540);
defs.ApiMuFifoEmpty							= int32(541);
defs.ApiMuFifoFull							= int32(542);
defs.ApiInvalidRegister						= int32(543);
defs.ApiDoorbellClearFailed					= int32(544);
defs.ApiInvalidUserPin						= int32(545);
defs.ApiInvalidUserState                     = int32(546);
defs.ApiEepromNotPresent                     = int32(547);
defs.ApiEepromTypeNotSupported				= int32(548);
defs.ApiEepromBlank							= int32(549);
defs.ApiConfigAccessFailed					= int32(550);
defs.ApiInvalidDeviceInfo					= int32(551);
defs.ApiNoActiveDriver						= int32(552);
defs.ApiInsufficientResources				= int32(553);
defs.ApiObjectAlreadyAllocated				= int32(554);
defs.ApiAlreadyInitialized					= int32(555);
defs.ApiNotInitialized						= int32(556);
defs.ApiBadConfigRegEndianMode				= int32(557);
defs.ApiInvalidPowerState					= int32(558);
defs.ApiPowerDown							= int32(559);
defs.ApiFlybyNotSupported					= int32(560);
defs.ApiNotSupportThisChannel				= int32(561);
defs.ApiNoAction                             = int32(562);
defs.ApiHSNotSupported						= int32(563);
defs.ApiVPDNotSupported						= int32(564);
defs.ApiVpdNotEnabled						= int32(565);
defs.ApiNoMoreCap							= int32(566);
defs.ApiInvalidOffset						= int32(567);
defs.ApiBadPinDirection						= int32(568);
defs.ApiPciTimeout							= int32(569);
defs.ApiDmaChannelClosed                     = int32(570);
defs.ApiDmaChannelError						= int32(571);
defs.ApiInvalidHandle						= int32(572);
defs.ApiBufferNotReady						= int32(573);
defs.ApiInvalidData							= int32(574);
defs.ApiDoNothing							= int32(575);
defs.ApiDmaSglBuildFailed					= int32(576);
defs.ApiPMNotSupported						= int32(577);
defs.ApiInvalidDriverVersion                 = int32(578);
defs.ApiWaitTimeout							= int32(579);
defs.ApiWaitCanceled                         = int32(580);
defs.ApiBufferTooSmall						= int32(581);
defs.ApiBufferOverflow						= int32(582);
defs.ApiInvalidBuffer						= int32(583);
defs.ApiInvalidRecordsPerBuffer				= int32(584);
defs.ApiDmaPending							= int32(585);
defs.ApiLockAndProbePagesFailed				= int32(586);
defs.ApiWaitAbandoned						= int32(587);
defs.ApiWaitFailed							= int32(588);
defs.ApiTransferComplete                     = int32(589);
defs.ApiPllNotLocked                         = int32(590);
defs.ApiNotSupportedInDualChannelMode        = int32(591);

%--------------------------------------------------------------------------
% Board types
%--------------------------------------------------------------------------

defs.ATS_NONE        = int32(0);
defs.ATS850          = int32(1);
defs.ATS310          = int32(2);
defs.ATS330          = int32(3);
defs.ATS855          = int32(4);
defs.ATS315          = int32(5);
defs.ATS335          = int32(6);
defs.ATS460          = int32(7);
defs.ATS860          = int32(8);
defs.ATS660          = int32(9);
defs.ATS665          = int32(10);
defs.ATS9462         = int32(11);
defs.ATS9434         = int32(12);
defs.ATS9870         = int32(13);
defs.ATS9350         = int32(14);
defs.ATS9325         = int32(15);
defs.ATS9440         = int32(16);
defs.ATS_LAST        = int32(17);

%--------------------------------------------------------------------------
% Clock Control
%--------------------------------------------------------------------------

% Clock sources
defs.INTERNAL_CLOCK              =	hex2dec('00000001');
defs.EXTERNAL_CLOCK              =	hex2dec('00000002');
defs.FAST_EXTERNAL_CLOCK         =	hex2dec('00000002')';
defs.MEDIMUM_EXTERNAL_CLOCK      =	hex2dec('00000003')';
defs.MEDIUM_EXTERNAL_CLOCK       =	hex2dec('00000003')';
defs.SLOW_EXTERNAL_CLOCK         =	hex2dec('00000004')';
defs.EXTERNAL_CLOCK_AC           =	hex2dec('00000005')';
defs.EXTERNAL_CLOCK_DC           =	hex2dec('00000006')';
defs.EXTERNAL_CLOCK_10MHz_REF    =	hex2dec('00000007')';
defs.INTERNAL_CLOCK_DIV_5        =	hex2dec('000000010')';
defs.MASTER_CLOCK                =	hex2dec('000000011')';

% Internal sample rates
defs.SAMPLE_RATE_1KSPS           =	hex2dec('00000001');
defs.SAMPLE_RATE_2KSPS           =	hex2dec('00000002');
defs.SAMPLE_RATE_5KSPS           =	hex2dec('00000004');
defs.SAMPLE_RATE_10KSPS          =	hex2dec('00000008');
defs.SAMPLE_RATE_20KSPS          =	hex2dec('0000000A');
defs.SAMPLE_RATE_50KSPS          =	hex2dec('0000000C');
defs.SAMPLE_RATE_100KSPS         =	hex2dec('0000000E');
defs.SAMPLE_RATE_200KSPS         =	hex2dec('00000010');
defs.SAMPLE_RATE_500KSPS         =	hex2dec('00000012');
defs.SAMPLE_RATE_1MSPS           =	hex2dec('00000014');
defs.SAMPLE_RATE_2MSPS           =	hex2dec('00000018');
defs.SAMPLE_RATE_5MSPS           =	hex2dec('0000001A');
defs.SAMPLE_RATE_10MSPS          =	hex2dec('0000001C');
defs.SAMPLE_RATE_20MSPS          =	hex2dec('0000001E');
defs.SAMPLE_RATE_25MSPS          =	hex2dec('00000021');
defs.SAMPLE_RATE_50MSPS          =	hex2dec('00000022');
defs.SAMPLE_RATE_100MSPS         =	hex2dec('00000024');
defs.SAMPLE_RATE_125MSPS         =   hex2dec('00000025');
defs.SAMPLE_RATE_160MSPS         =   hex2dec('00000026');
defs.SAMPLE_RATE_180MSPS         =   hex2dec('00000027');
defs.SAMPLE_RATE_200MSPS         =	hex2dec('00000028');
defs.SAMPLE_RATE_250MSPS         =   hex2dec('0000002B');
defs.SAMPLE_RATE_500MSPS         =   hex2dec('00000030');
defs.SAMPLE_RATE_1GSPS           =   hex2dec('00000035');
defs.SAMPLE_RATE_2GSPS           =   hex2dec('0000003A');
defs.SAMPLE_RATE_USER_DEF        =	hex2dec('00000040'); 

% Clock edges
defs.CLOCK_EDGE_RISING           =	hex2dec('00000000');
defs.CLOCK_EDGE_FALLING          =	hex2dec('00000001');

% Decimation
defs.DECIMATE_BY_8               =   hex2dec('00000008');
defs.DECIMATE_BY_64              =   hex2dec('00000040');

%--------------------------------------------------------------------------
% Input Control
%--------------------------------------------------------------------------

% Input channels
defs.CHANNEL_ALL                 =   hex2dec('00000000');
defs.CHANNEL_A                   =   hex2dec('00000001');
defs.CHANNEL_B                   =   hex2dec('00000002');
defs.CHANNEL_C                   =   hex2dec('00000004');
defs.CHANNEL_D                   =   hex2dec('00000008');
defs.CHANNEL_E                   =   hex2dec('00000010');
defs.CHANNEL_F                   =   hex2dec('00000012');
defs.CHANNEL_G                   =   hex2dec('00000014');
defs.CHANNEL_H                   =   hex2dec('00000018');

% Input ranges
defs.INPUT_RANGE_PM_20_MV        =   hex2dec('00000001');
defs.INPUT_RANGE_PM_40_MV        =   hex2dec('00000002');
defs.INPUT_RANGE_PM_50_MV        =   hex2dec('00000003');
defs.INPUT_RANGE_PM_80_MV        =   hex2dec('00000004');
defs.INPUT_RANGE_PM_100_MV       =   hex2dec('00000005');
defs.INPUT_RANGE_PM_200_MV       =   hex2dec('00000006');
defs.INPUT_RANGE_PM_400_MV       =   hex2dec('00000007');
defs.INPUT_RANGE_PM_500_MV       =   hex2dec('00000008');
defs.INPUT_RANGE_PM_800_MV       =   hex2dec('00000009');
defs.INPUT_RANGE_PM_1_V          =   hex2dec('0000000A');
defs.INPUT_RANGE_PM_2_V          = 	hex2dec('0000000B');
defs.INPUT_RANGE_PM_4_V          =	hex2dec('0000000C');
defs.INPUT_RANGE_PM_5_V          =	hex2dec('0000000D');
defs.INPUT_RANGE_PM_8_V          =	hex2dec('0000000E');
defs.INPUT_RANGE_PM_10_V         =	hex2dec('0000000F');
defs.INPUT_RANGE_PM_20_V         =	hex2dec('00000010');
defs.INPUT_RANGE_PM_40_V         =	hex2dec('00000011');
defs.INPUT_RANGE_PM_16_V         =   hex2dec('00000012');
defs.INPUT_RANGE_HIFI            = 	hex2dec('00000020');

% Input impedances
defs.IMPEDANCE_1M_OHM            =	hex2dec('00000001');
defs.IMPEDANCE_50_OHM            =	hex2dec('00000002');
defs.IMPEDANCE_75_OHM            =	hex2dec('00000004');
defs.IMPEDANCE_300_OHM           =	hex2dec('00000008');
defs.IMPEDANCE_600_OHM           =	hex2dec('0000000A');

% Input coupling 
defs.AC_COUPLING                 =   hex2dec('00000001');
defs.DC_COUPLING                 =	hex2dec('00000002');

%--------------------------------------------------------------------------
% Trigger Control
%--------------------------------------------------------------------------

% Trigger engines
defs.TRIG_ENGINE_J               =	hex2dec('00000000');
defs.TRIG_ENGINE_K               =	hex2dec('00000001');

% Trigger engine operations
defs.TRIG_ENGINE_OP_J            =   hex2dec('00000000');
defs.TRIG_ENGINE_OP_K            =	hex2dec('00000001');
defs.TRIG_ENGINE_OP_J_OR_K		=   hex2dec('00000002');
defs.TRIG_ENGINE_OP_J_AND_K		=   hex2dec('00000003');
defs.TRIG_ENGINE_OP_J_XOR_K		=   hex2dec('00000004');
defs.TRIG_ENGINE_OP_J_AND_NOT_K	=   hex2dec('00000005');
defs.TRIG_ENGINE_OP_NOT_J_AND_K	=   hex2dec('00000006');

% Trigger engine sources
defs.TRIG_CHAN_A                 =   hex2dec('00000000');
defs.TRIG_CHAN_B                 =   hex2dec('00000001');
defs.TRIG_EXTERNAL               =   hex2dec('00000002');
defs.TRIG_DISABLE                =   hex2dec('00000003');
defs.TRIG_CHAN_C                 =   hex2dec('00000004');
defs.TRIG_CHAN_D                 =   hex2dec('00000005');

% Trigger slopes
defs.TRIGGER_SLOPE_POSITIVE      =   hex2dec('00000001');
defs.TRIGGER_SLOPE_NEGATIVE      =   hex2dec('00000002');

% External trigger ranges
defs.ETR_DIV5                    =   hex2dec('00000000');
defs.ETR_X1                      =   hex2dec('00000001');
defs.ETR_5V                      =   hex2dec('00000000');
defs.ETR_1V                      =   hex2dec('00000001');

%--------------------------------------------------------------------------
% Auxiliary I/O and LED Control
%--------------------------------------------------------------------------

% AUX outputs
defs.AUX_OUT_TRIGGER             =	0;
defs.AUX_OUT_PACER               =	2;
defs.AUX_OUT_BUSY                =	4;
defs.AUX_OUT_CLOCK               =	6;
defs.AUX_OUT_RESERVED            =	8;
defs.AUX_OUT_CAPTURE_ALMOST_DONE	=	10;
defs.AUX_OUT_AUXILIARY			=	12;
defs.AUX_OUT_SERIAL_DATA			=	14;
defs.AUX_OUT_TRIGGER_ENABLE		=	16;

% AUX inputs
defs.AUX_IN_TRIGGER_ENABLE		=	1;
defs.AUX_IN_DIGITAL_TRIGGER		=	3;
defs.AUX_IN_GATE					=	5;
defs.AUX_IN_CAPTURE_ON_DEMAND	=	7;
defs.AUX_IN_RESET_TIMESTAMP		=	9;
defs.AUX_IN_SLOW_EXTERNAL_CLOCK	=	11;
defs.AUX_INPUT_AUXILIARY			=	13;
defs.AUX_INPUT_SERIAL_DATA		=	15;

% LED states
defs.LED_OFF                     =	hex2dec('00000000');
defs.LED_ON                      =	hex2dec('00000001');

%--------------------------------------------------------------------------
% Get/Set Parameters
%--------------------------------------------------------------------------

defs.NUMBER_OF_RECORDS           =   hex2dec('10000001');
defs.PRETRIGGER_AMOUNT           =   hex2dec('10000002');
defs.RECORD_LENGTH               =   hex2dec('10000003');
defs.TRIGGER_ENGINE              =   hex2dec('10000004');
defs.TRIGGER_DELAY               =   hex2dec('10000005');
defs.TRIGGER_TIMEOUT             =   hex2dec('10000006');
defs.SAMPLE_RATE                 =   hex2dec('10000007');
defs.CONFIGURATION_MODE          =   hex2dec('10000008'); 
defs.DATA_WIDTH                  =   hex2dec('10000009'); 
defs.SAMPLE_SIZE                 =   defs.DATA_WIDTH;
defs.AUTO_CALIBRATE              =   hex2dec('1000000A');
defs.TRIGGER_XXXXX               =   hex2dec('1000000B');
defs.CLOCK_SOURCE                =   hex2dec('1000000C');
defs.CLOCK_SLOPE                 =   hex2dec('1000000D');
defs.IMPEDANCE                   =   hex2dec('1000000E');
defs.INPUT_RANGE                 =   hex2dec('1000000F');
defs.COUPLING                    =   hex2dec('10000010');
defs.MAX_TIMEOUTS_ALLOWED        =   hex2dec('10000011');
defs.ATS_OPERATING_MODE          =   hex2dec('10000012'); 
defs.CLOCK_DECIMATION_EXTERNAL   =   hex2dec('10000013');
defs.LED_CONTROL                 =   hex2dec('10000014');
defs.ATTENUATOR_RELAY            =   hex2dec('10000018');
defs.EXT_TRIGGER_COUPLING        =   hex2dec('1000001A');
defs.EXT_TRIGGER_ATTENUATOR_RELAY    =  hex2dec('1000001C');
defs.TRIGGER_ENGINE_SOURCE       =   hex2dec('1000001E');
defs.TRIGGER_ENGINE_SLOPE        =   hex2dec('10000020');
defs.SEND_DAC_VALUE              =   hex2dec('10000021');
defs.SLEEP_DEVICE                =   hex2dec('10000022');
defs.GET_DAC_VALUE               =   hex2dec('10000023');
defs.GET_SERIAL_NUMBER           =   hex2dec('10000024');
defs.GET_FIRST_CAL_DATE          =   hex2dec('10000025');
defs.GET_LATEST_CAL_DATE         =   hex2dec('10000026');
defs.GET_LATEST_TEST_DATE        =   hex2dec('10000027');
defs.GET_LATEST_CAL_DATE_MONTH   =   hex2dec('1000002D');
defs.GET_LATEST_CAL_DATE_DAY     =   hex2dec('1000002E');
defs.GET_LATEST_CAL_DATE_YEAR    =   hex2dec('1000002F');
defs.GET_PCIE_LINK_SPEED         =   hex2dec('10000030');
defs.GET_PCIE_LINK_WIDTH         =   hex2dec('10000031');
defs.SETGET_ASYNC_BUFFCOUNT      =   hex2dec('10000040');
defs.SET_DATA_FORMAT             =   hex2dec('10000041');
defs.GET_DATA_FORMAT             =   hex2dec('10000042');
defs.DATA_FORMAT_UNSIGNED        = 0;
defs.DATA_FORMAT_SIGNED          = 1;
defs.SET_SINGLE_CHANNEL_MODE     =   hex2dec('10000043');
defs.MEMORY_SIZE                 =   hex2dec('1000002A');
defs.BOARD_TYPE                  =   hex2dec('1000002B');
defs.ASOPC_TYPE                  =   hex2dec('1000002C');
defs.GET_BOARD_OPTIONS_LOW       =   hex2dec('10000037');
defs.GET_BOARD_OPTIONS_HIGH      =   hex2dec('10000038');
defs.OPTION_STREAMING_DMA        =   uint32(2^0);
defs.OPTION_AVERAGE_INPUT        =   uint32(2^1);
defs.OPTION_EXTERNAL_CLOCK       =   uint32(2^1);
defs.OPTION_DUAL_PORT_MEMORY 	=   uint32(2^2);
defs.OPTION_180MHZ_OSCILLATOR    =   uint32(2^3);
defs.OPTION_LVTTL_EXT_CLOCK      =   uint32(2^4);
defs.OPTION_OEM_FPGA             =   uint32(2^47);
defs.TRANSFER_OFFET              =   hex2dec('10000030');
defs.TRANSFER_LENGTH             =   hex2dec('10000031');
defs.TRANSFER_RECORD_OFFSET      =   hex2dec('10000032');
defs.TRANSFER_NUM_OF_RECORDS     =   hex2dec('10000033');
defs.TRANSFER_MAPPING_RATIO      =   hex2dec('10000034');
defs.TRIGGER_ADDRESS_AND_TIMESTAMP = hex2dec('10000035');
defs.MASTER_SLAVE_INDEPENDENT    =   hex2dec('10000036');
defs.TRIGGERED                   =   hex2dec('10000040');
defs.BUSY                        =   hex2dec('10000041');
defs.WHO_TRIGGERED               =   hex2dec('10000042');
defs.SET_DATA_FORMAT				=   hex2dec('10000041');
defs.GET_DATA_FORMAT				=   hex2dec('10000042');
defs.DATA_FORMAT_UNSIGNED		=   0;
defs.DATA_FORMAT_SIGNED			=   1;
defs.SET_SINGLE_CHANNEL_MODE		=   hex2dec('10000043');
defs.GET_SAMPLES_PER_TIMESTAMP_CLOCK	=   hex2dec('10000044');
defs.GET_RECORDS_CAPTURED		=   hex2dec('10000045');
defs.GET_MAX_PRETRIGGER_SAMPLES	=   hex2dec('10000046');
defs.SET_ADC_MODE				=   hex2dec('10000047');
defs.ECC_MODE					=   hex2dec('10000048');
defs.ECC_DISABLE					=   0;
defs.ECC_ENABLE					=   1;
defs.GET_AUX_INPUT_LEVEL			=   hex2dec('10000049');
defs.AUX_INPUT_LOW				=   0;
defs.AUX_INPUT_HIGH				=   1;
defs.GET_ASYNC_BUFFERS_PENDING   =   hex2dec('10000050');
defs.GET_ASYNC_BUFFERS_PENDING_FULL =    hex2dec('10000051');
defs.GET_ASYNC_BUFFERS_PENDING_EMPTY =   hex2dec('10000052');
defs.ACF_SAMPLES_PER_RECORD      =   hex2dec('10000060');
defs.ACF_RECORDS_TO_AVERAGE      =   hex2dec('10000061');
defs.EXT_TRIGGER_IMPEDANCE		=   hex2dec('10000065');
defs.EXT_TRIG_50_OHMS			= 0;
defs.EXT_TRIG_300_OHMS			= 1;

% Master/Slave Configuration
defs.BOARD_IS_INDEPENDENT        =   hex2dec('00000000');
defs.BOARD_IS_MASTER             =	hex2dec('00000001');
defs.BOARD_IS_SLAVE              =	hex2dec('00000002');
defs.BOARD_IS_LAST_SLAVE         =	hex2dec('00000003');

% Attenuator Relay
defs.AR_X1                       =   hex2dec('00000000');
defs.AR_DIV40                    =   hex2dec('00000001');

% Device Sleep state
defs.POWER_OFF                   =   hex2dec('00000000');
defs.POWER_ON                    =   hex2dec('00000001');

% Software Events control
defs.SW_EVENTS_OFF               =   hex2dec('00000000');
defs.SW_EVENTS_ON                =   hex2dec('00000001');

% TimeStamp Value Reset Control
defs.TIMESTAMP_RESET_FIRSTTIME_ONLY	= hex2dec('00000000');
defs.TIMESTAMP_RESET_ALWAYS			= hex2dec('00000001');

% DAC Names used by API AlazarDACSettingAdjust 
defs.ATS460_DAC_A_GAIN			=   hex2dec('00000001');
defs.ATS460_DAC_A_OFFSET			=   hex2dec('00000002');
defs.ATS460_DAC_A_POSITION		=   hex2dec('00000003');
defs.ATS460_DAC_B_GAIN			=   hex2dec('00000009');
defs.ATS460_DAC_B_OFFSET			=   hex2dec('0000000A');
defs.ATS460_DAC_B_POSITION		=   hex2dec('0000000B');
defs.ATS460_DAC_EXTERNAL_CLK_REF	=   hex2dec('00000007');

% DAC Names Specific to the ATS660
defs.ATS660_DAC_A_GAIN			=   hex2dec('00000001');
defs.ATS660_DAC_A_OFFSET			=   hex2dec('00000002');
defs.ATS660_DAC_A_POSITION		=   hex2dec('00000003');
defs.ATS660_DAC_B_GAIN			=   hex2dec('00000009');
defs.ATS660_DAC_B_OFFSET			=   hex2dec('0000000A');
defs.ATS660_DAC_B_POSITION		=   hex2dec('0000000B');
defs.ATS660_DAC_EXTERNAL_CLK_REF	=   hex2dec('00000007');

% DAC Names Specific to the ATS665
defs.ATS665_DAC_A_GAIN			=   hex2dec('00000001');
defs.ATS665_DAC_A_OFFSET			=   hex2dec('00000002');
defs.ATS665_DAC_A_POSITION		=   hex2dec('00000003');
defs.ATS665_DAC_B_GAIN			=   hex2dec('00000009');
defs.ATS665_DAC_B_OFFSET			=   hex2dec('0000000A');
defs.ATS665_DAC_B_POSITION		=   hex2dec('0000000B');
defs.ATS665_DAC_EXTERNAL_CLK_REF	=   hex2dec('00000007');

% Error return values
defs.SETDAC_INVALID_SETGET       = 660;
defs.SETDAC_INVALID_CHANNEL      = 661;
defs.SETDAC_INVALID_DACNAME      = 662;
defs.SETDAC_INVALID_COUPLING     = 663;
defs.SETDAC_INVALID_RANGE        = 664;
defs.SETDAC_INVALID_IMPEDANCE    = 665;
defs.SETDAC_BAD_GET_PTR          = 667;
defs.SETDAC_INVALID_BOARDTYPE    = 668;

% Constants to be used in the Application when dealing with Custom FPGAs
defs.FPGA_GETFIRST               =   hex2dec('FFFFFFFF');
defs.FPGA_GETNEXT                =   hex2dec('FFFFFFFE');
defs.FPGA_GETLAST                =   hex2dec('FFFFFFFC');

%--------------------------------------------------------------------------
% AutoDMA Control 
%--------------------------------------------------------------------------

% AutoDMA flags 
defs.ADMA_EXTERNAL_STARTCAPTURE  =   hex2dec('00000001');
defs.ADMA_ENABLE_RECORD_HEADERS  =   hex2dec('00000008');
defs.ADMA_SINGLE_DMA_CHANNEL     =   hex2dec('00000010');
defs.ADMA_ALLOC_BUFFERS          =   hex2dec('00000020');
defs.ADMA_TRADITIONAL_MODE       =   hex2dec('00000000');
defs.ADMA_CONTINUOUS_MODE        =   hex2dec('00000100');
defs.ADMA_NPT                    =   hex2dec('00000200');
defs.ADMA_TRIGGERED_STREAMING    =   hex2dec('00000400');
defs.ADMA_FIFO_ONLY_STREAMING    =   hex2dec('00000800');
defs.ADMA_INTERLEAVE_SAMPLES     =   hex2dec('00001000');
defs.ADMA_GET_PROCESSED_DATA     =   hex2dec('00002000');

% AutoDMA header constants
defs.ADMA_CLOCKSOURCE            =   hex2dec('00000001');
defs.ADMA_CLOCKEDGE              =   hex2dec('00000002');
defs.ADMA_SAMPLERATE             =   hex2dec('00000003');
defs.ADMA_INPUTRANGE             =   hex2dec('00000004');
defs.ADMA_INPUTCOUPLING          =   hex2dec('00000005');
defs.ADMA_IMPUTIMPEDENCE         =   hex2dec('00000006');
defs.ADMA_EXTTRIGGERED           =   hex2dec('00000007');
defs.ADMA_CHA_TRIGGERED          =   hex2dec('00000008');
defs.ADMA_CHB_TRIGGERED          =   hex2dec('00000009');
defs.ADMA_TIMEOUT                =   hex2dec('0000000A');
defs.ADMA_THISCHANTRIGGERED      =   hex2dec('0000000B');
defs.ADMA_SERIALNUMBER           =   hex2dec('0000000C');
defs.ADMA_SYSTEMNUMBER           =   hex2dec('0000000D');
defs.ADMA_BOARDNUMBER            =   hex2dec('0000000E');
defs.ADMA_WHICHCHANNEL           =   hex2dec('0000000F');
defs.ADMA_SAMPLERESOLUTION       =   hex2dec('00000010');
defs.ADMA_DATAFORMAT             =   hex2dec('00000011');
