
S.card                  = 'ATS9440';
S.systemID              = 1;
S.FIFO                  = 1;
S.SampleRate            = 1e8;
S.RecordLength          = 100;
S.RecordCountPerChannel = 1;
S.numBuffers            = 1;
S.RecordsPerBuffer      = 80000;
S.numchannels2record    =4;
S.buffersPerAcquisition =10;


ch.A.InputRange = 0.4;
ch.A.Coupling   = 'dc';
ch.A.InputImped = 50;
ch.B.InputRange = 0.2;
ch.B.Coupling   = 'dc';
ch.B.InputImped = 50;
ch.C.InputRange = 2;
ch.C.Coupling   = 'dc';
ch.C.InputImped = 50;
ch.D.InputRange = 2;
ch.D.Coupling   = 'dc';
ch.D.InputImped = 1;

filename = 'C:\Data\TES_Data\Sep21\Test_pulses1'; 

 %filename = 'C:\Data\TES_Data\Sep21\pump_44dB'; 
 %filename = 'C:\Data\TES_Data\Sep20\Chip_pass_DCDC_3900nW_pumpATTp5dB_s_i';
%filename = 'power_28dB'; 
filename=[filename,'_',datestr(now, 'yyyymmdd-HHMM')];
%[test,data] = Alazar(S,filename,'calibrate','save','Analyze_realtime');
[test,data] = Alazar(S,filename,'save');
%[test,data] = Alazar(S,filename,'save','Analyze_realtime');
%[test,data] = Alazar(S,filename,'save','calibrate');
 
%  [test,data3D] = Alazar(S);