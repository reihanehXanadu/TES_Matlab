% example using default values, just chA (ext trigger):
% usage
% 'initialize_and_run',...
%            SampleRate,RecordLength,RecordCountPerChannel,Buffers,Repeats,...
%            A,B,T,[],'filename','templatename',[parameters],'options'


clear all

S.card                  = 'ATS9440';
S.systemID              = 1;
S.FIFO                  = 1;
S.numBuffers            = 4;

S.SampleRate            = 20e6;
S.RecordLength          = 2^6;%2^19
S.RecordCountPerChannel = 2^16;

ch.A.InputRange = 0.2;
ch.A.Coupling   = 'dc';
ch.A.InputImped = 50;

ch.B.record     = 1;
ch.B.InputRange = 0.2;
ch.B.Coupling   = 'dc';
ch.B.InputImped = 50;

ch.C.record     = 0;
ch.C.InputRange = 1;
ch.C.Coupling   = 'dc';
ch.C.InputImped = 50;

ch.D.record     = 0;
ch.D.InputRange = 1;
ch.D.Coupling   = 'dc';
ch.D.InputImped = 50;

% filename         = 'testpulses';	% auto filename

P.nametemplate	 = 'E:\data\phaseEstimation\TES\TESEfficiencies\template.mat';

% T.TrigEng = 'K';	% J seems to be bad
T.TrigEngSource  = 'external';
T.TrigEngLevel   = 0.5;
T.TrigEngSlope   = 1;
T.TriggerTimeOut = 10;	% don't want it to timeout - wait for trigger
T.PreDepth       = 0;
T.InputRange     = 5;
T.Coupling       = 'dc';
%T.Clock          = 'fast_external';
T.Clock          = 'internal';

P.hist_bins      = 200;


S.Repeats               = 64;
loops                   = 1000;


data = [];

FS1 = stoploop('stop');

for i = 1:loops
        if FS1.Stop()
            FS1.Clear()
            break;
        end
        
        filename = ['E:\data\phaseEstimation\TES\TESEfficiencies\ppop_efficiency_devices_A215_B175_at_HOM_dip.daq'];       
       
        % for template creation
%          [test] = Alazar('initialize_and_run',S,ch,T,filename,P,'show','nopopup','analyze','create_template','save'); % for template creation
     
        % for histogram, but no photon ditribution analysis
%         [test] = Alazar('initialize_and_run',S,ch,T,filename,P,'show','nopopup','analyze','histogram','save'); % for photon histogram
        
        % for threshold histogram
%         [test] = Alazar('initialize_and_run',S,ch,T,filename,P,'show','nopopup','analyze','probability_distribution'); % % for photon probability distribution and stats
        % for joint histogram
       [test] = Alazar('initialize_and_run',S,ch,T,filename,P,'show','nopopup','analyze','probability_distribution','3d','save'); % % for photon probability distribution and stats

        
 end
FS1.Clear()

