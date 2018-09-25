% example using default values, just chA (ext trigger):
% usage
% 'initialize_and_run',...
%            SampleRate,RecordLength,RecordCountPerChannel,Buffers,Repeats,...
%            A,B,T,[],'filename','templatename',[parameters],'options'


clear all

S.card                  = 'ATS9440';
S.systemID              = 1;
S.FIFO                  = 0;
S.numBuffers            = 4;

S.SampleRate            = 20e6;
S.RecordLength          = 2^7;%2^19
S.RecordCountPerChannel = 2^13;

ch.A.InputRange = 0.4; %0.4
ch.A.Coupling   = 'dc';
ch.A.InputImped = 50;

ch.B.record     = 0;
ch.B.InputRange = 0.4;
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

P.nametemplate	 = 'c:\data\template.mat';

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


S.Repeats               = 1;
loops                   = 1;


data = [];

FS1 = stoploop('stop');

ts = tic;

for i = 1:loops
        if FS1.Stop()
            FS1.Clear()
            break;
        end
        
%         filename         = ['C:\Data\1550TES\LNfun' '_' num2str(ts,'%d') '.daq'];	% auto filename
        filename         = '';	% auto filename
        % for DSO
		[test,data] = Alazar('initialize_and_run',S,ch,T,filename,P,'show','nopopup','analyze','DSO'); % for digital oscilloscope
        
        % for template creation
%          [test] = Alazar('initialize_and_run',S,ch,T,filename,P,'show','nopopup','analyze','create_template'); % for template creation
     
        % for histogram, but no photon ditribution analysis
%         [test] = Alazar('initialize_and_run',S,ch,T,filename,P,'show','nopopup','analyze','histogram','save'); % for photon histogram
        
        % for histogram, and photon ditribution analysis
        %[test] = Alazar('initialize_and_run',S,ch,T,filename,P,'show','nopopup','analyze','probability_distribution'); % % for photon probability distribution and stats
%         [test] = Alazar('initialize_and_run',S,ch,T,filename,P,'show','nopopup','analyze','probability_distribution','save'); % % for photon probability distribution and stats

        

        %for histogram, and save data
%          [test] = Alazar('initialize_and_run',S,ch,T,filename,P,'show','nopopup','analyze','histogram','save'); % for photon histogram
%        [test] = Alazar('initialize_and_run',S,ch,T,filename,P,'show','nopopup','analyze','probability_distribution','joint_histogram'); % % for photon probability distribution and stats
%       
        
        %[test] = Alazar('initialize_and_run',S,ch,T,filename,P,'show','nopopup','analyze','probability_distribution','3d'); % % for photon probability distribution and stats 3d
        %[test] = Alazar('initialize_and_run',S,ch,T,filename,P,'save','show','nopopup','analyze','probability_distribution'); % % for photon probability distribution and stats
%         [test] = Alazar('initialize_and_run',S,ch,T,filename,P,'show','nopopup','save'); % for save data
        %[test] = Alazar('initialize_and_run',S,ch,T,filename,P,'show','nopopup','analyze','save','get_variance');
        %[test] = Alazar('initialize_and_run',S,ch,T,filename,P,'show','nopopup','analyze','probability_distribution','3d','homdata'); % take data for photon probability distribution and stats 3d while reading stage position
        
end
FS1.Clear()

