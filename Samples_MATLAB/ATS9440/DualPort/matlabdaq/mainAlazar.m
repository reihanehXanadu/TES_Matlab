% example using default values, just chA (ext trigger):
% usage
% 'initialize_and_run',...
%            SampleRate,RecordLength,RecordCountPerChannel,Buffers,Repeats,...
%            A,B,T,[],'filename','templatename',[parameters],'options'


clear all

S.card                  = 'ATS9440';
S.systemID              = 1;
S.FIFO                  = 1;
S.SampleRate            = 20e6;
S.RecordLength          = 2^8;
S.RecordCountPerChannel = 2^14;
S.numBuffers            = 4;
S.Repeats               = 2;

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
ch.D.InputImped = 50;

% filename = '';	% auto filename
filename = 'testpulses';	% auto filename

% T.TrigEng = 'K';	% J seems to be bad
T.TrigEngSource  = 'external';
T.TrigEngLevel   = .1;
T.TrigEngSlope   = 1;
T.TriggerTimeOut = 1;	% don't want it to timeout - wait for trigger
T.PreDepth       = 0;
T.InputRange     = 5;
T.Coupling       = 'dc';
%T.Clock          = 'fast_external';
T.Clock          = 'internal';

P.hist_bins      = 50;
P.nametemplate	 = '';
data = [];

% DSO acquisition

% FS1 = stoploop('stop');
% while ~FS1.Stop()
%		[test] = Alazar('initialize_and_run',S,ch,T,filename,P,'show','nopopup','analyze','DSO');
% end


% template creation

% FS1 = stoploop('stop');
% while ~FS1.Stop()
%		[test] = Alazar('initialize_and_run',S,ch,T,filename,P,'show','nopopup','analyze','create_template');
% end


% histogram creation

% FS1 = stoploop('stop');
% while ~FS1.Stop()
%		[test] = Alazar('initialize_and_run',S,ch,T,filename,P,'show','nopopup','analyze','histogram');
% end

% save data

%filename = 'E:\2011\October\Oct21\signal_bck_K_1550_50kHz_reprate_multimode_fiber';


% FS1 = stoploop('stop');
% while ~FS1.Stop()
		%[test] = Alazar('initialize_and_run',S,ch,T,filename,P,'show','nopopup','analyze','save');
         [test,data3D] = Alazar('initialize_and_run',S,ch,T,filename,P,'show','nopopup','analyze','dso');
% end

