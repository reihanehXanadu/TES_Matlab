function [varargout,data] = Alazar(op,S,ch,T,varargin)
% 
% v0 Nathan stole python code from Sae Woo & Brice, and matlab from Dan
%		Schmidt
% v1 Thomas did big rewrite (110404)
% v2 Nathan slight changes
%		no setting of maxbuffersize, onboardmem, or buffersposted...
% v3 Thomas checked channel values (sum) - should work A = 1; B = 2; C = 4;
%                                                      D = 8
%    added 9440 option. def file needs to be made
%    checked and changed before async read parameters
%    added min sample size check - need to add in def file
%    improved analysis to work with fifo and no fifo
%    added a few more input parameter checks
%    stream to disk does not seem to work - maybe we need to buy it?
% v4 Nathan - changes to input checking
%	added minBitsPerSample to defs files
%	changed to have realtime & save
%	analyze
%		change name to analyze_realtime
%		auto find pretrigger 
%		normalized template
%		got rid of for loop in template (can't test fifo...)
%		add draggable lines to histogram
%	'analyze' is flag that changes check inputs
%		if not doing realtime, way less limits
%	added 'err' to EVERYTHING! it's ugly, but crashes are handled well now
% v5 Nathan added to SVN - not keeping track of changes here anymore...
% 
% to do:
% 


% defaults
warning off all
graph       = 0;
dosave      = 0;
show        = 0;	% don't print status messages (not verbose)
analyzeData = 0;
% FIFO    = 0;	% only valid for 9462 and 9440 
return_traces  = 0;	% don't return data in varargout
append         = 0;	% new file, don't append to current
err = 0;
popup = 1;	% popup to stop acquisition

% don't know why this is needed, but
% varargin gets changed when passed thru function
if numel(varargin) == 1 && numel(varargin{1}) > 1
	varargin = varargin{1};
end
numarg = numel(varargin);

% check for options
for m = 1:numarg
	if ischar(varargin{m})	% varargin commands are always char
		switch lower(varargin{m})
            case 'analyze'
                analyzeData = 1;
            case 'graph'
				graph = 1;
			case 'save'
				dosave = 1;
			case 'show'
				show = 1;
			case 'return_traces'
				return_traces = 1;
			case 'append'
				append = 1;
			case 'nopopup'
				popup = 0;
            case 'dispData'
                dispdata = 1;
		end
	end
end


if ~exist('S','var') || isempty(S) || ~isfield(S,'card')
	S.card = 'ATS460';	% default card
else
	switch lower(S.card)	% make sure input correctly
		case {'ats9462','9462'}
			S.card = 'ATS9462';
		case {'ats9440','9440'}
			S.card = 'ATS9440';
		case {'ats660','660'}
			S.card = 'ATS660';
		case {'ats460','460'}
			S.card = 'ATS460';
        case {'ats9440','9440'}
			S.card = 'ATS9440';
		otherwise
			error('S.card has wrong value')
	end
end

if ~isfield(S,'systemID')
	S.systemID = 1;		% default if only one board, set to 2 for 2nd board
end

if ~isfield(S,'boardID')
	S.boardID = 1;	% default, set to something else for slave maybe?
end

% load correct definition file for card
% calls an mfile with definitions in 'defs' structure
defs = AlazarDefs(S.card);


% check whether device initialized
% just setup to talk to one board at a time as written
global OPEN_Alazar
global dev_Alazar
lib = 'ATSApi';	% dll library
comp = computer;	% info about whether 32 bit or 64 bit
if isempty(OPEN_Alazar) || OPEN_Alazar ~= 1, 
	% Load driver library 
	if ~alazarLoadLibrary()
		disp('Error: ATSApi.dll not loaded');
		return
	end
	% libfunctionsview(lib)	% list of dll functions
    
% 	% Find the number of boards in this board system
% 	boardCount = calllib(lib, 'AlazarBoardsInSystemBySystemID', systemID);
% 	if boardCount < 1
% 		fprintf('Error: No boards found in system\n');
% 		return
% 	end

	% Get a handle to the board
	dev = calllib(lib, 'AlazarGetBoardBySystemID', S.systemID, S.boardID);
	setdatatype(dev, 'voidPtr', 1, 1);
	if dev.Value == 0
		fprintf('Error: Unable to open board system ID %u board ID %u\n', S.systemID, S.boardID);
		return
	end

	% Get the board type in this board system
	boardTypeId = calllib('ATSApi', 'AlazarGetBoardKind', dev);
	err = CheckBoardTypeId(defs,boardTypeId);
	if err ~= 0, return, end
	
	dev_Alazar = dev;
	OPEN_Alazar = 1;
end
dev = dev_Alazar;

    

switch lower(op)
    
	case 'close'
		unloadlibrary(lib);
		OPEN_Alazar = 0;
		
	case 'initialize_and_run'
		numvar			= numel(varargin);
		if numvar > 0
	        filename		= varargin{1};
		else
			filename = '';
		end
		if numvar > 1
	        P				= varargin{2};
		else
			P = [];
		end
		if numvar > 2
% 	        args			= varargin{6}
			args			= varargin{3:end};
% not used?
        end
        
        
        
		[err,S,ch] = ChannelDefaults(err,defs,S,ch,show);
             		
		[err,T] = TriggerDefaults(err,defs,T,ch); % subfunc
       
        % check input parameters and determine some params
		[err,S] = CheckInputs(lib,dev,err,defs,S,T,show,analyzeData);
        if ~isfield(S,'BytesPerBuffer')
            data = [];
            varargout{1} = err;
            return; 
        end;
		                
		% set record size
		err = AlazarSetRecordSize(lib,dev,err,S,T,show);	% subfunc
	        
		% setup clock
		err = AlazarSetCaptureClock(lib,dev,err,defs,S,T,show);	% subfunc

		% setup channels
		err = AlazarInputControlEX(lib,dev,err,defs,S,ch,show);	% subfunc
        
		% setup trigger
        disp('setting trigger')
		err = setupTrigger(lib,dev,err,T,show); % subfunc
		disp('done setting trigger')
		if show, disp(['creating buffers, BytesPerBuffer: ',num2str(S.BytesPerBuffer)]), end
		for m = 1:S.numBuffers
			if ~err
% 				pBufferData(m) = libpointer('uint8Ptr',zeros(1,S.BytesPerBuffer+16));
				pBufferData(m) = libpointer('uint8Ptr',zeros(1,S.BytesPerBuffer));
% 				pBufferData(m) = libpointer('uint16Ptr',zeros(1,S.BytesPerBuffer/S.BytesPerSample));
			else
				break
			end
		end
		
        if dosave
			if isempty(filename)	% generate default filename
				filename = getvalfile('daq');
			elseif strcmp(filename(end-3),'.')	% stip .ext if exists
				filename = filename(1:end-4);
			end
			
			% open file
            if ~append		% overwrite
%                 if show, disp('opening new data file'), end
                fid = fopen([filename,'.daq'], 'w');
			else			% append
%                 if show, disp('opening existing data file'), end
                fid = fopen([filename,'.daq'], 'a');
            end
		end
        
		if popup, FS = stoploop('stop acquisition'); end
        if show, disp('starting acquisition'), end
        
        first_call = 1;
    
        for reps = 1:S.Repeats
            
            
			if ~err
	            err = AlazarBeforeAsyncRead(lib,dev,err,defs,S,T,show);
			else
				break
			end

            % send buffer addresses to alazar
            if show, disp(['sending buffer addresses to alazar']), end
            for m = 1:S.numBuffers
				if ~err
					err = AlazarPostAsyncBuffer(lib,dev,err,pBufferData(m),S);
				else
					break
				end
            end		

            if show, disp('acquiring data ...'), end
            
            tic;
            
% could move this into Buffer loop to stop more quickly...
            if popup && FS.Stop() && err == 0
				err = AlazarAbortAsyncRead(lib,dev,err);
                disp('acquistion stopped by user')
				break
            end
            
            % starting async capture
            err = AlazarStartCapture(lib,dev,err,show);	% subfunc
            
            waittime         = 50000;	% ms 
% 			% calculate waittime?
% 			freqlow = S.SampleRate / S.RecordLength;
% 			t = S.RecordCountPerChannel / freqlow / S.numBuffers;
% 			waittime = t * 5;	% wait 5 times as long as needed
			
            BuffersCompleted = 0;
                                                            
            for BufferIndex = 1:S.numBuffers
                
                % stream to disk should be here. However, did not get it to
                % run - do we have to buy it?
                %[sts,ErrString,A] = AlazarCreateStreamFileA(dev,'data.bin');
                %AlazarErrorCheck(lib,errorCode);
                
				% wait for a buffer to complete
				err = AlazarWaitAsyncBufferComplete(lib,dev,err,pBufferData(BufferIndex),waittime);
				
				% check whether any samples overrange
				err = AlazarGetStatus(lib,dev,err);
				
				if ~err
                    BuffersCompleted = BuffersCompleted + 1;
                    if show, disp(['Buffers completed: ',num2str(BuffersCompleted),' of ',num2str(S.numBuffers)]), end
                    if analyzeData
						[err,P] = Analyze_realtime(err,first_call, pBufferData(BufferIndex), P, S, T, ch, varargin); 
                    end
                    if dosave
						% if user cancels, then correct number of repeats
						% want to write current Repeats to file each time
						actRepeats = S.Repeats;
						S.Repeats = reps;	
						err = SaveRawData(fid,filename,S,ch,T,P,pBufferData(BufferIndex));
						% not passing err through...
						S.Repeats = actRepeats;	% reset to correct number
% 					elseif return_traces
					end
					% tell Alazar this buffer is now free
					err = AlazarPostAsyncBuffer(lib,dev,err,pBufferData(BufferIndex),S);
				else
                    disp('ERROR: problem with buffer acquisition')
					break
                end
                first_call = 0;
         end
            
            tm = toc;
			if show
	            disp(sprintf('%3d/%d; average rate: %.3f MS/sec',reps,S.Repeats,S.RecordLength*S.RecordCountPerChannel*S.numchannels2record/tm/1e6))
			end
            
			err = AlazarAbortAsyncRead(lib,dev,err);
   
		end
        
		if popup
	        FS.Clear() ; % Clear up the box
		    clear FS;
		end
        
%         err = AlazarAbortAsyncRead(lib,dev,err);
		   
		if dosave
			fclose(fid);
		end
		
		if err
			disp('unloading library')
	 		unloadlibrary(lib);
			OPEN_Alazar = 0;
		end
		
		if graph && dosave && ~err
			load_Alazar(filename,'graph');
		end
		
		% varargout
		varargout{1} = err;
		varargout{2} = S;
		varargout{3} = ch;
		varargout{4} = T;
		if dosave
			varargout{5} = filename;
		else
			varargout{5} = '';
		end
		varargout{6} = P;
		if return_traces
			varargout{7} = load_Alazar(filename);
		else
			varargout{7} = [];
		end


	case {'voltmeter','dmm'}
		% inputs [op,S,ch,T,varargin]
		% for this 
		ch2use = ch;
		ch = [];
		numarg = numel(varargin);
		if ~isempty(varargin{1})
			InputRange = varargin{1};
			ch.(ch2use).InputRange = InputRange;
		end
		if ~isempty(varargin{2})
			S.RecordLength  = varargin{2};
		end
		
		ch.(ch2use).Coupling = 'dc';
		ch.(ch2use).InputImped = 1e6;

		% just want trigger to timeout (free run)
		T.TrigEngSource = ch2use;
		T.TriggerTimeOut = 1;

		% combine all varargins
		filename = 'temp_volt_file';
		varargin2 = {filename,[],'save','return_traces','nopopup'};
		for m = 1:numel(varargin)
			varargin2 = [varargin2,varargin{m}];
		end
		[err,S,ch,T,filename,P,data3D] = Alazar('initialize_and_run',...
            S,ch,T,varargin2);
		
		% delete temp file
		delete([filename,'.daq'])
		delete([filename,'.mat'])
				
		% convert data to 1D and find mean
		[s1,s2,s3] = size(data3D);
		data1D = reshape(data3D,1,s1*s2*s3);
		varargout{1} = mean(data1D);
		
    		if graph
				figure; 
				plot(data1D)
			end

	case {'autovoltmeter','autodmm'}
		ch2use = ch;
		ch = [];
		numarg = numel(varargin);
		if ~isempty(varargin{1})
			S.RecordLength  = varargin{1};
		end
		
		ch.(ch2use).Coupling = 'dc';
		ch.(ch2use).InputImped = 1e6;
		
		% start at max range, then decrease after get a guess
		ch.(ch2use).InputRange = max(defs.inputranges_1000000);

% 		rangegood = 0
% 		while rangegood == 0
		
		% just want trigger to timeout (free run)
		T.TrigEngSource = ch2use;
		T.TriggerTimeOut = 1;

		% combine all varargins
		filename = 'temp_volt_file';
		varargin2 = {filename,[],'save','return_traces','nopopup'};
		for m = 1:numel(varargin)
			varargin2 = [varargin2,varargin{m}];
		end
		[err,S,ch,T,filename,P,data3D] = Alazar('initialize_and_run',...
            S,ch,T,varargin2);
		
		% convert data to 1D and find min & max
		[s1,s2,s3] = size(data3D);
		data1D = reshape(data3D,1,s1*s2*s3);
		vmin = abs(min(data1D));
		vmax = abs(max(data1D));
		minrange = max([vmin,vmax]);

		if graph
			figure; hold all
			plot(data1D)
		end
		
		% decrease input range close to measured value
		iranges = find(defs.inputranges_1000000 > minrange);
		if ~isempty(iranges)
			irange = iranges(1);
% 		else
% 			irange = [];
		end
% 		if ~isemtpy(irange) && irange ~= numel(defs.inputranges_1000000)	
		if irange ~= numel(defs.inputranges_1000000)	
			% less than max range, so set lower
			
			ch.(ch2use).InputRange = defs.inputranges_1000000(iranges(1));

			% take data again with new input range
			[err,S,ch,T,filename,P,data3D] = Alazar('initialize_and_run',...
				S,ch,T,varargin2);

			% convert data to 1D and find mean
			[s1,s2,s3] = size(data3D);
			data1D = reshape(data3D,1,s1*s2*s3);
			varargout{1} = mean(data1D);
		else	% max range - already took at max range
			irange = numel(defs.inputranges_1000000);
			varargout{1} = mean(data1D);
		end
		varargout{2} = defs.inputranges_1000000(irange);

		% delete temp file
		delete([filename,'.daq'])
		delete([filename,'.mat'])
				
		if graph
% 			figure; 
			plot(data1D)
		end
	otherwise            
		error(['unknown operation (op)']);
end    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%        subfunctions      %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ---------------------
function [err,S] = CheckInputs(lib,dev,err,defs,S,T,show,analyzeData)
% 
% 
%! thomas 
%	need to modify checks for more than 1 channel
%		for some reason, can only handle 2^18 with 2 channels
%		but it works for 2^29 for just 1!!

if ~err

	
	if ~isfield(S,'SampleRate')	% default to max samplerate
		S.SampleRate = max(defs.SampleRates);
	end
	if ~isfield(S,'RecordLength')	% default recordlength
		S.RecordLength = 2^8;
	end
	if ~isfield(S,'RecordCountPerChannel') % default recordcountperchannel
		S.RecordCountPerChannel = 2^4;
	end
	if ~isfield(S,'Repeats')	% default - don't repeat
		S.Repeats = 1;
	end
	if ~isfield(S,'numBuffers')	% default number of buffers
		S.numBuffers = 4;
	end
	
	% enable or disable 20MHz BW filter
	if isempty(defs.BWLimit)	% card does not support BW limit
		S.BWLimit = -2;
	elseif ~isfield(S,'BWLimit')	% default BWLimit auto (-1=auto, 0=off, 1=on)
		S.BWLimit = -1;
	end

%  		if show
% 			disp([	'RecordLength ',num2str(S.RecordLength),...
% 					', predepth ',num2str(T.PreDepth),...
% 					', buffers ',num2str(S.numBuffers),...
% 					', BytesPerBuf ',num2str(S.BytesPerBuffer),...
% 					', reccountperch ',num2str(S.RecordCountPerChannel)]);
% 		end
 
   
	
	% get Alazar channel info
	[retCode,pdummy,MemSizePerChannel,S.BitsPerSample] = calllib(lib,'AlazarGetChannelInfo',dev,1,1);
	S.BitsPerSample = double(S.BitsPerSample);
	err = AlazarErrorCheck(lib,retCode,'AlazarGetChannelInfo');
% 	S.BytesPerSample = double(ceil(S.BitsPerSample/8));
	S.BytesPerSample = floor((double(S.BitsPerSample) + 7) / double(8));
		% from Alazar sample code
	if show
		disp([	'MemSizePerChannel: ',num2str(MemSizePerChannel),...
				', BitsPerSample: ',num2str(S.BitsPerSample)])
	end

	
	
	% need to determine RecsPerChannelPerBuffer
	% also check whether realtime analysis OK if selected
	% check numBuffers is enough
	maxRecsPerChannelPerBuffer = defs.maxbuffersize/S.RecordLength/S.BytesPerSample/S.numchannels2record;
	if maxRecsPerChannelPerBuffer < 1
		disp(['ERROR: need to decrease RecordLength by factor of ',num2str(1/maxRecsPerChannelPerBuffer)])
		err = 1;
		return
	end
	if analyzeData
		maxRecsPerChannelPerBuffer = maxRecsPerChannelPerBuffer / 2;
		% need twice as much memory to copy data to other variables
	end
	fullRecsPerChannelPerBuffer = S.RecordCountPerChannel / S.numBuffers;
		% see whether can fit all RecordCountPerChannel in one go (realtime
		% check)
	if fullRecsPerChannelPerBuffer <= maxRecsPerChannelPerBuffer
		% realtime analyze ok with input #Buffers
		S.RecsPerChannelPerBuffer = fullRecsPerChannelPerBuffer;
	else
		S.RecsPerChannelPerBuffer = maxRecsPerChannelPerBuffer;
		if analyzeData
			disp(['ERROR: too much data for realtime analysis. Reduce data by factor of ',...
				num2str(fullRecsPerChannelPerBuffer/maxRecsPerChannelPerBuffer)])
			err = 1;
			return
		else
			S.numBuffers = S.RecordCountPerChannel / maxRecsPerChannelPerBuffer;
			disp(['warning: had to increase numBuffers to ',num2str(S.numBuffers)])
		end
	end
	if show, disp(['RecsPerChannelPerBuffer ',num2str(S.RecsPerChannelPerBuffer)]), end

	S.BytesPerBuffer = S.RecsPerChannelPerBuffer * S.RecordLength * S.numchannels2record * S.BytesPerSample;
	if show, disp(['BytesPerBuffer ',num2str(S.BytesPerBuffer)]), end


	if S.BytesPerBuffer > defs.maxbuffersize
        disp(['ERROR: max buffer size overrun! maxbuffersize:',num2str(defs.maxbuffersize),', BytesPerBuffer:',num2str(S.BytesPerBuffer)])
        err = 1;
        return
    end
	if S.BytesPerBuffer > defs.maxbuffersize
        disp(['ERROR: max buffer size overrun! maxbuffersize:',num2str(defs.maxbuffersize),', BytesPerBuffer:',num2str(S.BytesPerBuffer)])
        err = 1;
        return
    end
	
	
% 	log2(S.BytesPerBuffer*S.numBuffers*2)
% 	if S.BytesPerBuffer*S.numBuffers*2 > sV.PhysicalMemory.Available
%         disp(sprintf('System memory overload! Max Memory Available: %d bytes. Bytes per buffer: %d',sV.PhysicalMemory.Available,S.BytesPerBuffer))
%         err = 1;
%         return
%     end


	% get matlab & computer memory info
	[uV sV] = memory;

	% only realtime limit
% Thomas - need a x2 in here? 
%	also switched from mem available to maxarraybytes
%     if analyzeData && S.RecordLength*S.RecordCountPerChannel*S.numchannels*S.BytesPerSample*2 > sV.PhysicalMemory.Available
% 	 if S.BytesPerBuffer*S.numBuffers*2 > sV.PhysicalMemory.Available
	 if S.BytesPerBuffer*S.numBuffers > uV.MaxPossibleArrayBytes
		disp(['ERROR: System memory overload. ',...
				'Max Memory Available: ',...
				num2str(uV.MaxPossibleArrayBytes),' bytes. ',...
				'Reduce pts by factor of ',...
				num2str(S.BytesPerBuffer*S.numBuffers/uV.MaxPossibleArrayBytes)]);
        err = 1;
        return
   end
        
 	% only realtime limit
   if analyzeData && S.RecordLength*S.RecordCountPerChannel > double(MemSizePerChannel)
        disp('ERROR: board memory overload!')
        err = 1;
        return
    end

    if S.numBuffers < 2 % don't know why yet, but seems to be the case
        disp('ERROR: need at least 2 buffers per channel!')
        err = 1;
        return
	end

    if T.PreDepth >= S.RecordLength
        disp('ERROR: predepth (pretrigger) must be less than RecordLength')
        err = 1;
        return
	end

	if defs.maxPreTrigger < 0	% code for maxPreTrigger = RecLen - #
		lesspts = abs(defs.maxPreTrigger);
		defs.maxPreTrigger = S.RecordLength - lesspts;
			% set maxPreTrigger based on RecordLength
		if T.PreDepth > S.RecordLength + defs.maxPreTrigger
			disp(['ERROR: predepth (pretrigger) must be less than RecordLength - ',num2str(lesspts),' for this card'])
			err = 1;
		end
	else
	    if T.PreDepth > defs.maxPreTrigger
			disp(['ERROR: predepth (pretrigger) must be less than ',num2str(defs.maxPreTrigger),' for this card'])
			err = 1;
			return
		end
	end

%!!thomas!! I don't think this is correct - I thought that was just min
%recordlength...
%     if floor(S.RecordLength/double(BitsPerSample))*double(BitsPerSample) ~= S.RecordLength
%         disp(sprintf('record length must be a multiple of %d!',BitsPerSample))
%         err = 1;
%         return
%     end

%     if S.RecordCountPerChannel/S.numBuffers < 2 % seems to be the case
%         disp('ERROR: #Records per channel must be at 2x #Buffers!')
%         err = 1;
%         return
% 	end
   

    % check minsample for all cards (in manual)
    if S.RecordLength * S.RecordCountPerChannel < defs.minsample
        disp(sprintf('ERROR: sum of samples must be bigger than %d for each acquistion! is: %d',defs.minsample,S.RecordLength * S.RecordCountPerChannel))
        err = 1;
        return
	end
	
	% check whether FIFO mode input. if not, set to default
	if defs.FIFO == 0 % FIFO not supported for current card
		S.FIFO = 0;
	elseif ~isfield(S,'FIFO') % default FIFO on if supported
		S.FIFO = 1;
	end
end
end

% ---------------------
function [err,S,ch] = ChannelDefaults(err,defs,S,ch,show)
% parse through  struct 'ch' and get values for each channel
% setup all channels even if not recorded - avoids overrange errors
% check for whether fields were input - if not, use defaults
% 
% ex: ch.A.InputRange = 1; ch.C.InputRange = 2;
%		setup channels A & C
% (using dynamic field names for chname)
% 
% if want to setup channel but NOT record, set .record = 0
%	ex: ch.B.record = 0 to setup B for trigger, but not save traces
% 
% 
% !! not sure how to deal with channelsum (SDK calls it uChannelMode)
% !! the SDK manual I have only has A&B options - don't know what to do
%	about C&D

if ~err
    chfieldnames = fieldnames(ch);
	inputchannelnames = '';
	for m = 1:numel(chfieldnames)
		inputchannelnames = [inputchannelnames,chfieldnames{m}];
	end
		% ex: makes inputchannelnames = 'AB'
% 	channelnames = {'A','B','C','D'};
	channelnames = 'ABCD';
	S.numchannels2record = 0;
	S.channels2record = [];
	S.channelsum = 0;
	for m = 1:defs.numchannels
		chname = channelnames(m);

		% figure out whether to record this channel or not
		if ~isempty(strfind(inputchannelnames,chname))
			% current ch (chname) was in input ch. struct
			% need to find out whether ch should be recorded or not
			if ~isfield(ch.(chname),'record')
				ch.(chname).record = 1;	% setup channel and also record traces
			end
		else	% current ch (chname) not input in ch. struct - don't record
			ch.(chname).record = 0;
		end
		
		% default to 1 MOhm impedance
        if ~isfield(ch.(chname),'InputImped')
  			ch.(chname).InputImped = 1e6;	% default 1MOhm
		end
		[err,ch.(chname).InputImpedConv] = ConvertInputImped(err,defs,ch.(chname).InputImped);	% convert input impedance

		% default to largest input range
		if ~isfield(ch.(chname),'InputRange')
			switch ch.(chname).InputImped
				case 50
					ch.(chname).InputRange = max(defs.inputranges_50);	
				case 1e6
					ch.(chname).InputRange = max(defs.inputranges_1000000);	
			end
		end	
		[err,ch.(chname).InputRangeConv] = ConvertInputRange(err,defs,ch.(chname).InputRange,ch.(chname).InputImped);	% convert input range
		
		if ~isfield(ch.(chname),'Coupling')
			ch.(chname).Coupling = 'dc';	% default DC coupling
		end
		[err,ch.(chname).CouplingConv] = ConvertCoupling(err,defs,ch.(chname).Coupling);	% convert coupling
		
		switch chname
			case 'A'
				ch.(chname).chvalue = defs.CHANNEL_A;
			case 'B'
				ch.(chname).chvalue = defs.CHANNEL_B;
			case 'C'
				ch.(chname).chvalue = defs.CHANNEL_C;
			case 'D'
				ch.(chname).chvalue = defs.CHANNEL_D;
		end

		
		% channels2record examples
		%	ch A&B, [1,2]
		%	ch A&C, [1,3]
		if ch.(chname).record == 1
			S.numchannels2record = S.numchannels2record + 1;
			switch chname
				case 'A'
					S.channels2record = [S.channels2record,1];
				case 'B'
					S.channels2record = [S.channels2record,2];
				case 'C'
					S.channels2record = [S.channels2record,3];
				case 'D'
					S.channels2record = [S.channels2record,4];
			end
			
			% also find channel sum value 
			% (code to alazar which channels to record)
			S.channelsum = S.channelsum + ch.(chname).chvalue;
		end
		
	end
else
% 	ch = [];
% 	numchannels2record = [];
end
end

% %-------------------------
% function channelsum = findchannelsum(S,ch)
% 
% % channel values: A = 1; B = 2; C = 4; D = 8; see def file
% % for channels that need to be setup (for triggering) but not recorded,
% %	ch.x.record = 0 so do not add to channelsum
% 
% 	channelnames = fieldnames(ch);
% % 	numchannels2setup = numel(channelnames);
% 	channelsum = 0;
% 	for m = 1:S.numchannels2record
% 		chname = channelnames{S.channels2record(m)};
% 	for m = 1:numchannels2setup
% 		chname = channelnames{m};
% 		if ch.(chname).record == 1;	% only add to sum if want to record traces
% 			channelsum = channelsum + ch.(chname).chvalue;
% 		end
% 	end
% end


%-------------------------
function [err,T] = TriggerDefaults(err,defs,T,ch)
% parse through  struct 'T' and setup trigger

	if isempty(T)
		T=struct; 	% just use all defaults
	end

	if ~isfield(T,'PreDepth')
		T.PreDepth = 0;
    end
    
    if ~isfield(T,'Clock')
		T.Clock = 'internal'; 
	end

	
	if ~isfield(T,'ClockEdge')
		T.ClockEdge = defs.CLOCK_EDGE_RISING; 
	end

	if ~isfield(T,'TrigEngSource')
		T.TrigEngSource = 'external';	% default external 
	end
	
	% trigger source
	switch lower(T.TrigEngSource)
		case {'external','ext'}
			if ~isfield(T,'InputRange')
				T.InputRange = 5; 	% default 5V (other choice 1)
			end
			[err,T.InputRangeConv] = ConvertTriggerInputRange(err,defs,T.InputRange);
			
			if ~isfield(T,'Coupling')
				T.Coupling = 'dc'; 	% default DC coupling
			end
			[err,T.CouplingConv] = ConvertCoupling(err,defs,T.Coupling);	% convert coupling
			
			% for external, both engines equally bad, so just pick K
			if ~isfield(T,'TrigEngOperation')
				T.TrigEngOperation = 'K'; 
			end
		case 'a'
			T.TrigEngOperation = 'J';	% engine J works best for ChA triggering
			if ~isfield(ch.A,'InputRange')
				disp('ERROR: ch.A.InputRange not specified - needed to set trigger level')
				err = 1;
			else
				T.InputRange = ch.A.InputRange;
			end
		case 'b'
			T.TrigEngOperation = 'K';	% engine K works best for ChB triggering
			if ~isfield(ch.B,'InputRange')
				disp('ERROR: ch.B.InputRange not specified - needed to set trigger level')
				err = 1;
			else
				T.InputRange = ch.B.InputRange;
			end
	end
	[err,T.TrigEngOperationConv] = ConvertTrigEngOperation(err,defs,T.TrigEngOperation);
	
	% 1 is always j and 2 is always k
	T.TrigEng1 = defs.TRIG_ENGINE_J;
	T.TrigEng2 = defs.TRIG_ENGINE_K;

	if ~isfield(T,'TrigEngSlope')
		T.TrigEngSlope = 1;		% positive
	end
	
	if ~isfield(T,'TrigEngLevel')
		T.TrigEngLevel = 0.5; 	% default 0.5 Volts
	end

	% always have to setup both trig1 & trig2 
	switch lower(T.TrigEngOperation)
		case 'j'	% setup J, and use defaults for K (not used)
			% trigger source
			[err,T.TrigEngSource1Conv] = ConvertTrigEngSource(err,defs,T.TrigEngSource);	%convert
			[err,T.TrigEngSource2Conv] = ConvertTrigEngSource(err,defs,'disable');

			% trigger slope
			[err,T.TrigEngSlope1Conv] = ConvertTrigEngSlope(err,defs,T.TrigEngSlope);	%convert
			[err,T.TrigEngSlope2Conv] = ConvertTrigEngSlope(err,defs,1);	% dummy

			% trigger level
			[err,T.TrigEngLevel1Conv] = ConvertTrigEngLevel(err,T.TrigEngLevel,T.InputRange);
			T.TrigEngLevel2Conv = 128;	% dummy value
		case 'k'
			% trigger source
			[err,T.TrigEngSource1Conv] = ConvertTrigEngSource(err,defs,'disable');
			[err,T.TrigEngSource2Conv] = ConvertTrigEngSource(err,defs,T.TrigEngSource);	%convert

			% trigger slope
			[err,T.TrigEngSlope1Conv] = ConvertTrigEngSlope(err,defs,1);	% dummy
			[err,T.TrigEngSlope2Conv] = ConvertTrigEngSlope(err,defs,T.TrigEngSlope);	%convert

			% trigger level
			T.TrigEngLevel1Conv = 128;	% dummy value
			[err,T.TrigEngLevel2Conv] = ConvertTrigEngLevel(err,T.TrigEngLevel,T.InputRange);
	end
	
	if ~isfield(T,'TriggerTimeOut')
		T.TriggerTimeOut = 100; 	% 100kHz clock
	end

	if ~isfield(T,'TriggerDelay')
		T.TriggerDelay = 0; 	% units of samples
	end
end

%-------------------------
function err = setupTrigger(lib,dev,err,T,show)
% parse through  struct 'T' and setup trigger

	% setup trigger
	err = AlazarSetTriggerOperation(lib,dev,err,T,show);	%subfunc
    
	% setup trigger delay
	err = AlazarSetTriggerDelay(lib,dev,err,T,show);	% subfunc.

	% setup external trigger if set to external
	switch lower(T.TrigEngSource)
		case {'ext','external'}
			err = AlazarSetExternalTrigger(lib,dev,err,T,show);	% subfunc
	end

	% setup trigger timeout, multiple of 100 kHz
	err = AlazarSetTriggerTimeOut(lib,dev,err,T,show);	% subfunc
end

%-------------------------
function err = AlazarSetRecordSize(lib,dev,err,S,T,show)
	% set record size
if ~err	
	if show
		disp(['set record size to ',num2str(S.RecordLength)])
	end
	retCode = calllib(lib,'AlazarSetRecordSize',dev,T.PreDepth,S.RecordLength-T.PreDepth);
	err = AlazarErrorCheck(lib,retCode,'AlazarSetRecordSize');
end
end

%-------------------------
function err = AlazarSetCaptureClock(lib,dev,err,defs,S,T,show)
	% setup clock
    
if isfield(T,'Clock')
      switch T.Clock
            case 'fast_external'
                T.ClockSource = defs.FAST_EXTERNAL_CLOCK;
                if show, disp('clock set to fast external clock'), end
            case 'medium_external'
                T.ClockSource = defs.MEDIUM_EXTERNAL_CLOCK;
                if show, disp('clock set to medium external clock'), end
            case 'slow_external'
                T.ClockSource = defs.SLOW_EXTERNAL_CLOCK;
                if show, disp('clock set to slow external clock'), end
            case 'internal'
                T.ClockSource = defs.INTERNAL_CLOCK;
                if show, disp('clock set to internal clock'), end
    end
else
    T.ClockSource = defs.INTERNAL_CLOCK;
    if show, disp('clock set to internal clock'), end
end
  
if ~err	
	[err,SampleRateConv] = ConvertSampleRate(err,defs,S.SampleRate);	% convert sample rate
	%if show, disp(['setup clock']), end
	retCode = calllib(lib,'AlazarSetCaptureClock',dev,...
					   T.ClockSource,...
					   SampleRateConv,...	% converted
					   T.ClockEdge,...
					   0);	% 0 means disable decimation
	err = AlazarErrorCheck(lib,retCode,'AlazarSetCaptureClock');
end
end

%-------------------------
function err = AlazarInputControlEX(lib,dev,err,defs,S,ch,show)
if ~err
	channelnames = fieldnames(ch);
	% setup each channel one by one
	for m = 1:defs.numchannels
		chname = channelnames{m};

		if ~err	
			if show, disp(['setup channel ',chname]), end
			% setup ch - use converted values
			retCode = calllib(lib,'AlazarInputControl',dev,...
								ch.(chname).chvalue,...
								ch.(chname).CouplingConv,...
								ch.(chname).InputRangeConv,...
								ch.(chname).InputImpedConv);
			err = AlazarErrorCheck(lib,retCode,'AlazarInputControl');
		else
			return
		end
		
		% enable or disable 20MHz BW filter
		if S.BWLimit == -1	% auto set = BWlimit on for samplerate <= 40 MHz
			if S.SampleRate > 2*defs.BWLimit	% do not want to limit BW
				S.BWLimit = 0;
			else	% enable bandwidth limit
				S.BWLimit = 1;
			end				
		end
		switch S.BWLimit
			case 0
				if S.SampleRate < max(defs.SampleRates)
					disp(['warning: when samplerate below ',num2str(max(defs.SampleRates)*1e-6),' MHz can get aliasing'])
				end
			case 1
				if S.SampleRate < 2*defs.BWLimit
					disp(['warning: when samplerate below ',num2str(2*defs.BWLimit*1e-6),' MHz can get aliasing'])
				end
        end
        err = AlazarSetBWLimit(lib,dev,err,ch.(chname).chvalue,S.BWLimit,show);
	end
end
end

%-------------------------
function err = AlazarSetTriggerOperation(lib,dev,err,T,show)
	% setup trigger
if ~err
	if show
		disp(['setup trigger level'])
% 		disp(['setup trigger level ',num2str(T.TrigEngLevel1Conv)])
	end
	retCode = calllib(lib,'AlazarSetTriggerOperation',dev,...
						T.TrigEngOperationConv,...
						T.TrigEng1,...
						T.TrigEngSource1Conv,...
						T.TrigEngSlope1Conv,...
						T.TrigEngLevel1Conv,...
						T.TrigEng2,...
						T.TrigEngSource2Conv,...
						T.TrigEngSlope2Conv,...
						T.TrigEngLevel2Conv) ;              
	err = AlazarErrorCheck(lib,retCode,'AlazarSetTriggerOperation');	
end
end

%-------------------------
function err = AlazarSetTriggerDelay(lib,dev,err,T,show)
if ~err
	if show
		disp(['setup trigger delay = ',num2str(T.TriggerDelay),' clock cycles'])
	end
	retCode = calllib(lib,'AlazarSetTriggerDelay',dev,T.TriggerDelay);              
	err = AlazarErrorCheck(lib,retCode,'AlazarSetTriggerDelay');
end
end

%-------------------------
function err = AlazarSetExternalTrigger(lib,dev,err,T,show)
if ~err
	% setup external trigger
	if show, disp(['set external trigger']), end
	retCode = calllib(lib,'AlazarSetExternalTrigger',dev,...
						T.CouplingConv,...
						T.InputRangeConv) ;              
	err = AlazarErrorCheck(lib,retCode,'AlazarSetExternalTrigger');
end
end

%-------------------------
function err = AlazarSetTriggerTimeOut(lib,dev,err,T,show)
	% setup trigger timeout, multiple of 100 kHz
if ~err
	if show
		disp(['set timeout ',num2str(T.TriggerTimeOut/100000),' sec'])
	end
	retCode = calllib(lib,'AlazarSetTriggerTimeOut',dev,T.TriggerTimeOut);
	err = AlazarErrorCheck(lib,retCode,'AlazarSetTriggerTimeOut');
end
end

%-------------------------
function err = AlazarBeforeAsyncRead(lib,dev,err,defs,S,T,show)
% setting up before async read

if ~err
    ADMA_EXTERNAL_STARTCAPTURE = 1;
	ADMA_TRADITIONAL_MODE = 0;
    ADMA_FIFO_ONLY_STREAMING = hex2dec('00000800');
    
    if S.FIFO == 1 && defs.FIFO == 1
        flags = ADMA_EXTERNAL_STARTCAPTURE + ADMA_TRADITIONAL_MODE + ADMA_FIFO_ONLY_STREAMING;
		if show
			disp(['setting up before async read with FIFO'])
		end
    else
        flags = ADMA_EXTERNAL_STARTCAPTURE + ADMA_TRADITIONAL_MODE;
		if show
			disp(['setting up before async read without FIFO'])
		end
	end
    	
	retCode = calllib(lib,'AlazarBeforeAsyncRead',dev,...
						S.channelsum,...                        
						-T.PreDepth,...
						S.RecordLength,...
						S.RecsPerChannelPerBuffer*S.numchannels2record,...
						S.RecordCountPerChannel,...
						flags);
	err = AlazarErrorCheck(lib,retCode,'AlazarBeforeAsyncRead');
end
end

%-------------------------
function err = AlazarStartCapture(lib,dev,err,show)
if ~err
	if show
		disp(['starting async capture'])
	end
	retCode = calllib(lib,'AlazarStartCapture',dev);
	err = AlazarErrorCheck(lib,retCode,'AlazarStartCapture');
end
end

%-------------------------
function err = AlazarWaitAsyncBufferComplete(lib,dev,err,pData,waittime)
if ~err
	retCode =  calllib(lib,'AlazarWaitAsyncBufferComplete',dev,pData,waittime);
	err = AlazarErrorCheck(lib,retCode,'AlazarWaitAsyncBufferComplete');
end
end

%-------------------------
function err = AlazarGetStatus(lib,dev,err)
if ~err
	[status,retCode] =  calllib(lib,'AlazarGetStatus',dev);
	% ! retCode doesn't seem to work, so disable errorcheck
% 	err = AlazarErrorCheck(lib,retCode,'AlazarGetStatus');
	switch status
		case 0	% everything OK
		case 1	% trigger timeout - maybe desired, so don't do anything
		case 2
			disp('ERROR: At least one chA sample was out of range during the last acquisition.')
			err = 1;
% 		case 3
% 			disp('ERROR: not sure, but think: Trigger timeout and at least one chA sample was out of range during the last acquisition.')
% 			err = 1;
		case 4
			disp('ERROR: At least one chB sample was out of range during the last acquisition.')
			err = 1;
% 		case 5
% 			disp('ERROR: not sure, but think: Trigger timeout and at least one chB sample was out of range during the last acquisition.')
% 			err = 1;
% 		case 6
% 			disp('ERROR: not sure, but think: At least one chA & one chB sample was out of range during the last acquisition.')
% 			err = 1;
% 		case 7
% 			disp('ERROR: not sure, but think: Trigger timeout & at least one chA & one chB sample was out of range during the last acquisition.')
% 			err = 1;
		case 8
			disp('ERROR: PLL locked (660 only).')
			err = 1;
		otherwise
			disp(['WARNING: unknown status code from AlazarGetStatus ',num2str(status)])
% 			err = 1;
	end
end
end

%-------------------------
function err = AlazarPostAsyncBuffer(lib,dev,err,pData,S)
if ~err
	retCode = calllib(lib,'AlazarPostAsyncBuffer',dev,pData,S.BytesPerBuffer);
	err = AlazarErrorCheck(lib,retCode,'AlazarPostAsyncBuffer');
end
end

%-------------------------
function err = AlazarAbortAsyncRead(lib,dev,err)
if ~err
	retCode = calllib(lib,'AlazarAbortAsyncRead',dev);
	err = AlazarErrorCheck(lib,retCode,'AlazarAbortAsyncRead');
end
end

%-------------------------
function err = AlazarSetBWLimit(lib,dev,err,chvalue,enableBWfilter,show)
% Set the bandwidth for a given channel to either its maximum value
% or the reduced value. Only the ATS460, ATS660, ATS860 and ATS9462
% support this functionality – reduced bandwidth is approximately 20 MHz
% ATS 9440 does not support this feature.
if ~err
	switch enableBWfilter
		case 0
			if show, disp(['disabling BW filter']), end
		case 1
			if show, disp(['enabling BW filter']), end
        case -2
            if show, disp(['ATS9440 does not support BW limit']), return, end
    end
    retCode = calllib(lib,'AlazarSetBWLimit',dev,chvalue,enableBWfilter);
    err = AlazarErrorCheck(lib,retCode,'AlazarSetBWLimit');
end
end

%-------------------------
function [err,RateIndex] = ConvertSampleRate(err,defs,SampleRateNum)
% 	AlazarDefs
try
	% check whether SampleRate valid for card
	if isempty(find(defs.SampleRates == SampleRateNum))
		% input SampleRate not in list of acceptable rates
		err = 1;
		disp(['ERROR: invalid sample rate for selected card.'])
		disp(['accepted samplerates are:',num2str(defs.SampleRates)])
	end
	switch SampleRateNum
        case 180e6
			RateIndex = defs.SAMPLE_RATE_180MSPS;
        case 180e6
			RateIndex = defs.SAMPLE_RATE_180MSPS;
		case 125e6
			RateIndex = defs.SAMPLE_RATE_125MSPS;
		case 100e6
			RateIndex = defs.SAMPLE_RATE_100MSPS;
		case 50e6
			RateIndex = defs.SAMPLE_RATE_50MSPS;
		case 25e6
			RateIndex = defs.SAMPLE_RATE_25MSPS;
		case 20e6
			RateIndex = defs.SAMPLE_RATE_20MSPS;    
		case 10e6
			RateIndex = defs.SAMPLE_RATE_10MSPS;    
		case 5e6
			RateIndex = defs.SAMPLE_RATE_5MSPS;  
		case 2e6
			RateIndex = defs.SAMPLE_RATE_2MSPS;   
		case 1e6
			RateIndex = defs.SAMPLE_RATE_1MSPS;  
		case 500e3
			RateIndex = defs.SAMPLE_RATE_500KSPS;  
		case 200e3
			RateIndex = defs.SAMPLE_RATE_200KSPS;
		case 100e3
			RateIndex = defs.SAMPLE_RATE_100KSPS;
		case 50e3
			RateIndex = defs.SAMPLE_RATE_50KSPS; 
		case 20e3
			RateIndex = defs.SAMPLE_RATE_20KSPS; 
		case 10e3
			RateIndex = defs.SAMPLE_RATE_10KSPS;      
		case 5e3
			RateIndex = defs.SAMPLE_RATE_5KSPS;      
		case 2e3
			RateIndex = defs.SAMPLE_RATE_2KSPS;      
		case 1e3
			RateIndex = defs.SAMPLE_RATE_1KSPS;      
		otherwise
			disp('ERROR: sample rate not valid')
			RateIndex = [];
			err = 1;
	end
catch
	disp('ERROR: check defs file - probably picked a bad sample rate')
	RateIndex = [];
	err = 1;
end
end


%-------------------------
function [err,InputRangeIndex] = ConvertInputRange(err,defs,InputRangeNum,InputImped)
%  	AlazarDefs
try
	% check whether valid input range for card
	switch InputImped
		case 50
			if isempty(find(defs.inputranges_50 == InputRangeNum))
				err = 1;
				disp(['ERROR: invalid input range for 50 Ohm coupling'])
				disp(['accepted input ranges are:',num2str(defs.inputranges_50)])
			end
		case 1e6
			if isempty(find(defs.inputranges_1000000 == InputRangeNum))
				err = 1;
				disp(['ERROR: invalid input range for 1 MOhm coupling'])
				disp(['accepted input ranges are:',num2str(defs.inputranges_1000000)])
			end
	end
			
	switch(InputRangeNum)
		case 40
			InputRangeIndex = defs.INPUT_RANGE_PM_40_V;
		case 20
			InputRangeIndex = defs.INPUT_RANGE_PM_20_V;
		case 16
			InputRangeIndex = defs.INPUT_RANGE_PM_16_V;
		case 10
			InputRangeIndex = defs.INPUT_RANGE_PM_10_V;
		case 8
			InputRangeIndex = defs.INPUT_RANGE_PM_8_V;    
		case 5
			InputRangeIndex = defs.INPUT_RANGE_PM_5_V;  
		case 4
			InputRangeIndex = defs.INPUT_RANGE_PM_4_V;   
		case 2
			InputRangeIndex = defs.INPUT_RANGE_PM_2_V;  
		case 1
			InputRangeIndex = defs.INPUT_RANGE_PM_1_V;  
		case 800e-3
			InputRangeIndex = defs.INPUT_RANGE_PM_800_MV;
		case 500e-3
			InputRangeIndex = defs.INPUT_RANGE_PM_500_MV;
		case 400e-3
			InputRangeIndex = defs.INPUT_RANGE_PM_400_MV; 
		case 200e-3
			InputRangeIndex = defs.INPUT_RANGE_PM_200_MV; 
		case 100e-3
			InputRangeIndex = defs.INPUT_RANGE_PM_100_MV;  
		case 80e-3
			InputRangeIndex = defs.INPUT_RANGE_PM_80_MV;     
		case 50e-3
			InputRangeIndex = defs.INPUT_RANGE_PM_50_MV;   
		case 40e-3
			InputRangeIndex = defs.INPUT_RANGE_PM_40_MV;    
		case 20e-3
			InputRangeIndex = defs.INPUT_RANGE_PM_20_MV;      
		case 0
			InputRangeIndex = 0;      % channel not used
		otherwise
			disp('ERROR: input range not valid')
			InputRangeIndex = [];
			err = 1;
	end
catch
	disp('ERROR: check defs file - probably picked a bad input range')
	InputRangeIndex = [];
	err = 1;
end
end

%-------------------------
function [err,out] = ConvertTriggerInputRange(err,defs,in)
	switch(in)
		case 5
			out = defs.ETR_X1;
		case 1
			out = defs.ETR_DIV5;
		otherwise
			disp('ERROR: trigger input range not valid')
			err = 1;
	end
end


%-------------------------
function [err,out] = ConvertInputImped(err,defs,in)
	switch(in)
		case 1e6
			out = defs.IMPEDANCE_1M_OHM;
		case 600
			out = defs.IMPEDANCE_600_OHM;
		case 300
			out = defs.IMPEDANCE_300_OHM;
		case 75
			out = defs.IMPEDANCE_75_OHM;
		case 50
			out = defs.IMPEDANCE_50_OHM;
		otherwise
			disp('ERROR: input impedance not valid')
			err = 1;
	end
end


%-------------------------
function [err,out] = ConvertCoupling(err,defs,in)
	switch lower(in)
		case 'dc'
			out = defs.DC_COUPLING;
		case 'ac'
			out = defs.AC_COUPLING;
		otherwise
			disp('ERROR: input coupling not valid')
			err = 1;
	end
end

%-------------------------
function [err,out] = ConvertTrigEngOperation(err,defs,in)
% written for j or k only
	switch lower(in)
		case 'j'	% J only
			out = defs.TRIG_ENGINE_OP_J;
		case 'k'	% K only
			out = defs.TRIG_ENGINE_OP_K;
% 		case 'j or k'
% 			out = defs.TRIG_ENGINE_OP_J_OR_K;
% 		case 'j and k'
% 			out = defs.TRIG_ENGINE_OP_J_AND_K;
% 		case 'j xor k'
% 			out = defs.TRIG_ENGINE_OP_J_XOR_K;
% 		case 'j and not k'
% 			out = defs.TRIG_ENGINE_OP_J_AND_NOT_K;
% 		case 'not j and k'
% 			out = defs.TRIG_ENGINE_OP_NOT_J_AND_K;
		otherwise
			disp('ERROR: trigger eng operation not valid')
			err = 1;
	end
end

%-------------------------
function [err,out] = ConvertTrigEngSource(err,defs,in)
	switch lower(in)
		case 'a'
			out = defs.TRIG_CHAN_A;
		case 'b'
			out = defs.TRIG_CHAN_B;
		case 'c'
			out = defs.TRIG_CHAN_C;
		case 'd'
			out = defs.TRIG_CHAN_D;
		case {'ext','external'}
			out = defs.TRIG_EXTERNAL;
		case 'disable'
			out = defs.TRIG_DISABLE;
		otherwise
			disp('ERROR: trigger eng source not valid')
			err = 1;
	end
end

%-------------------------
function [err,out] = ConvertTrigEngSlope(err,defs,in)
	switch in
		case 1	% positive
			out = defs.TRIGGER_SLOPE_POSITIVE;
		case -1	% negative
			out = defs.TRIGGER_SLOPE_NEGATIVE;
		otherwise
			disp('ERROR: T.TrigEngSlope not valid (must be +1 or -1')
			err = 1;
	end
end

%-------------------------
function [err,out] = ConvertTrigEngLevel(err,level,range)
% level depends on range, so have to know
%	1. which channel triggering off
%	2. input range of triggering channel

	if level == -range
		out = 0;
	else
		out = round(128 + 127 * level / range);	
			% doesn't seem right, never get 0...
	end
end


%-------------------------
% function [ret,str] = AlazarErrorCheck(lib,retCode)
function err = AlazarErrorCheck(lib,retCode,func)
	errorText = calllib(lib,'AlazarErrorToText',retCode);

	if retCode == 512	% no error
		err = 0;
	else			% got an error
		err = 1;
        disp(['ERROR in ',func,': ',errorText])
	end

end

%------------------------------
function err = CheckBoardTypeId(defs,boardTypeId)
switch boardTypeId
	case defs.ATS850
		boardTypeText = 'ATS850';
	case defs.ATS310
		boardTypeText = 'ATS310';
	case defs.ATS330
		boardTypeText = 'ATS330';
	case defs.ATS855
		boardTypeText = 'ATS855';
	case defs.ATS315
		boardTypeText = 'ATS315';
	case defs.ATS335
		boardTypeText = 'ATS335';
	case defs.ATS460
		boardTypeText = 'ATS460';
	case defs.ATS860
		boardTypeText = 'ATS860';
	case defs.ATS660
		boardTypeText = 'ATS660';
	case defs.ATS9462
		boardTypeText = 'ATS9462';
	case defs.ATS9870
		boardTypeText = 'ATS9870';
	case defs.ATS9350
		boardTypeText = 'ATS9350';
	case defs.ATS9325
		boardTypeText = 'ATS9325';
	case defs.ATS9440
		boardTypeText = 'ATS9440';
	otherwise
		boardTypeText = '?';
end          

	% check whether user input 'card' matches actual card
	if ~strcmp(defs.card,boardTypeText)
		disp(['ERROR: defs.card says ',defs.card,' while card reports ',boardTypeText])
		err = 1;
	else
		err = 0;
	end
end



