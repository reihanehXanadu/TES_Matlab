function [data,t] = load_Alazar_part(filename,nevents2load,skipevents,varargin)
% Usage: [data,t] = load_Alazar_part(filename,nevents2load,skipevents,varargin)
% 


% defaults
graph = 0;
maxsize2open = 2^24;	% matlab crashes if load more data than this
	% maybe can get this from matlab 'memory' command...
ftype = 'alazar';

%check for options
for i = 1:length(varargin),
	if ischar(varargin{i})	% varargin commands are always char
		switch lower(varargin{i})
			case 'graph'
				graph = 1;
			case 'tds'
				ftype = 'tds';
		end
	else % not char - don't check for varargin
	end
end


% stip .ext if exists
if strcmp(filename(end-3),'.')	
	filename = filename(1:end-4);
end

% read settings from filename.mat
load([filename,'.mat'])

if S.RecordLength * nevents2load > maxsize2open
	error('too much data to load for matlab to handle - reduce nevents2load')
end

%%
fid = fopen([filename,'.daq'], 'rb');
% 2*RecordLength*skipevents
switch ftype
	case 'alazar'
		status = fseek(fid,2*S.RecordLength*skipevents,'bof');	
			% x2 b/c data is 16 bit, which is 2 bytes
		if status ~= 0, error(['fseek error ',num2str(status)]), end
% 		data = fread(fid, [S.RecordLength nevents2load], 'uint16','ieee-le');
		if S.FIFO	% when FIFO, structure: ABCDABCDABCD
			data3D = reshape(double(fread(fid, S.RecordLength*nevents2load,'uint16','ieee-le')),...
				S.numchannels2record,S.RecordLength,nevents2load/S.numchannels2record);
		else	% when not FIFO, structure AAAAABBBBBCCCCCDDDDD
	%  		data3D = fread(fid, numpts,'uint16','ieee-le');
	% 		save('test','data3D')
			data3D = reshape(double(fread(fid, S.RecordLength*nevents2load,'uint16','ieee-le')),...
				S.RecordLength,S.numchannels2record,nevents2load/S.numchannels2record);
	% 		data3D = reshape(data3D,S.RecordLength,S.numchannels,S.RecordCountPerChannel*S.Repeats);
		end
	case 'tds'
		status = fseek(fid,4*S.RecordLength*skipevents,'bof');	
			% x4 b/c data is double (32 bit), which is 4 bytes
		if status ~= 0, error(['fseek error ',num2str(status)]), end
		data3D = fread(fid, [S.RecordLength nevents2load], 'double');
end
if (size(data3D,3)<nevents2load/S.numchannels2record)
	data3D = [];
	fclose(fid);
	error('not enough events to load')
	return
end
fclose(fid);

%% convert from uint16 to volts
switch ftype
	case 'alazar'
		% convert from uint16 to signed value from -1 to 1
		switch S.BitsPerSample
			case 14
				data3D = data3D/4;	
					% /4 converts from 16 bit to 14 bit
			case 16	% do nothing
		end
		data3D = (data3D - (2^(S.BitsPerSample-1) - 0.5)) / (2^(S.BitsPerSample-1) - 0.5);	
				% used to be (y14bit - 2^15) * range/2^15

		% figure out how many channels and which channels recorded			
		channelnames = fieldnames(ch);
		for m = 1:S.numchannels2record
			name = channelnames{S.channels2record(m)};
			if S.FIFO
				data3D(m,:,:) = data3D(m,:,:) * ch.(name).InputRange;
		% 		S.numchannels,S.RecordLength,S.RecordCountPerChannel*S.Repeats)
			else
				data3D(:,m,:) = data3D(:,m,:) * ch.(name).InputRange;
			end
		end

% 	case 'tds'	% already in volts
% 		Chselect = 1;
end

% %% some settings set to zero is changed in file (like complex impedance data)
% if SampleRate == 0, SampleRate = 1; end
% % if InputRangeA == 0, InputRangeA = 1; end
% % if InputRangeB == 0, InputRangeB = 1; end

% time
t = (0:S.RecordLength-1) / S.SampleRate;

%%
if graph
	for m = 1:S.numchannels2record
		figure
		if S.FIFO
			plot(t,squeeze(data3D(m,:,:)))
		else
			size(t)
			size(data3D)
			size(squeeze(data3D(:,m,:)))
			plot(t,squeeze(data3D(:,m,:)))
		end
		title(['Ch',channelnames{m}])
		xlabel('time')
		ylabel('V')
	end
end

data = squeeze(data3D);