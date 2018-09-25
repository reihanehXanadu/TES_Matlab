function [data3D,t] = load_Alazar(filename,varargin)
%
% all filenames need associated .mat file with all settings


% defaults
graph = 0;
maxsize2open = 2^24;	% matlab chokes if read in larger than that
	% maybe can get this from matlab 'memory' command...
	
%check for options
for i = 1:length(varargin),
	if ischar(varargin{i})	% varargin commands are always char
		switch lower(varargin{i})
			case 'graph'
				graph = 1;
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

fid = fopen([filename,'.daq'], 'rb');

% if S.RecordCountPerChannel*S.RecordLength > maxsize2open	% need to open in parts
numpts = S.RecordCountPerChannel*S.Repeats*S.RecordLength*S.numchannels2record;
% numpts = 2 * S.RecsPerChannelPerBuffer * S.numchannels * S.RecordLength;
if numpts > maxsize2open	
	% need to open in parts
	error('file too large to read in whole - try load_Alazar_part.m')
else
    if S.FIFO	% when FIFO, structure: ABCDABCDABCD
		data3D = reshape(double(fread(fid, numpts,'uint16','ieee-le')),...
			S.numchannels2record,S.RecordLength,S.RecordCountPerChannel*S.Repeats);
	else	% when not FIFO, structure AAAAABBBBBCCCCCDDDDD
%  		data3D = fread(fid, numpts,'uint16','ieee-le');
		data3D = reshape(double(fread(fid, numpts,'uint16','ieee-le')),...
			S.RecordLength,S.numchannels2record,S.RecordCountPerChannel*S.Repeats);
% 		data3D = reshape(data3D,S.RecordLength,S.numchannels,S.RecordCountPerChannel*S.Repeats);
	end
end

fclose(fid);

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


% time
t = (0:S.RecordLength-1) / S.SampleRate;

  
if graph
	for m = 1:S.numchannels2record
		figure
		if S.FIFO
			plot(t,squeeze(data3D(m,:,:)))
		else
% 			size(t)
% 			size(data3D)
% 			size(squeeze(data3D(:,m,:)))
			plot(t,squeeze(data3D(:,m,:)))
		end
		title(['Ch',channelnames{m}])
		xlabel('time')
		ylabel('V')
	end
end
