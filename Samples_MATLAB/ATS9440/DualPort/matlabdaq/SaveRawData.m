function err = SaveRawData(fid,filename,S,ch,T,P,p_rawdata,varargin)



% defaults
err = 0;
graph       = 0;
traditional = 1; % if 0 use Alazar write to disk function

%check for options
for i = 1:length(varargin),
	if ischar(varargin{i})	% varargin commands are always char
		switch lower(varargin{i})
			case 'graph'
				graph = 1;
            case 'stream_to_disk'
                traditional = 0;
		end
	else % not char - don't check for varargin
	end
end
    

if traditional == 1

%     numchannels = numel(fieldnames(ch));
    numpts = S.BytesPerSample * S.RecsPerChannelPerBuffer * S.numchannels2record * S.RecordLength;
%     numpts = S.RecsPerChannelPerBuffer * S.numchannels * S.RecordLength;

    % convert from uint8 to uint16. (reduce num pts by half)
    %	ex 76,127 -> 4C,7F -> 7F4C -> 32588/4 -> 8147
    data = typecast(uint8(p_rawdata.value(1:numpts)),'uint16');
%     data = p_rawdata.value(1:numpts);

    % save data from current buffer to file as uint16
    fwrite(fid,data,'uint16');    
elseif traditional == 0
    error(' stream to disk does not work')   
end

% save header file with each .daq file
save(filename,'S','ch','T','P')

if graph
	data = load_Alazar(filename,'graph');
end
