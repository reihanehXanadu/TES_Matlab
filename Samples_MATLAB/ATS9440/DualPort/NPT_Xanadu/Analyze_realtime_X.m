function [err,P] = Analyze_realtime(err,first_call, data, P, S, T, ch, varargin)

% funtions so far:
% create_template
% histogram


create_template = 0;
histogram       = 0;
joint_histogram = 0;
DSO             = 0;
get_variance    = 0;
trace_variance  = 0;
prob_dist       = 0;
plot_3d         = 0;
HOM             = 0;


numpts          = 2 * S.RecsPerChannelPerBuffer * S.numchannels2record * S.RecordLength;

if ~isfield(P,'nametemplate') || isempty(P.nametemplate)
    P.nametemplate = 'template.mat';
end
% if ~isfield(P,'figH') || isempty(S.FigH)
%     S.FigH = figure;
% end

for m = 1:length(varargin{1}),
    if ischar(varargin{1}{m})	% varargin commands are always char
		switch lower(varargin{1}{m})
            case 'create_template'
                create_template = 1;
            case 'create_template_and_run'
                create_template = 1;
            case 'histogram'
                histogram = 1;
            case 'probability_distribution'
                prob_dist = 1;
            case 'joint_histogram'
                joint_histogram = 1;
            case 'dso'
                DSO = 1;
            case 'get_variance'
                get_variance = 1;
            case 'trace_variance'
                trace_variance = 1;
            case 'homdata'
                HOM = 1;
            case '3d'
                plot_3d = 1;
        end
    end
end


%size(double((typecast(uint8(data.Value(1:numpts)),'uint16'))))

% create template %
if create_template
    if first_call == 1
        template = zeros(S.RecordLength, S.numchannels2record);
        disp(sprintf('creating template ...'))		
	else
        load(P.nametemplate)
    end
    
    % when FIFO, structure: ABCDABCDABCD
    % otherwise AAAAABBBBBCCCCCDDDDD

	% read data from buffer
    if S.FIFO
        data3D = reshape(data, S.numchannels2record, S.RecordLength, S.RecordsPerBuffer );
		template = template + mean(data3D(1,:,:),3);
	else
       data3D = reshape(data,  S.RecordLength,S.numchannels2record, S.RecordsPerBuffer );
		template = template + mean(data3D(:,1,:),3);
    end
    
     t=template./repmat( abs(sum(template,2)),1,S.RecordLength);
     