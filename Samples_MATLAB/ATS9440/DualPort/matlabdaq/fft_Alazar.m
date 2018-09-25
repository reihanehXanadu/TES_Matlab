function [f,Vnoise] = fft_Alazar(filename,varargin)
% function [f,Vnoise] = fft_Alazar(filename,varargin)
% Usage: [f, Vnoise] = plot_fft_Alazar(filename, RecordCountPerChannel, RecordLength, fs, rangePM_V);
% Vnoise: V/sqrt(Hz),  f: Hz
%
% all filenames need associated .mat file with all settings
% 
% steps to convert from time trace to V/rtHz
%	FFT
%	abs
%	square
%	mean
%	only use 1st half of FFT data (positive freqs)
%	multiply by 2/SampleRate/RecordLength (no 2 for 0 freq pt)
%	sqrt to get from PSD to V/rtHz
% 
% maybe slight error from lack of windowing


% defaults
graph = 0;
maxsize = 2^23;	% matlab tends to crash if set larger
% get matlab & computer memory info
% [uV sV] = memory;


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

if S.numchannels2record > 1
	err('not written to handle more than 1 channel yet')
end

numsubsets = max(round(S.RecordCountPerChannel*S.Repeats*S.RecordLength/maxsize),1);
FFTmean1 = zeros(S.RecordLength,1);
TimeData = zeros(S.RecordLength,S.RecordCountPerChannel*S.Repeats/numsubsets);
for m = 1:numsubsets
    TimeData = load_Alazar_part(filename,S.RecordCountPerChannel*S.Repeats/numsubsets,(m-1)*S.RecordCountPerChannel*S.Repeats/numsubsets);
	FFTmean1 = FFTmean1 + mean((abs(fft(TimeData))).^2,2);
		% quadrature average - divide by total number after
end
FFTmean1 = FFTmean1/numsubsets;	% divide by total number to get mean

halfpts = S.RecordLength / 2;
f = S.SampleRate * (0:halfpts) / S.RecordLength;	% positive frequencies
FFTpos = FFTmean1(1:(halfpts+1));
clear FFTmean1

Vnoise = FFTpos * (2/S.SampleRate/S.RecordLength);	% don't understand this normalization
	% Brice figured this out though, so he's the expert
clear FFTpos
Vnoise(1) = Vnoise(1)/2;	% DC component doesn't need to be doubled
Vnoise = sqrt(Vnoise);


if graph == 1
	figure, 
	loglog(f, Vnoise);
	axis tight
 	ylim([min(Vnoise),max(Vnoise(2:end))])
	grid on;
	xlabel('Frequency [Hz]');  
	ylabel('Output Voltage Noise [V/\surd{Hz}]');
	title(strrep(filename, '_', '\_'));
end

