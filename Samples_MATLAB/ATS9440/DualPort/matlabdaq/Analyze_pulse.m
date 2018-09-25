function [data,data1] = pulse_Alazar(filename,varargin)
% function [out] = pulse_Alazar(filename,varargin)
% 
% steps are:
%	- make cuts based on offsets
%	- average all remaining traces together to get average pulse (multiples)
%	-	fit average pulse to theory - use this as filter
%	- filter all traces (matched filter) to get histogram of peak heights
%	- pick out just single pulses
%	- average all singles to get correct average pulse
%	-	fit average pulse to theory - use this as filter
%	- filter all traces (matched filter) to get final histogram
% 
% cuts based on offsets make low pts 'slightly' lower
% cuts based on peak position are bad - gets rid of lots of zeroes
% matched filter based on data vs fit is basically same

warning off all

% defaults
graph = 0;
ipretrigger = 0;	% # of pts before pulse to average for offset
nbins = 200;	% # of histogram bins
ftype = 'alazar';
partial_load = 0;	% load all traces
verbose = 0;
filter_theory = 0;	% whether to use average pulse or fit to average pulse
					% for matched filter
					% 0 = aver pulse, 1 = theory fit to aver pulse
filtermode = 'filter_1s';	% what filter to use
		% 'filter_1s' - use average ones pulse for filter
		% 'filter_ns' - use average n pulse for n pulse filter
quicktest = 0;	% if quick = 1, then numsubsets = 1
maxsize2open = 2^23;	% matlab crashes if load more data than this
	% maybe can get this from matlab 'memory' command...
dojitter = 0;	% also analyze jitter
data = []; data1 = []; place = 1;	% goes with dostruct

%check for options
for m = 1:length(varargin),
	if ischar(varargin{m})	% varargin commands are always char
		switch lower(varargin{m})
			case 'graph'
				graph = 1;
			case 'pretrigger'
				ipretrigger = varargin{m+1};
			case 'tds'
				ftype = 'tds';
			case 'partial_load'
				partial_load = varargin{m+1};
			case 'verbose'
				verbose = 1;
			case 'theory'
				filter_theory = 1;
			case 'filtermode'
				filtermode = varargin{m+1};
			case 'quicktest'
				quicktest = 1;
			case 'jitter'
				dojitter = 1;
			case {'struct'}
				dostruct = 1;
				data = varargin{m+1};
				place = varargin{m+2};
		end
	else % not char - don't check for varargin
	end
end

% read settings from filename.mat
load([filename,'.mat'])
SampleRate = S.SampleRate;
RecordLength = S.RecordLength;
RecordCountPerChannel = S.RecordCountPerChannel;
Repeats = S.Repeats;
% if Repeats ~= 1
% 	error('not written for >1 repeats yet')
% end
if S.numchannels2record > 1
	error('not written for >1 channels')
end

filename = [filename,'.daq'];

t = (0:RecordLength-1) / SampleRate;

% option to just load subset of traces
if partial_load ~= 0
	RecordCountPerChannel = partial_load;
	Repeats = 1;
end

numsubsets = max(round(RecordCountPerChannel*Repeats*RecordLength/maxsize2open),1);
	% limited by matlab memory - number of partial loads to perform
RecordCountPerLoad = floor(RecordCountPerChannel*Repeats/numsubsets);
if quicktest == 1
	numsubsets = 2;
	RecordCountPerChannel = RecordCountPerLoad*numsubsets;
end
if verbose
	disp(['numsubsets ',num2str(numsubsets)])
	disp(['RecordCountPerLoad ',num2str(RecordCountPerLoad)])
end


%% 1st load - find all offsets & get average pulse from all pulses  
% also some optional cuts 
%-----------------------------------------------------
tic
tempdata = zeros(RecordLength,RecordCountPerLoad);
[offsets] = deal(zeros(1,RecordCountPerChannel*Repeats));
meantrace = 0; 
% [means,maxs,imaxs,mins,imins,offsets] = deal(zeros(1,RecordCountPerChannel));
A = 1;
for m = 1:numsubsets
% 	if verbose, disp(['load 1: ',num2str(m),' of ',num2str(numsubsets)]), end
	tempdata = A*load_Alazar_part(filename, RecordCountPerLoad, (m-1)*RecordCountPerLoad,ftype);


	% pretrigger
	if m == 1 && ipretrigger == 0	% select pretrigger pt if not specified
		% get ipretrigger from average pulse of 1st chunk of data
% 		size(tempdata)
		meandata = mean(tempdata,2);
		
		% check whether pos or neg - want to flip negative traces to make
		% positive
		if abs(max(meandata)) < abs(min(meandata))	% neg pulses
			A = -1;
			tempdata = A*tempdata;
		end

		meandata = abs(meandata + min(meandata));	% make sure not crossing zero
		[dummy,imaxdata] = max(meandata);
		[dummy,i_dmax] = max(diff(meandata));	% look for biggest jump in pulse
			% offset to make sure not crossing zero
		ipretrigger = i_dmax - round(i_dmax*0.10);	% play it safe and back off 10%
% 		%plot pretrigger selection
% 		figure; hold on
% 		plot(meandata(1:imaxdata))
% 		plotline(gca,'vert',ipretrigger)
% % 		title('select last pretrigger point')
% % 		ipretrigger = selectdata2('selectionmode','closest','pointer','crosshair')
% % 		close(gcf)
% 		title('auto selected pretrigger for average trace')
% 		drawnow
		clear meandata
	end
	i1 = 1 + (m-1) * RecordCountPerLoad;
	i2 = m * RecordCountPerLoad;
	
	% find offset for each trace
	offsets(i1:i2) = mean(tempdata(1:ipretrigger,:),1);	% mean of each pretrigger
	
	% offset all data to zero pretrigger level
	tempdata = tempdata - repmat(offsets(i1:i2),RecordLength,1);
	
	meantrace = meantrace + sum(tempdata,2);	% sum of all traces
	
% 	means(i1:i2) = mean(tempdata,1);	% mean of each trace
% 	sums(i1:i2) = sum(tempdata,1);	% sum of each trace
% 	[maxs(i1:i2),imaxs(i1:i2)] = max(tempdata,[],1);	% max of each trace
% 	[mins(i1:i2),imins(i1:i2)] = min(tempdata,[],1);	% min of each trace

end
% meantrace = meantrace/numtraces_keep;	% divide by total number to get mean
meantrace = meantrace/RecordCountPerChannel/Repeats;	% divide by total number to get mean
meantrace_norm = meantrace / sum(meantrace);	% normalize

if verbose, disp(['load 1 took ',num2str2(toc,1),' sec (find offsets, get average pulse)']), end


% pulse cuts based on offsets & peak position
% make pulse cuts - offsets and max indices are best cuts
CI.default = [];	% default to keep all traces
% CI.offsets = selectdata_cuts(offsets,'(pulse) offset','hist',nbins,'ylog');
% CI.max = selectdata_cuts(maxs,'(pulse) max','hist',nbins,'ylog');
% CI.maxi = selectdata_cuts(imaxs,'(pulse) max index','hist',nbins,'ylog');
% CI.sum = selectdata_cuts(sums,'(pulse) sum','hist',nbins,'ylog');
% CI.mean = selectdata_cuts(means,'(pulse) mean','hist',nbins,'ylog');

% turn pulse cuts into KeepIndexes
KeepVec = ones(1,RecordCountPerChannel*Repeats);
CutNames = fieldnames(CI);
for q = 1:length(CutNames)
	CutName = char(CutNames(q));
	KeepVec(CI.(CutName)) = 0;
end
KeepIndex = find(KeepVec);	% find all nonzero elements (keepers)
numtraces_keep = numel(KeepIndex);	% num traces to keep


% if filter_theory == 1	% fit mean pulse with theory
	[max_mean,imax_mean] = max(meantrace);
	tfall = 2e-6;
	trise = 100e-9;
	% ideal pulse = A * ( e^(-(t-t0)/tfall) - e^(-(t-t0)/trise) ) * step(t-t0)
	%	p = [A,tfall,trise,t0]
	pinit = [max_mean,tfall,trise,t(imax_mean-5)];	% pulse fit guesses - [A,t1,t2,t0]
	func_pulse = @(p,t) p(1) .* (exp(-(t-p(4))./p(2)) - exp(-(t-p(4))/p(3))) .* step(t-p(4));
	[pbest,deltap,ybest] = easyfit(t,meantrace,pinit,func_pulse,[],[],'np');
	if verbose, disp(['  load 1 fit took ',num2str(toc),' sec']), end
	% 	disp(pinit)
	% 	disp(pbest)

	if graph && filter_theory == 1
		figure; hold on; 
		plot(t,meantrace); 
% 		plot(t,func_pulse(pinit,t),'r'); 
		plot(t,ybest,'k')
		% plot(xthry2,'k')
% 		legend('average pulse','guess','fit')
		legend('average pulse',...
			['fit, rise ',num2str(pbest(3)*1e9),' ns, fall ',num2str(pbest(2)*1e6),' us'])
		title('average of all pulses and theory fit')
	end
% else	% don't fit mean pulse with theory
% 	if graph
% 		figure; hold on; 
% 		plot(t,meantrace); 
% 		title('average of all pulses')
% 	end
% end


%% 2nd load - filter all pulses with rough mean pulse
%-----------------------------------------------------
tic
tempdata = zeros(RecordLength,RecordCountPerLoad);
[Afilter] = deal(zeros(1,RecordCountPerChannel*Repeats));
for m = 1:numsubsets
% 	if verbose, disp(['load 3: ',num2str(m),' of ',num2str(numsubsets)]), end
	tempdata = A*load_Alazar_part(filename, RecordCountPerLoad, (m-1)*RecordCountPerLoad,ftype);

	i1 = 1 + (m-1) * RecordCountPerLoad;
	i2 = m * RecordCountPerLoad;
	
	% offset all data to zero pretrigger level
	tempdata = tempdata - repmat(offsets(i1:i2),RecordLength,1);	
	
	% apply rough matched filter (aver all pulses) to all traces
	if filter_theory == 0	% use mean trace for filter
	 	Afilter(i1:i2) = dot(tempdata,repmat(meantrace_norm',RecordCountPerLoad,1)');
	else	% use theory for filter
		Afilter(i1:i2) = dot(tempdata,repmat(ybest',RecordCountPerLoad,1)');
	end
end
if verbose, disp(['load 2 took ',num2str2(toc,1),' sec (filter with average pulse)']), end

% figure, hold on
% [yhist,xhist] = hist(Afilter,nbins);
% [yhist2,xhist2] = hist(Afilter(KeepIndex),nbins);
% [yhist3,xhist3] = hist(Afilter_thry,nbins);
% [yhist4,xhist4] = hist(Afilter_thry(KeepIndex),nbins);
% plot(xhist,yhist,'.-')
% plot(xhist2,yhist2,'.-r')
% plot(xhist3,yhist3,'.-g')
% plot(xhist4,yhist4,'.-k')
% legend('mean data, all pulses','mean data, cuts','mean fit, all pulses','mean fit, cuts')


% user picks 0s, 1s, 2s,...
% switch filtermode
% 	case 'filter_1s'	% use average ones as filter for all pulses
% 		numlines = 4;	% just find 0,1,2 peaks
% 	case 'filter_ns'	% use n filter for n pulse
% 		numlines = 20;	% find as many peaks as wanted
% end
[KeepIndices,nLabel] = select_nphotons(Afilter,nbins);
numfilters = numel(KeepIndices);
for m = 1:numfilters
	numtraces(m) = numel(KeepIndices{m});
	if verbose
		disp(['  num ',num2str(m-1),'s traces ',num2str(numtraces(m))])
	end
end


%% 3rd load - get average pulse(s) with a better filter
%	filter using either singles only or n-photon
%-----------------------------------------------------
tic
tempdata = zeros(RecordLength,RecordCountPerLoad);
meantrace = zeros(numfilters,RecordLength);
for m = 1:numsubsets
% 	if verbose, disp(['load 3: ',num2str(m),' of ',num2str(numsubsets)]), end
	tempdata = A*load_Alazar_part(filename, RecordCountPerLoad, (m-1)*RecordCountPerLoad,ftype);

	i1 = 1 + (m-1) * RecordCountPerLoad;
	i2 = m * RecordCountPerLoad;
	
	% offset all data to zero pretrigger level
	tempdata = tempdata - repmat(offsets(i1:i2),RecordLength,1);
		
	% find mean trace for specific photon number pulses, average 1s,
	% average 2s,...
	for p = 2:numfilters
		KeepIndex_subset = KeepIndices{p}([KeepIndices{p} >= i1 & KeepIndices{p} <= i2]);
		KeepIndex_subset = KeepIndex_subset - (m-1) * RecordCountPerLoad;
		if numel(KeepIndex_subset) > 1
			meantrace(p,:) = meantrace(p,:) + sum(tempdata(:,KeepIndex_subset),2)';
		elseif numel(KeepIndex_subset) == 1
			meantrace(p,:) = meantrace(p,:) + tempdata(:,KeepIndex_subset)';
		else	% no traces - do nothing
		end		
	end
end
meantrace(1,:) = ones(1,RecordLength) / RecordLength;	% filter for 0s
for m = 2:numfilters
	meantrace(m,:) = meantrace(m,:)/numtraces(m);	% divide by total number to get mean
end

if verbose, disp(['load 3 took ',num2str(toc),' sec (finding better filter(s))']), end


% if graph
% 	figure, %hold all
% 	plot(t,meantrace(2:end,:))
% 	switch filtermode
% 		case 'filter_1s'	
% 			legend('1s')
% 		case 'filter_ns'
% 			addnum2legend([1:numfilters-2],0,'','s')
% 	end
% 	title('average pulse filter (after normalization')
% end

% if filter_theory == 1 %&& numfilters == 2
% fit mean pulse (just 1s) with theory
	for m = 2:numfilters
		[max_mean,imax_mean] = max(meantrace(m,:));
		% ideal pulse = A * ( e^(-(t-t0)/tfall) - e^(-(t-t0)/trise) ) * step(t-t0)
		%	p = [A,tfall,trise,t0]
		pinit = [max_mean,pbest(2),pbest(3),pbest(4)];	% pulse fit guesses - [A,t1,t2,t0]
		func_pulse = @(p,t) p(1) .* (exp(-(t-p(4))./p(2)) - exp(-(t-p(4))/p(3))) .* step(t-p(4));
		[pbest_ns(m,:),deltap,ybest_ns(m,:)] = easyfit(t,meantrace(m,:),pinit,func_pulse,[],[],'np');
	end
% 	if graph
% 		figure; hold on; 
% 		plot(t,meantrace(2:end,:),'.'); 
% 		shiftcolororder(1);
% 		plot(t,ybest_ns(2:end,:))
% 		legend('average 1s pulse','guess','fit')
% 		legend('average pulse',...
% 			['fit, rise ',num2str(pbest_ns(m,3)*1e9),' ns, fall ',num2str(pbest_ns(m,2)*1e6,' us'])
% 		title('average of 1s and theory fit')
% 	end	
	if verbose, disp(['  load 3 with fit took ',num2str(toc),' sec']), end
% end


if graph
	figure, hold on
	plot(t,meantrace(2:end,:),'.')
% 	addnum2legend([1:numfilters-1],0,'','s')
% 	add2legend('pretrigger')
% 	addnum2legend([1:numfilters-1],0,'','s')
% 	if filter_theory == 0
% 		title('average pulses')
% 	else
		plot(t,ybest_ns(2:end,:))
		for m = 2:numfilters
			add2legend([num2str(m-1),'s fit, rise ',...
				num2str2(pbest_ns(m,3)*1e9,1),' ns, fall ',...
				num2str2(pbest_ns(m,2)*1e6,3),' us'])

		end
		title('average pulses and fits')
% 	end
	plotline(gca,'vert',t(ipretrigger))	% pretrigger line (for offsets)
end
meantrace_norm = meantrace ./ repmat(sum(meantrace,2),1,RecordLength);	% normalize


if dojitter
	% want mean pulse derivative for jitter
	if filter_theory == 0	% take smoothed derivative of data pulse
		pts_sm = round(RecordLength*0.01);	% set to 1% of total pts
		for m = 2:numfilters
			pulse_deriv(m,:) = smooth_deriv(t,meantrace(m,:),pts_sm);
			pulse_deriv(m,1:pts_sm) = 0;	% 1st pts bad
		end
	else	% take derivative of theory pulse
		% ideal pulse = A * ( e^(-(t-t0)/tfall) - e^(-(t-t0)/trise) ) * step(t-t0)
		% derivative = A * ( -1/tfall * e^(-(t-t0)/tfall) + 1/trise * e^(-(t-t0)/trise) ) * step(t-t0)
		%	p = [A,tfall,trise,t0]
		func_pulse_deriv = @(p,t) p(1) .* (-1./p(2).*exp(-(t-p(4))./p(2)) + 1./p(3).*exp(-(t-p(4))/p(3))) .* step(t-p(4));
		for m = 2:numfilters
			pulse_deriv(m,:) = func_pulse_deriv(pbest_ns(m,:),t);
		end
	end
% 	if graph
% 		% plot deriv for jitter filter
% 		figure; 
% 		switch filtermode
% 			case 'filter_1s'	% use average ones as filter for all pulses
% 				plot(t,pulse_deriv(2,:))
% 				legend('1s')
% 			case 'filter_ns'	% use average ones as filter for all pulses
% 				plot(t,pulse_deriv(2:end,:))
% 				addnum2legend([1:numfilters-1],0,'','s')
% 		end
% % 		plot(t,pulse_deriv(2:end,:))
% 		title('average pulse derivative (for jitter filter)')
% % 		addnum2legend([1:numfilters-2],0,'','s')
% 	end
end

%% 4th load (final) - filter all pulses with good filter
%	also filter just singles with derivative average for jitter
%-----------------------------------------------------
tic
tempdata = zeros(RecordLength,RecordCountPerLoad);
Aconv_jitter = zeros(1,RecordLength);
Afilter = zeros(1,RecordCountPerChannel*Repeats);
k(1:numfilters) = 1;
for m = 1:numsubsets
% 	if verbose, disp(['load 4: ',num2str(m),' of ',num2str(numsubsets)]), end
	tempdata = A*load_Alazar_part(filename, RecordCountPerLoad, (m-1)*RecordCountPerLoad,ftype);

	i1 = 1 + (m-1) * RecordCountPerLoad;
	i2 = m * RecordCountPerLoad;
	
	% offset all data to zero pretrigger level
	tempdata = tempdata - repmat(offsets(i1:i2),RecordLength,1);	
		
	switch filtermode
		case 'filter_1s'	% use average ones as filter for all pulses
			% apply 1s matched filter to all traces
			if filter_theory == 0	% use mean trace for filter
				Afilter(i1:i2) = dot(tempdata,repmat(meantrace_norm(2,:),RecordCountPerLoad,1)');
			else	% use theory for filter
				Afilter(i1:i2) = dot(tempdata,repmat(ybest_ns(2,:),RecordCountPerLoad,1)');
			end
		case 'filter_ns'	% use n filter for n pulse
			% have to figure out which filter to use for each pulse (find n bin)
			if filter_theory == 0	% use mean trace for filter
				for p = i1:i2
					Afilter(p) = dot(tempdata(:,p-i1+1),meantrace_norm(nLabel(p)+1,:)');
						% multiply by appropriate filter
						% nLabel is 1 for 1s, 2 for 2s, ...
				end
			else	% use theory for filter
				for p = i1:i2
					Afilter(p) = dot(tempdata(:,p-i1+1),ybest_ns(nLabel(p)+1,:)');
						% multiply by appropriate filter
						% nLabel is 1 for 1s, 2 for 2s, ...
				end
			end
% 			njitter = [1,2];	% peaks to look at jitter (1 = 1s, 2 = 2s,...)
% 			njitter = [1:round(numfilters/2)];	% peaks to look at jitter 
	end
	njitter = [1:round(numfilters/2)];	% peaks to look at jitter - don't do all
	
	if dojitter
		for m1 = njitter+1	% analyze jitter for only a few peaks
			% for each single, apply jitter filter, then threshold at zero
			KeepIndex_subset = KeepIndices{m1}([KeepIndices{m1} >= i1 & KeepIndices{m1} <= i2]);
			KeepIndex_subset = KeepIndex_subset - (m-1) * RecordCountPerLoad;
			for m2 = KeepIndex_subset	% look at pulses one at a time
				switch filtermode
					case 'filter_1s'	% use average ones as filter for all pulses
						% conv 1s derivative filter to each trace
						Aconv_jitter = fftshift(ifft(fft(tempdata(:,m2)).*conj(fft(pulse_deriv(2,:)'))));
					case 'filter_ns'	% use n filter for n pulse
						% conv ns derivative filter to each trace
						Aconv_jitter = fftshift(ifft(fft(tempdata(:,m2)).*conj(fft(pulse_deriv(m1,:)'))));
				end
				% gives a negative peak, then positive peak 
				% want to fit zero crossing between the two peaks
				% can't just find zero of whole trace b/c of noise
		% 		[dummy,imin_jitter] = min(Aconv_jitter(:,k));
		% 		[dummy,imax_jitter] = max(Aconv_jitter(:,k));
				% know the zero crossing is near middle of data, so start there
				% with really noisy data, cannot use global min and max
				% just use few pts away from middle
				ilow = RecordLength/2 - 5;
				ihigh = RecordLength/2 + 5;
				done1 = 0; done2 = 0;
				while done1+done2 ~= 2	% in case not far enough from middle, move pts further away
					if Aconv_jitter(ilow) >= 0	% want it to be negative
						ilow = ilow - 5;
					else
						done1 = 1;
					end
					if Aconv_jitter(ihigh) <= 0	% want it to be positive
						ihigh = ihigh + 5;
					else
						done2 = 1;
					end
				end
				try
			% 		izero_jitter(k) = interp1_ptout(Aconv_jitter(imin_jitter:imax_jitter,k),0) + imin_jitter - 1;
					izero_jitter{m1}(k(m1)) = interp1_ptout(Aconv_jitter(ilow:ihigh),0) + ilow - 1;
				catch
					Aconv_jitter(ilow:ihigh)
					figure; plot(Aconv_jitter(ilow:ihigh))
					error('jitter error');
				end
				k(m1) = k(m1) + 1;
			end
		end
	end
end
if verbose, disp(['load 4 took ',num2str(toc),' sec (filtered data)']), end

if dojitter
	% find jitter
	t_half = t(1:(RecordLength/2+1));
	t_conv = [fliplr(-t_half),t_half(2:end-1)];
	for m1 = njitter+1	% analyze jitter for only a few peaks
		trise = interp1_ptin(t_conv,izero_jitter{m1});
		[yhist_trise(m1,:),xhist_trise(m1,:)] = hist(trise,nbins);
		try		% fit gaussian to find FWHM
			[pgauss_trise(m1,:),ygauss_trise(m1,:)] = easyfitgauss(xhist_trise(m1,:),yhist_trise(m1,:));	
				% pbest = [A,mu,sigma]
			center_trise(m1) = pgauss_trise(m1,2);
			FWHM_trise(m1) = 2*sqrt(2*log(2))*pgauss_trise(m1,3);

			if verbose
				disp([num2str(m1-1),'s jitter FWHM ',num2str(FWHM_trise(m1)*1e9),' ns'])
			end
		catch
			disp([num2str(m1-1),'s jitter - gaussian fit error'])
			ygauss_trise(m1,:) = yhist_trise(m1,:)*0;
			center_trise(m1) = 0;
			FWHM_trise(m1) = 0;
		end
	end

	if graph
		figure; hold on
		plot(xhist_trise(njitter+1,:)'*1e9,yhist_trise(njitter+1,:)','.-')
		shiftcolororder(1);
		plot(xhist_trise(njitter+1,:)'*1e9,ygauss_trise(njitter+1,:)')
		xlabel('arrival time [ns]')
		title('jitter histogram')
		for m1 = njitter+1	% analyze jitter for only a few peaks
			add2legend([num2str(m1-1),'s, FWHM ',num2str(FWHM_trise(m1)*1e9),' ns'])
		end
	end
end
	

%% correct histogram of peak height/area
% also calculate resolution (FWHM of 1s peak) / (separation between 1s and 0s/2s)

% get final histogram using correct filter
[yhist_A,xhist_A] = hist(Afilter,nbins);

clear ybest pbest
for m = 1:numfilters
	xmin = min(Afilter(KeepIndices{m}));
	xmax = max(Afilter(KeepIndices{m}));
	is = find(xhist_A > xmin & xhist_A <= xmax);
	xpeak{m} = xhist_A(is); 
	ypeak{m} = yhist_A(is);
% 	ypeak{m}
	try
		% fit gaussian to find center & FWHM 
		[pbest(m,:),ybest{m}] = easyfitgauss(xpeak{m},ypeak{m});	% pbest = [A,mu,sigma]
		center(m) = pbest(m,2);
		FWHM(m) = 2*sqrt(2*log(2))*pbest(m,3);
	catch
		disp([num2str(m-1),'s peak - could not fit gaussian'])
		figure; plot(xpeak{m},ypeak{m})
		ybest{m} = ypeak{m}*0;
		center(m) = 0;
		FWHM(m) = 0;
	end
end
% save('temp','xpeak','ypeak')

if graph == 1
	figure, hold all
	plot(xhist_A',yhist_A','.-k')
	currylim = ylim;
	for m = 1:numfilters
		plot(xpeak{m},ybest{m})
	end
	set(gca,'yscale','Log')
	ylim([1,currylim(2)])	% fits mess up y axis
	legend('data')
	for m = 1:numfilters
		add2legend([num2str(m-1),'s fit, normFWHM ',num2str(FWHM(m)/FWHM(1))])
	end
	switch filtermode
		case 'filter_1s'	% use average ones as filter for all pulses
			title('pulse area/height histogram - filtered using 1s')
		case 'filter_ns'	% use average ones as filter for all pulses
			title('pulse area/height histogram - filtered using n-photon pulse')
	end
end

% if verbose
% 	disp(['FWHM of 0s peak ',num2str(FWHM_0s)])
% 	disp(['FWHM of 1s peak ',num2str(FWHM_1s)])
% 	disp(['FWHM of 2s peak ',num2str(FWHM_2s)])
% end

% energy resolution of 1s peak
dE = mean([center(3)-center(2),center(2)-center(1)]);
	% average of difference between 0s & 1s and between 1s & 2s
resolution = dE / FWHM(2);


xhist_A2 = xhist_A / center(2);	% scale x histogram by 1s position

% % output
% out.numfilters = numfilters;
% out.meantrace = meantrace;
% out.numtraces = numtraces;
% out.pulseheights = Afilter;
% out.xhist = xhist_A;
% out.yhist = yhist_A;
% out.center = center;
% out.FWHM = FWHM;
% out.dE = dE;
% out.resolution = resolution;
% out.t = t;
% % 
% if dojitter
% 	out.njitter = njitter;
% 	out.risetimes = trise;
% 	out.xhist_trise = xhist_trise;
% 	out.yhist_trise = yhist_trise;
% 	out.center_trise = center_trise;
% 	out.FWHM_trise = FWHM_trise;
% end

%---------------------------------------------
% add data to structure for output
%---------------------------------------------

data1 = [];
[data,data1] = add2structB(data,data1,'numfilters',numfilters,place);
[data,data1] = add2structB(data,data1,'meantrace',meantrace,place);
[data,data1] = add2structB(data,data1,'meantrace_norm',meantrace_norm,place);
[data,data1] = add2structB(data,data1,'numtraces',numtraces,place);
[data,data1] = add2structB(data,data1,'pulseheights',Afilter,place);
[data,data1] = add2structB(data,data1,'xhist',xhist_A,place);
[data,data1] = add2structB(data,data1,'xhist2',xhist_A2,place);
[data,data1] = add2structB(data,data1,'yhist',yhist_A,place);
[data,data1] = add2structB(data,data1,'center',center,place);
[data,data1] = add2structB(data,data1,'FWHM',FWHM,place);
[data,data1] = add2structB(data,data1,'dE',dE,place);
[data,data1] = add2structB(data,data1,'resolution',resolution,place);
[data,data1] = add2structB(data,data1,'t',t,place);
[data,data1] = add2structB(data,data1,'tfall',pbest_ns(:,2),place);
[data,data1] = add2structB(data,data1,'trise',pbest_ns(:,3),place);
[data,data1] = add2structB(data,data1,'meanfit',ybest_ns,place);
[data,data1] = add2structB(data,data1,'histfit',ybest,place);

if dojitter == 1
	[data,data1] = add2structB(data,data1,'njitter',njitter,place);
	[data,data1] = add2structB(data,data1,'risetimes',trise,place);
	[data,data1] = add2structB(data,data1,'xhist_trise',xhist_trise,place);
	[data,data1] = add2structB(data,data1,'yhist_trise',yhist_trise,place);
	[data,data1] = add2structB(data,data1,'center_trise',center_trise,place);
	[data,data1] = add2structB(data,data1,'FWHM_trise',FWHM_trise,place);
	[data,data1] = add2structB(data,data1,'jitterfit',ygauss_trise,place);
end


