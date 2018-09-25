function [KeepIndices,nLabel] = select_nphotons(ys,nbins,varargin)
% function ipeak_limits = select_nphotons(ys,nbins)
% 
% see if it's possible to pick out peaks and valleys automatically
% pass photon number resolved histogram
% 
% attempts to guess peaks and valleys - only works for large S/N
% user can drag vertical lines around to correct errors
% 
% options:
%	numlines,#		
% 
% needed subfunctions:
%	draggable, interp1_ptout, plotline

% defaults
graph = 0;
numlines = 20;	% default # of vertical lines to draw on histogram

%check for options
for m = 1:length(varargin),
	if ischar(varargin{m})	% varargin commands are always char
		switch lower(varargin{m})
			case 'graph'
				graph = 1;
			case 'numlines'
				numlines = varargin{m+1};
		end
	else % not char - don't check for varargin
	end
end


[yhist,xhist]=hist(ys,nbins); 

% peaklimit = 100;	% # counts to distinguish between peak & valley
peaklimit = max(yhist)/10;	% # counts to distinguish between peak & valley
[ipeak,ivalley] = deal([]);
lastpt = 1;	% 0 = valley, 1 = peak
% go through all pts one by one
try
	for m = 1:nbins
		if yhist(m) < peaklimit		% currently at valley
			currpt = 0;
			if lastpt == 1
				ivalley = [ivalley,m];
				lastpt = 0;
			else	% valley to valley (keep lowest)
				[dummy,imin] = min([yhist(ivalley(end)),yhist(m)]);
				if imin == 2	% current peak highest
					ivalley(end) = m;
				end
			end
		else						% currently at peak
			if lastpt == 0	% went from valley to peak
				ipeak = [ipeak,m];
				lastpt = 1;
			else	% peak to peak (keep highest)
				[dummy,imax] = max([yhist(ipeak(end)),yhist(m)]);
				if imax == 2	% current peak highest
					ipeak(end) = m;
				end
			end
		end
	end
catch
	ipeak = [];
	ivalley = [];
end
ivalley(1) = 1;	% make sure lowest line at lowest value

% plot histogram with lots of movable vertical lines
figH = figure; 
uibutton = uicontrol('Style','pushbutton','string','done',...
	'Callback','uiresume(gcbf)');
	% add 'done' button
hold on
plot(xhist,yhist)
plot(xhist(ivalley),yhist(ivalley),'ob')
plot(xhist(ipeak),yhist(ipeak),'or')
set(gca,'yscale','log')
title('move vertical lines to valley positions between peaks (extra lines @ right)')

for m = 1:min([numel(ivalley),numlines])
	valleylines(m) = line([xhist(ivalley(m)),xhist(ivalley(m))],ylim,...
			'linestyle','--','color','k');
	draggable(valleylines(m),'h');
end
xmax = max(xhist);

numlines = numlines - numel(ivalley);	% numlines remaining to draw
if numlines > 0
	for m = 1:numlines
		extravalleylines(m) = line([xmax,xmax],ylim,...
				'linestyle','--','color','k');
		draggable(extravalleylines(m),'h');
	end
else
	extravalleylines = [];
end

% wait for user to move lines around and then push button
uiwait(figH);
delete(uibutton);

% need to get data from plot lines
xlines = []; xvalleys = [];
alllines = [valleylines,extravalleylines];
for m = 1:numel(alllines)
	x1 = get(alllines(m),'XData');
		% x position of each line
	if x1 < xmax % only keep if not at far right of plot
		xlines = [xlines,x1(1)];
	end
end
xlines = sort(xlines);

% convert lines positions to histogram indices
clear ivalley
ivalley = round(interp1_ptout(xhist,xlines));
ivalley(ivalley<1) = [];	% get rid of negative

% clear draggable plot and replot lines
cla
hold on
plot(xhist,yhist)
% plot(xhist(ivalley),yhist(ivalley),'ob')

plotline(gca,'vert',xhist(ivalley))
set(gca,'yscale','log')
title('selected valley positions')

% convert positions to min/max for each peak
% then keepindices for each peak
% nLabel is 0 for zeros, 1 for ones,...
nLabel = ys*0+1;
for m = 1:numel(ivalley)-1
% 	ipeak_limits(m,1:2) = ivalley(m:m+1);
% 	ipeak_limits(m,1:2) = ivalley(m:m+1);
	ymin = xhist(ivalley(m));
	ymax = xhist(ivalley(m+1));
	KeepIndices{m} = find(ys > ymin & ys <= ymax);
	nLabel(KeepIndices{m}) = m-1;
% 	switch m
% 		case 1
% 			KeepIndices{m} = find(ys <= ymax);
% 		otherwise
% 			KeepIndices{m} = find(ys > ymin & ys <= ymax);
% 	end
end
KeepIndices{m+1} = find(ys > ymax);	% lump all higher n photons together
nLabel(KeepIndices{m+1}) = m;


close(figH)

