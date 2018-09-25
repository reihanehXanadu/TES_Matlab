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
        data3D = double(reshape(typecast(uint8(data.Value(1:numpts)),'uint16'), S.numchannels2record, S.RecordLength, S.RecsPerChannelPerBuffer))-2^15;
		template = template + mean(data3D,3)'/S.numBuffers/S.Repeats;
	else
        data3D = double(reshape(typecast(uint8(data.Value(1:numpts)),'uint16'), S.RecordLength, S.numchannels2record, S.RecsPerChannelPerBuffer))-2^15;
		template = template + mean(data3D,3)/S.numBuffers/S.Repeats;
	end
	template_norm = template./repmat(abs(sum(template,1)),S.RecordLength,1);

	% get ipretrigger from template - look for biggest jump in pulse
	[dummy,i_dmax] = max(diff(abs(template),1),[],1);
	% play it safe and back off 10%
	ipretrigger = i_dmax - round(i_dmax*0.10);
	
	
	% plot some traces of temp
	numtrace2plot = min([100,S.RecordCountPerChannel/S.numBuffers]);
	subplot(2,1,1)
	if S.FIFO
	 	plot(squeeze(data3D(end,:,1:numtrace2plot)))	% just plot 100 traces
		% ! not tested, might have to rotate
	else
	 	plot(squeeze(data3D(:,end,1:numtrace2plot)))	% just plot 100 traces
	end
%! problem - only plots last channel 100 traces
%!	not really made for plotting traces for multiple channels
	title(['last ',num2str(numtrace2plot),' traces'])
	
	% plot template and ipretrigger
	subplot(2,1,2), hold off
	plot(template), hold on
	for n = 1:S.numchannels2record
	  	plot([ipretrigger(n),ipretrigger(n)],[min(template_norm(:,n)),max(template_norm(:,n))],'--'); 
	end
	title('template')
    drawnow;
	
	% save template to .mat file
    save(P.nametemplate, '-mat', 'template')
end

if prob_dist
    
    chnames = 'ABCDEFGH';

    if first_call == 1
		disp(sprintf('analyzing histogram'))
	end
    load(P.nametemplate)
	template = template./repmat(abs(sum(template,1)),S.RecordLength,1);	% normalize

	% get ipretrigger from template - look for biggest jump in pulse
	[dummy,i_dmax] = max(diff(abs(template),1),[],1);
	% play it safe and back off 50%
	ipretrigger = i_dmax - round(i_dmax*0.50);
	

    [dim1,dim2] = size(template);
    
    % when FIFO, structure: ABCDABCDABCD
    % otherwise AAAAABBBBBCCCCCDDDDD

	% read data from buffer 

    
    if S.FIFO
        data3D = double(reshape(typecast(uint8(data.Value(1:numpts)),'uint16'), S.numchannels2record, S.RecordLength, S.RecsPerChannelPerBuffer))-2^15;
        if dim1 ~= size(data3D(:,:,1),2) || dim2 ~= size(data3D(:,:,1),1)
            disp('could not match template - create or load new template')\
            err = 1;
            return
        end
        for chns = 1:S.numchannels2record
            template_offset = mean(template(1:ipretrigger,chns));
            for recpc = 1:S.RecsPerChannelPerBuffer
               PH(recpc,chns) =  sum((data3D(chns,:,recpc)-template_offset)'.*((template(:,chns))/abs(sum(template(:,chns)))-template_offset)); % standard convolution
            end
        end
    else
        data3D = double(reshape(typecast(uint8(data.Value(1:numpts)),'uint16'), S.RecordLength, S.numchannels2record, S.RecsPerChannelPerBuffer))-2^15;
        if dim1 ~= size(data3D(:,:,1),1) || dim2 ~= size(data3D(:,:,1),2)
            disp('could not match template - create template')
            err = 1;
            return
        end
        for chns = 1:S.numchannels2record
            template_offset = mean(template(1:ipretrigger,chns));
            for recpc = 1:S.RecsPerChannelPerBuffer
               PH(recpc,chns) =  sum((data3D(:,chns,recpc)-template_offset).*((template(:,chns))-template_offset)); % standard convolution
            end
        end
    end
    load hist_vector.txt
    load photon_bounds.txt    
    pnum = zeros(length(PH'),2);
    nbar = zeros(1,2);
    g2 = zeros(1,2);
%     PH
    
    idx = zeros(S.numchannels2record,length(photon_bounds)+2);
    for chns = 1:S.numchannels2record

        if plot_3d
            figure(1)
        end
        if ~joint_histogram
            switch S.numchannels2record
                case 1
                    subplot(1,1,chns)
                case 2
                    subplot(2,1,chns)
                case {3,4}
                    subplot(2,2,chns)
            end
        end
        
        bins = linspace(hist_vector(chns,1),hist_vector(chns,2),P.hist_bins);
        bins3d(chns,:) = bins;
        disp(sprintf('channel: %d; min PH: %.2e; max PH: %.2e',chns, min(PH(:,chns)),max(PH(:,chns))))
        Y = hist(PH(:,chns),bins);
        if ~joint_histogram
            semilogy(bins,Y,'linewidth',3);
        end
        tempidx =[];
        tempidx(1) = 1;
        hold on
        for i = 1:length(photon_bounds)
            if ~joint_histogram
                plot([photon_bounds(chns,i) photon_bounds(chns,i)], [min(Y)+1 max(Y)],'r','linewidth',3)
            end
            tmp = find(bins < photon_bounds(chns,i));
            tempidx(i+1) = tmp(end);
        end
        tempidx(end+1) = length(bins);
        idx(chns,:) = tempidx;
        hold off
       
        [np,p] = get_probabilities(idx(chns,:),Y);
        for i=1:length(PH')
            pnum(i,chns) = length(photon_bounds(chns,:))-length(find(photon_bounds(chns,:)>PH(i,chns)));
        end
        %Nmean  = sum(np.*p);
        nbar(chns) = mean(pnum(:,chns));
        nvar(chns) = var(pnum(:,chns));
        g2(chns) = 1+((var(pnum(:,chns))-nbar(chns))/(nbar(chns)^2));

         formatSpec = 'nbar: %.3f, g2: %.3f';
         if ~joint_histogram
             title(sprintf(formatSpec,nbar(chns),g2(chns)),'fontsize',30);  
         end
    end
    pnummean = mean(pnum(:,1))+mean(pnum(:,2));
    varpnum = var(pnum(:,1)-pnum(:,2));

    g2tot = 1+((varpnum-pnummean)/(pnummean^2));
    hold on
    text(5000,1000,['g2:' num2str(g2tot)],'fontsize',20); 
    hold off

    
    pnummean-varpnum;
    drawnow
    if plot_3d


        figure(3)
 %       Y3 = hist3(PH,[length(bins),length(bins)]);
        edges{1} = bins3d(1,:);
        edges{2} = bins3d(2,:);
        Y3 = hist3(PH,'Edges',edges);
        pcolor(log10(Y3));
        shading flat;
              % axis([min(hist_vector(:,1)) max(hist_vector(1,:)) min(hist_vector(2,:)) max(hist_vector(2,:))]);
        Int_3d = zeros(length(photon_bounds),length(photon_bounds));
        
        hold off
       
        [np3,p3] = get_probabilities3d(idx,Y3);
        Nmean3  = sum(sum(np3.*p3));
        figure(4)
        p3(1,1) = 0;
        bar3(p3)
        
        eta1 = p3(2,2)/(p3(2,2)+p3(2,1));
        eta2 = p3(2,2)/(p3(2,2)+p3(1,2));
        eta1HOM = 2/(2+(p3(2,1)/p3(3,1)));
        eta2HOM = 2/(2+(p3(1,2)/p3(1,3)));
        drawnow;
        title(sprintf('eta1/2: %.3f/%.3f; HOM1/2: %.3f/%.3f',eta1,eta2,eta1HOM,eta2HOM),'fontsize',20)
        drawnow
        
        figure(2)
        % create 2d matrix for ChA and ChB
        plot(PH(:,1),PH(:,2),'.')
        hold on
        for i = 1:length(photon_bounds)
            plot([photon_bounds(1,i) photon_bounds(1,i)], [min(PH(:,2)) max(PH(:,2))],'r','linewidth',3)
            plot([min(PH(:,1)) max(PH(:,1))], [photon_bounds(2,i) photon_bounds(2,i)],'r','linewidth',3)
        end
        hold off
        title(sprintf('0/1-1/0-1/1: %d-%d-%d',round(p3(1,2)*S.RecordCountPerChannel),round(p3(2,1)*S.RecordCountPerChannel),round(p3(2,2)*S.RecordCountPerChannel)),'fontsize',20)
        drawnow;
        
    end
    
    if HOM
        stageposn = [];
        while isempty(stageposn)
            load stageposn.txt 
        end
        StgPsnTemp = stageposn(1);
        FinPsn = stageposn(3);
        fn = (['C:\daq\Data\130304_WFH_Vaccum_PH_HOMdataStgPsn' num2str(StgPsnTemp,'%.6f') '.mat']);
        fn2 = (['C:\daq\Data\130304_WFH_Vaccum_NBAR_HOMdataStgPsn' num2str(StgPsnTemp,'%.6f') '.mat']);
        if exist(fn) == 2 
            load([fn])
            load([fn2])
        else
            PHtemp2 = [];
            nbar2 = [];
        end
        PHtemp2 = [PHtemp2 PH];
        nbar2 = [nbar2 nbar];
        save(fn,'PHtemp2')
        save(fn2,'nbar2')
     end
    
    
     if joint_histogram
        
        % create 2d matrix for ChA and ChB
        plot(PH(:,1),PH(:,2),'.')
        hold on
        for i = 1:length(photon_bounds)
            plot([photon_bounds(1,i) photon_bounds(1,i)], [min(PH(:,2)) max(PH(:,2))],'r','linewidth',3)
            plot([min(PH(:,1)) max(PH(:,1))], [photon_bounds(2,i) photon_bounds(2,i)],'r','linewidth',3)
        end
        hold off        
        edges{1} = bins3d(1,:);
        edges{2} = bins3d(2,:);
        Y3 = hist3(PH,'Edges',edges);
        [np3,p3] = get_probabilities3d(idx,Y3);
        Nmean3  = sum(sum(np3.*p3));
        p3(1,1) = 0;
        bar3(p3)
        eta1 = p3(2,2)/(p3(2,2)+p3(2,1));
        eta2 = p3(2,2)/(p3(2,2)+p3(1,2));
        drawnow;
        title(sprintf('eta1/2: %.3f/%.3f',eta1,eta2),'fontsize',28)
     end
     
     
end

if histogram
    
    
    chnames = 'ABCDEFGH';
    
    if first_call == 1
		disp(sprintf('analyzing histogram'))
	end
    load(P.nametemplate)
	template = template./repmat(abs(sum(template,1)),S.RecordLength,1);	% normalize

	% get ipretrigger from template - look for biggest jump in pulse
	[dummy,i_dmax] = max(diff(abs(template),1),[],1);
	% play it safe and back off 50%
	ipretrigger = i_dmax - round(i_dmax*0.50);
	

    [dim1,dim2] = size(template);
    
    % when FIFO, structure: ABCDABCDABCD
    % otherwise AAAAABBBBBCCCCCDDDDD
    
	% read data from buffer
    if S.FIFO
        data3D = double(reshape(typecast(uint8(data.Value(1:numpts)),'uint16'), S.numchannels2record, S.RecordLength, S.RecsPerChannelPerBuffer))-2^15;
        if dim1 ~= size(data3D(:,:,1),2) || dim2 ~= size(data3D(:,:,1),1)
            disp('could not match template - create or load new template')
            err = 1;
            return
        end
        for chns = 1:S.numchannels2record
            template_offset = mean(template(1:ipretrigger,chns));
            for recpc = 1:S.RecsPerChannelPerBuffer
               PH(recpc,chns) =  sum((data3D(chns,:,recpc)-template_offset)'.*((template(:,chns))/abs(sum(template(:,chns)))-template_offset)); % standard convolution
            end
        end
    else
        data3D = double(reshape(typecast(uint8(data.Value(1:numpts)),'uint16'), S.RecordLength, S.numchannels2record, S.RecsPerChannelPerBuffer))-2^15;
        if dim1 ~= size(data3D(:,:,1),1) || dim2 ~= size(data3D(:,:,1),2)
            disp('could not match template - create template')
            err = 1;
            return
        end
        for chns = 1:S.numchannels2record
            template_offset = mean(template(1:ipretrigger,chns));
            for recpc = 1:S.RecsPerChannelPerBuffer
               PH(recpc,chns) =  sum((data3D(:,chns,recpc)-template_offset).*((template(:,chns))-template_offset)); % standard convolution
            end
        end
    end
    for chns = 1:S.numchannels2record
        
		switch S.numchannels2record
			case 1
				subplot(1,1,chns)
			case 2
				subplot(2,1,chns)
			case {3,4}
				subplot(2,2,chns)
        end

     
% -----------------------------------------------
% 		% draw draggable vertical lines
% 		if first_call == 1
% % 			figH = figure; 
% 			set(get(gca,'parent'),'DoubleBuffer','on');	% makes plot smoother
% 			
% 			% plot hist, don't specify x-values
% 			hold off
% 	        [yhist,xhist] = hist(PH(:,chns),P.hist_bins);
% 			P.hLine.(chnames(chns)) = plot(xhist,yhist,'k');		
% 			hold on
% 			set(gca,'yscale','log')
% 			
% 			xlim1 = xlim;
% 			xdiff = diff(xlim1);
% 			% draw 6 vlines, 1=xmin, 2=xmax,3=0s-1s,4=1s-2s,5=2smax
% 			% need initial xposition for vlines 
% 			%	set max/min lines at max/min
% 			%	set other lines at guesses - user fixes
% 			xpos = [xlim1(1),xlim1(2),...
% 				.2*xdiff + xlim1(1),.35*xdiff + xlim1(1),...
% 				.5*xdiff + xlim1(1)];
% 			% extend xlim a bit past limit lines
% 			xlim([xlim1(1)-0.05*xdiff,xlim1(2)+0.05*xdiff])
% 			lcolor = {'r','r','b','b','b'};
% 			for m = 1:numel(xpos)	% draw vertical lines
% 				vlines(m) = line([xpos(m),xpos(m)],ylim,...
% 						'linestyle','--','color',lcolor{m});
% 			end
% 			% make vlines draggable - only horizontally
% 			for m = 1:numel(vlines)
% 				draggable(vlines(m),'h');
% 			end	
% 			P.vlines.(chnames(chns)) = vlines;	% add to struct
% % 			xpos = [xmin,x2,x3,xmax];
%         
%         else
%    
%   			for m = 1:numel(P.vlines.(chnames(chns)))
% 				xpos2 = get(P.vlines.(chnames(chns))(m),'XData');	% 2 equal values since vertical
% 				xpos(m) = xpos2(1);	% just keep one value
%             end
%             
% 			% set xlim based on red lines
% 			xdiff = xpos(2) - xpos(1);
% 			% extend xlim a bit past limit lines
% 			xlim([xpos(1)-0.05*xdiff,xpos(2)+0.05*xdiff])
% 			
% 			% update hist with specified xlim
% 			xhist = [xpos(1):(xdiff+1)/P.hist_bins:xpos(2)];
% 	        yhist = hist(PH(:,chns),xhist);
% 			% updata plot with new data (better than replotting)
% 			set(P.hLine.(chnames(chns)),'XData',xhist) 
% 			set(P.hLine.(chnames(chns)),'YData',yhist) 
% 			drawnow                             % refresh display
% 
% 			% make sure 0smin is on screen
% 			if xpos(3) < xpos(1)
% 				set(P.vlines(3).(chnames(chns)),'XData',[xpos(1)+0.05*xdiff,xpos(1)+0.05*xdiff]);
% 			end
% 
% 			% make sure peak lines in correct order
% % 			xpos(3:5) = sort(xpos(3:5));
% % 
% % 			% find number of 0s counts
% % 			i0s = find(xhist < xpos(3));
% % 			num0s = sum(yhist(i0s));
% % 			
% % 			% find number of 1s counts
% % 			i1s = find(xpos(3) < xhist & xhist < xpos(4));
% % 			num1s = sum(yhist(i1s));
% % 			
% % 			% find number of 2s counts
% % 			i2s = find(xpos(4) < xhist & xhist < xpos(5));
% % 			num2s = sum(yhist(i2s));
% % 			
% % 			% use number of 0s & 1s to find mu01 and dmu01
% % 			mu01 = num1s / num0s;
% % 			dmu01 = mu01 * sqrt((sqrt(num0s)/num0s)^2 + (sqrt(num1s)/num1s)^2);
% % 			
% % 			% use number of 1s & 2s to find mu12 and dmu12
% % 			mu12 = 2 * num2s / num1s;
% % 			dmu12 = 2 * mu12 * sqrt((sqrt(num1s)/num1s)^2 + (sqrt(num2s)/num2s)^2);
% % 			
% % 			title(['0s=',num2str(num0s),', 1s=',num2str(num1s),', 2s=',num2str(num2s),...
% % 				', \mu_{01}=',num2str2(mu01,2),'+-',num2str2(dmu01,2),...
% % 				', \mu_{12}=',num2str2(mu12,2),'+-',num2str2(dmu12,2)]);
%  		end

% ---------------------------------------------------------------------------

        Y = hist(PH(:,chns),P.hist_bins);
        plot(Y);
        drawnow
        %sum(hist(PH(:,chns)));
    end
     

end

% DSO %
if DSO
    
    % when FIFO, structure: ABCDABCDABCD
    % otherwise AAAAABBBBBCCCCCDDDDD

	% read data from buffer
    chname = 'ABCD';
    if S.FIFO
        data3D = double(reshape(typecast(uint8(data.Value(1:numpts)),'uint16'), S.numchannels2record, S.RecordLength, S.RecsPerChannelPerBuffer))-2^15;
        if strcmp(T.Clock,'internal')
            plot((1:S.RecordLength)/S.SampleRate,data3D(:,:,1)')
            %axis([0 S.RecordLength/S.SampleRate ylim])
            axis([0 S.RecordLength/S.SampleRate -2^15 2^15])
        else
            plot(data3D(:,:,1)')
            
        end
        drawnow;
	else
        data3D = double(reshape(typecast(uint8(data.Value(1:numpts)),'uint16'), S.RecordLength, S.numchannels2record, S.RecsPerChannelPerBuffer))-2^15;
		template = template + mean(data3D,3);
        if strcmp(T.Clock,'internal')
            plot((1:S.RecordLength)/S.SampleRate,data3D(:,:,1)')
            %axis([0 S.RecordLength/S.SampleRate ylim])
            axis([0 S.RecordLength/S.SampleRate -2^15 2^15])
        else
            plot(data3d(:,:,1)')
            
        end
        drawnow;
	end
end


if trace_variance
    
    bins = 2^3;
    data3D     = double(reshape(typecast(uint8(data.Value(1:numpts)),'uint16'), S.numchannels2record, S.RecordLength, S.RecsPerChannelPerBuffer))-2^15;
    for i = 1:S.numchannels2record
%         for k = 1:S.RecordLength/bins
%             for m = 1:S.RecsPerChannelPerBuffer
%                 %smooth(i,(k-1)*bins+1:k*bins,m) = mean(data3D(i,(k-1)*bins+1:k*bins,m));
%             end
%         end
        %variance(i)   = var(reshape(data3D(i,:,:),1,S.RecordLength * S.RecsPerChannelPerBuffer)-reshape(smooth(i,:,:),1,S.RecordLength * S.RecsPerChannelPerBuffer));
        variance(i)   = var(reshape(data3D(i,:,:),1,S.RecordLength * S.RecsPerChannelPerBuffer));
    end
    %plot(data3D(1,:,1)'-smooth(1,:,1)')
    plot(data3D(1,1:200,1)')
    title(sprintf('Variance A: %.3e', variance(1)),'fontsize',22);
    drawnow;
    
    
end

if get_variance % dont know what this was supposed to be doing
    
    chnames = 'ABCDEFGH';
    
    if first_call == 1
		disp(sprintf('analyzing variance'))
	end
    load(P.nametemplate)
	template = template./repmat(abs(sum(template,1)),S.RecordLength,1);	% normalize

	% get ipretrigger from template - look for biggest jump in pulse
	[dummy,i_dmax] = max(diff(abs(template),1),[],1);
	% play it safe and back off 10%
	ipretrigger = i_dmax - round(i_dmax*0.10);
    ipretrigger = 25;
	

    [dim1,dim2] = size(template);
    
    % when FIFO, structure: ABCDABCDABCD
    % otherwise AAAAABBBBBCCCCCDDDDD
    
	% read data from buffer
    if S.FIFO
        data3D = double(reshape(typecast(uint8(data.Value(1:numpts)),'uint16'), S.numchannels2record, S.RecordLength, S.RecsPerChannelPerBuffer))-2^15;
        if dim1 ~= size(data3D(:,:,1),2) || dim2 ~= size(data3D(:,:,1),1)
            disp('could not match template - create or load new template')
            err = 1;
            return
        end
        for chns = 1:S.numchannels2record
            template_offset = mean(template(1:ipretrigger,chns));
            for recpc = 1:S.RecsPerChannelPerBuffer
               PH(recpc,chns) =  sum((data3D(chns,:,recpc)-template_offset)'.*((template(:,chns))/abs(sum(template(:,chns)))-template_offset)); % standard convolution
            end
        end
        plot(data3D(1,:,1)')
        axis([xlim -2^15 2^15])
        title(sprintf('variance: %.2f',var(PH)),'fontsize',22)
        drawnow;
    else
        data3D = double(reshape(typecast(uint8(data.Value(1:numpts)),'uint16'), S.RecordLength, S.numchannels2record, S.RecsPerChannelPerBuffer))-2^15;
        if dim1 ~= size(data3D(:,:,1),1) || dim2 ~= size(data3D(:,:,1),2)
            disp('could not match template - create template')
            err = 1;
            return
        end
        for chns = 1:S.numchannels2record
            template_offset = mean(template(1:ipretrigger,chns));
            for recpc = 1:S.RecsPerChannelPerBuffer
               PH(recpc,chns) =  sum((data3D(:,chns,recpc)-template_offset).*((template(:,chns))-template_offset)); % standard convolution
            end
        end
    end
    

end

