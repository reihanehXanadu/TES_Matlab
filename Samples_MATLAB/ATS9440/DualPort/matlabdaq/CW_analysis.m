clear all


load c:\daq\data\10-13-11_dual_TES_coincidence_detection_80MHz.mat
fn = 'c:\daq\data\10-13-11_dual_TES_coincidence_detection_80MHz.daq';
fid = fopen(fn);
data = fread(fid,[S.numchannels2record * S.RecordLength],'uint16','ieee-le');
fclose(fid);
chnames = 'ABCDEFGH';

P.threshold.A = 3200;
P.invert.A    = 1;
P.zero.A      = 10;

P.threshold.B = 500;
P.invert.B    = 0;
P.zero.B      = 10;
peak_width    = 100;


%%
S.RecsPerChannelPerBuffer = 1;
data3D = double(reshape(data, S.numchannels2record, S.RecordLength, S.RecsPerChannelPerBuffer))-2^15;

load(P.nametemplate)
template = template./repmat(abs(sum(template,1)),size(template,1),1);	% normalize


for chns = 1:S.numchannels2record
    
    switch S.numchannels2record
			case 1
				subplot(1,1,chns)
			case 2
				subplot(2,1,chns)
			case {3,4}
				subplot(2,2,chns)
    end
    
    if isfield(P.threshold,chnames(chns))

        % invert channels
        if isfield(P.invert,chnames(chns))
            if P.invert.(chnames(chns)) == 1 
                data3D(chns,:) = -data3D(chns,:); 
                template(:,chns) = -template(:,chns);
            end
        end
        
        % convolve trace with template
        data2D{chns} = conv(data3D(chns,:)-mean(data3D(chns,P.zero.(chnames(chns)))),template(:,chns) - mean(template(1:P.zero.(chnames(chns)),chns)) );
        
        % find peak indeces, pulse height and pulse position
        idxs   = find(data2D{chns} < P.threshold.(chnames(chns)));
        trgpos = idxs(find(diff(idxs)>1));
        for k = 1:size(trgpos,2)
            pkpos{chns}(k) = trgpos(k) + find(data2D{chns}(trgpos(k):trgpos(k)+peak_width) == max(data2D{chns}(trgpos(k):trgpos(k)+peak_width)));
        end
        pulseheight{chns}          = data2D{chns}(pkpos{chns});
        H                          = hist(pulseheight{chns},100);
        plot(H)
    end
end

%% get coincidence plot

for chns = 1:2
    coincPH{chns} = zeros(1,S.RecordLength);
    for m = 1:size(pkpos{chns},2)
        coincPH{chns}(pkpos{chns}(m)-80:pkpos{chns}(m)+80) = data2D{chns}(pkpos{chns}(m));
    end
    
    
    
    
end
plot(coincPH{1},coincPH{2},'*')