

clear F


data3D = double(reshape(data, S.numchannels2record, S.RecordLength, S.RecordsPerBuffer )); 

NumberOfChannelsToProcess=S.numchannels2record-2;
F=zeros(S.RecordsPerBuffer,NumberOfChannelsToProcess);
%for i=1:S.numchannels2record-2
 
data2D=squeeze(data3D(1,:,:));
% Base=data2D(250:490,:)';
 X=data2D(100:300,:)';
 X=X-repmat(mean(X),S.RecordsPerBuffer,1);
hold off

 [PC,Lambda]=eig(transpose(X)*X);
 [maxL,Index_maxL]=max(diag(Lambda));
 PC1=PC(:,Index_maxL);
 F=sum(X.*repmat(PC1',S.RecordsPerBuffer,1),2);
 
 %F=

 clear p
 figure(3)
 hold off
 fprintf('Click at the boundry of peaks. Press enter when finished.\n')
 %subplot(2,1,i)
 hist(F,1000);
 %[PB(:,i),PBy(:,i)]=ginput;
 
%end

save('First_Principal_Component','PC1')
   load('photon_limits')
    PB=reshape(PB,length(PB(:))/(NumberOfChannelsToProcess),NumberOfChannelsToProcess);
   save('photon_limits','PB' )
   
%    maxpnumber=size(PB,2)+1;
%    for chns=1:NumberOfChannelsToProcess
%       % chns=1;
%    p(:,chns)=polyfit(PB(:,chns),(0:(size(PB,1)-1))',min(5,size(PB,1)-1));
%   
%    %plot scaled_peaks to see the scale 
%    %for chns=1:NumberOfChannelsToProcess
%     scaled_F(:,chns)=polyval(p(:,chns),F(:,chns));
%     subplot(2,1,chns)
%     hist(scaled_F(:,chns),1000)
%    end
%   save('scale_factors','p')
 