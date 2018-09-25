
function [varargout,F] = Calibrate_max(S,data,varargin)
clear F


data3D = double(reshape(data, S.numchannels2record, S.RecordLength, S.RecordsPerBuffer )); 

NumberOfChannelsToProcess=S.numchannels2record-2;
F=zeros(S.RecordsPerBuffer,NumberOfChannelsToProcess);
for i=1:S.numchannels2record-2
    
data2D=squeeze(data3D(i,:,:));

 X=data2D(100:280,:)';
 X=X-repmat(mean(X),S.RecordsPerBuffer,1);
hold off

 [PC,Lambda]=eig(transpose(X)*X);
 [maxL,Index_maxL]=max(diag(Lambda));
 PC1(:,i)=PC(:,Index_maxL);
% F(:,i)=sum(X.*repmat(PC1(:,i)',S.RecordsPerBuffer,1),2);
 F(:,i)=max(X');
 %F=

 clear p
 figure(8)
 hold off
 fprintf('Click at the boundry of peaks. Press enter when finished.\n')
 subplot(2,1,i)
 hist(F(:,i),1000);
 [PB(:,i),PBy(:,i)]=ginput;
 
end

save('First_Principal_Component','PC1')
   %load('photon_limits')
    PB=reshape(PB,length(PB(:))/(NumberOfChannelsToProcess),NumberOfChannelsToProcess);
   save('photon_limits','PB' )
   
   maxpnumber=size(PB,2)+1;
   for chns=1:NumberOfChannelsToProcess
      
   p(:,chns)=polyfit(PB(:,chns),(0:(size(PB,1)-1))',min(5,size(PB,1)-1));
  
   %plot scaled_peaks to see the scale 
   %for chns=1:NumberOfChannelsToProcess
   figure(9)
    scaled_F(:,chns)=polyval(p(:,chns),F(:,chns));
    subplot(2,1,chns)
    hist(scaled_F(:,chns),1000)
   end
  save('scale_factors','p')
end
 