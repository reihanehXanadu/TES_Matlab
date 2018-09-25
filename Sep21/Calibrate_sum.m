
function [varargout,F] = Calibrate_sum(S,data,varargin)
clear F


data3D = double(reshape(data, S.numchannels2record, S.RecordLength, S.RecordsPerBuffer )); 

NumberOfChannelsToProcess=S.numchannels2record-2;
F=zeros(S.RecordsPerBuffer,NumberOfChannelsToProcess);
for i=1:S.numchannels2record-2
    
data2D=squeeze(data3D(i,:,:));
% Base=data2D(520:550,:)';
 X=data2D(30:65,:)';
% X=X-repmat(mean(X),S.RecordsPerBuffer,1);
hold off

 %[PC,Lambda]=eig(transpose(X)*X);
 %[maxL,Index_maxL]=max(diag(Lambda));
 %PC1(:,i)=PC(:,Index_maxL);
 F(:,i)=sum(X,2);
 F_max(:,i)=max(X');
 
 %F=

 clear p
 figure(4)
 hold off
 fprintf('Click at the boundry of peaks. Press enter when finished.\n')
 subplot(2,1,i)
 hist(F(:,i),1000);
 title('Hist Sum')
 figure(5)
 subplot(2,1,i)
 hist(F_max(:,i),1000);
 title('Hist Max')
 %[PB(:,i),PBy(:,i)]=ginput;
 
end

end
 