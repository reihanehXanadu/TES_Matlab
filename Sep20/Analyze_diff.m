
function [N,Photon_number,Photon_number_diff] = Analyze_diff(S,data,varargin)


clear F
clear scaled_F

data3D = double(reshape(data, S.numchannels2record, S.RecordLength, S.RecordsPerBuffer*S.buffersPerAcquisition ));



NumberOfChannelsToProcess=S.numchannels2record-2;  
Photon_number=zeros(NumberOfChannelsToProcess,S.RecordsPerBuffer*S.buffersPerAcquisition); 
    load('photon_limits')
     maxpnumber=size(PB,1)-1;
        N=zeros(NumberOfChannelsToProcess, maxpnumber);
        N2=zeros(NumberOfChannelsToProcess, maxpnumber);
        N_Average=zeros(NumberOfChannelsToProcess);
        N2_Average=zeros(NumberOfChannelsToProcess);
        
    for chns=1:NumberOfChannelsToProcess
        
data2D=squeeze(data3D(chns,:,:));
%data2D=squeeze(data3D(1,:,:));

 X=data2D';
 X=data2D(100:250,:)';
 X=X-repmat(mean(X),S.RecordsPerBuffer*S.buffersPerAcquisition,1);

 load('First_Principal_Component')
 F(chns,:)=sum(X.*repmat(PC1(:,chns)',S.RecordsPerBuffer*S.buffersPerAcquisition,1),2);
 
 
        load('scale_factors')
        scaled_F(chns,:)=polyval(p(:,chns),F(chns,:));
        
 figure(3)
 hold off
 subplot(NumberOfChannelsToProcess,1,chns)
 hist(scaled_F(chns,:),1000)
  

    end
    scaled_Ff=floor(scaled_F);  %scaled Feature floor
    scaled_Ff(scaled_Ff<0)=0;
    Photon_number=scaled_Ff;
    Photon_number_diff=(Photon_number(1,:)-Photon_number(2,:));
  
    
    for i=1:maxpnumber
    N(:,i)=sum(Photon_number==i-1,2);
    %Delta_N(i)=sum(Photon_number_diff==i-1);
    end
    N_Avg=sum(N.*(0:maxpnumber-1),2)./sum(N,2);    
    N_Avg_diff=sum(Photon_number_diff)/size(Photon_number_diff,2);
end   
    