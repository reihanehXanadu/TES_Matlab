

function [N_Avg,N2_Avg,N,N2,C11] = Cross_correlation(S,data,varargin)

global Nh
clear F
clear scaled_F
data3D = double(reshape(data, S.numchannels2record, S.RecordLength, S.RecordsPerBuffer ));


NumberOfChannelsToProcess=S.numchannels2record-2;  
   
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
 X=data2D(50:200,:)';
 X=X-repmat(mean(X),S.RecordsPerBuffer,1);
 
 load('First_Principal_Component')
 F(chns,:)=sum(X.*repmat(PC1(:,chns)',S.RecordsPerBuffer,1),2);
 
 
        load('scale_factors')
        scaled_F(chns,:)=polyval(p(:,chns),F(chns,:));
        scaled_Fc=ceil(scaled_F);  %scaled Feature floor
        
 figure(3)
 hold off
 subplot(NumberOfChannelsToProcess,1,chns)
 hist(scaled_F(chns,:),1000)
 
       
        
        
        
         
         
        for i=1:maxpnumber
            N(chns,i)= sum(F(chns,:)>PB(i,chns)& F(chns,:)<PB(i+1,chns));
            N2(chns,i)=sum(scaled_Fc(chns,:)==i-1);
%             qp(chns,(photon(chns,:)>PB(chns,i)& photon(chns,:)<PB(chns,i+1)))=i-1;
%             N(chns,i)=sum(qp(chns,:)==i-1);
        end
        
             N_Avg(chns,1)=sum(N(chns,:).*(0:maxpnumber-1)/sum(N(chns,:)));
             N2_Avg(chns,1)=sum(N2(chns,:).*(0:maxpnumber-1)/sum(N2(chns,:)));
        
        %varargout=N_Average;
        %varargout=N2_Average;
    end
    
      
      channelstoshow=NumberOfChannelsToProcess;
    
    
    
    
    
    
    for k=0:100  
    p1=circshift(scaled_Fc(1,:),k);    
    p2=scaled_Fc(2,:);  
    C11(k+1)=(sum(p2(p1==1)==1)/S.RecordsPerBuffer);   
    end
    
end   
    