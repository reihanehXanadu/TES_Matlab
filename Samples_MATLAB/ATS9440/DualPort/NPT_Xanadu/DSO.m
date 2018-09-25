
function [N_Avg,N2_Avg,c] = DSO(S,data,varargin)

global Nh
clear F
clear scaled_F
data3D = double(reshape(data, S.numchannels2record, S.RecordLength, S.RecordsPerBuffer ));


NumberOfChannelsToProcess=S.numchannels2record-2;  
   
    load('photon_limits')
    
    for chns=1:NumberOfChannelsToProcess
        
data2D=squeeze(data3D(chns,:,:));
%data2D=squeeze(data3D(1,:,:));

 X=data2D';
 X=data2D(450:700,:)';
 X=X-repmat(mean(X),S.RecordsPerBuffer,1);
 %figure(1)
 %hold off
 % plot(data3D(2,:,1))
  %hold on
 %plot(data2D(:,:))
 
% [PC,Lambda]=eig(transpose(X)*X);
% [maxL,Index_maxL]=max(diag(Lambda));
% PC1=PC(:,Index_maxL);
 load('First_Principal_Component')
 F(chns,:)=sum(X.*repmat(PC1(:,chns)',S.RecordsPerBuffer,1),2);
 
 %figure(2)
% hold off
 %hist(F,100)
 

        load('scale_factors')
        scaled_F(chns,:)=polyval(p(:,chns),F(chns,:));
        scaled_Ff=floor(scaled_F);  %scaled Feature floor
        
% figure(3)
 hold off
 subplot(NumberOfChannelsToProcess,1,chns)
 %hist(scaled_F(chns,:),100)
 
        maxpnumber=size(PB,1)-1;
        
        
        N=zeros(NumberOfChannelsToProcess, maxpnumber);
        N2=zeros(NumberOfChannelsToProcess, maxpnumber);
        N_Average=zeros(NumberOfChannelsToProcess);
        N2_Average=zeros(NumberOfChannelsToProcess);
         
         
        for i=1:maxpnumber
            N(chns,i)= sum(F(chns,:)>PB(i,chns)& F(chns,:)<PB(i+1,chns));
            N2(chns,i)=sum(scaled_Ff(chns,:)==i-1);
%             qp(chns,(photon(chns,:)>PB(chns,i)& photon(chns,:)<PB(chns,i+1)))=i-1;
%             N(chns,i)=sum(qp(chns,:)==i-1);
        end
        
             N_Avg(chns,1)=sum(N(chns,:).*(0:maxpnumber-1)/sum(N(chns,:)));
             N2_Avg(chns,1)=sum(N2(chns,:).*(0:maxpnumber-1)/sum(N2(chns,:)));
        
        %varargout=N_Average;
        %varargout=N2_Average;
    end
    
      Nh=cat(1,Nh,N_Avg(:,1));
  %figure(4)
      channelstoshow=NumberOfChannelsToProcess;
    for chns=1:channelstoshow
      
            subplot(channelstoshow,1,chns)
            plot(Nh(:))
          %  title(['P', num2str(chns), '=', num2str(Nm(chns,end)),'x 10(-19)'],'FontSize',20)
            %subplot(channelstoshow,2,chns+1)
            %plot(Navg(chns,:))
            %title(['Pavg', num2str(chns), '=', num2str(Navg(chns,end)),'x 10(-19)'],'FontSize',20)
        end
    drawnow;
    
    
    p1=scaled_Ff(1,:);    
    p2=scaled_Ff(2,:);
    c=zeros(maxpnumber,maxpnumber);
    for i=1:maxpnumber-1
        for j=1:maxpnumber
            
            c(i,j)=(sum(p2(p1==i)==j)/sum(p1==i));
        end
    end
    
    
end   
    