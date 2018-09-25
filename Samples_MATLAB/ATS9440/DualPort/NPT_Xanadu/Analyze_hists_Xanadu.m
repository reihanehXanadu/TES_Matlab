
global Navg
clear p
hist(F,100)
fprintf('Click at the boundry of peaks. Press enter when finished.\n')
   %[PB,PBy]=ginput
    PB=reshape(PB,length(PB)/(S.numchannels2record-1),S.numchannels2record-1)';
   %save('photon_limits','PB' )
   maxpnumber=size(PB,2)+1;
   for chns=1:S.numchannels2record-1
   p(chns,:)=polyfit(PB(chns,1:size(PB,2)),(0:(size(PB,2)-1)),min(5,size(PB,2)-1));
     end
    save('scale_factors','p')
   %plot scaled_peaks to see the scale 
   for chns=1:S.numchannels2record-1
    scaled_F(chns,:)=polyval(p(chns,:),F(chns,:));
    subplot(1,1,chns)
    hist(scaled_F(chns,:),100)
   end
    
   
   
   numchan=S.numchannels2record-1;
   
    load('scale_factors')
    
    for chns=1:numchan
       
        scaled_F(chns,:)=polyval(p(chns,:),F(chns,:));
        scaled_Ff=floor(scaled_F);  %scaled Feature floor
        
        maxpnumber=size(PB,2)-1;
        
        for i=1:maxpnumber
            N(chns,i)= sum(F(chns,:)>PB(chns,i)& F(chns,:)<PB(chns,i+1));
            N2(chns,i)=sum(scaled_Ff(chns,:)==i-1);
%             qp(chns,(photon(chns,:)>PB(chns,i)& photon(chns,:)<PB(chns,i+1)))=i-1;
%             N(chns,i)=sum(qp(chns,:)==i-1);
        end
    end
    
    N_Average(chns,1)=sum(N.*(0:size(PB,2)-2))/sum(N)
    N2_Average(chns,1)=sum(N2.*(0:size(PB,2)-2))/sum(N)
    
    
    