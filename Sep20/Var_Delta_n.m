kik?%%
clear N_avg 
clear Variance
for i=1:2;
   
filename=['C:\Data\TES_Data\Sep20\Chip_ring_DCDC_',num2str(i)];

fid = fopen([filename,'.daq'], 'rb');
data = fread(fid,[1,S.RecordsPerBuffer*S.numchannels2record*S.RecordLength],'uint16','ieee-le');

Calibrate(S,data);
%Calibrate_max(S,data);
 [N,Photon_number,Photon_number_diff] = Analyze_diff(S,data);
 
maxpnumber=size(N,2);
P_n=N./repmat(sum(N,2),1,maxpnumber)
n=[0:maxpnumber-1];
N_Avg=mean(Photon_number,2);
N_Avg_diff=mean(Photon_number_diff);
%save(['C:\Data\TES_Data\Probability_distributions_',num2str(m),'dB'],'P_n','n')
figure(6)
%plot(n,P_n,'*')
hold off
plot(n,P_n,'*')
hold on
plot(n,(exp(-N_Avg).*((N_Avg).^n))./factorial(n))
xlabel('n')
ylabel('P_n')

Var_Delta=std(Photon_number_diff)^2;
Var= std(Photon_number').^2;


N_avg_ring(:,i)=N_Avg
N_avg_diff_ring(i)=N_Avg_diff;
%Prob_dist(:,i)=P_n(:);
Variance_ring(:,i)=Var';
Variance_Delta_ring(i)=Var_Delta
end

%%
save('C:\Data\TES_Data\Sep20\Ring','N_avg_ring','N_avg_diff_ring','Variance_ring','Variance_Delta_ring')


%% This part plot and fit the Variances

figure(8) 
hold off
%plot(sum(N_avg_ring),Variance_Delta_ring,'*')
hold on
%plot(sum(N_avg_ring),Variance_Delta_ring)
 f_ring=fit(sum(N_avg_ring)',Variance_Delta_ring','poly1');
hold on
plot(f_ring,sum(N_avg_ring),Variance_Delta_ring,'*' )
xlabel('<N_S-N_I>')
ylabel('Variance(N_S-N_I)')
legend(['slope=',num2str(f_ring.p1)])




%%

for i=1:11;
    m=i+44;
load(['C:\Data\TES_Data\Probability_distributions_',num2str(m),'dB'])
g20(:,i)=sum(n.*(n-1).*P_n,2)./sum(n.*P_n,2).^2;
end

figure(9)
plot([n;n])