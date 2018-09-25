clear N_avg N_avg_diff
clear Variance Variance_Delta
PB_whole=[];
for i=1:10;
   m=28+i;
filename=['C:\Data\TES_Data\Sep20\Easy_name_data\power_',num2str(m),'dB'];
fid = fopen([filename,'.daq'], 'rb');
data = fread(fid,[1,S.RecordsPerBuffer*S.numchannels2record*S.RecordLength],'uint16','ieee-le');


Calibrate(S,data);
%Calibrate_max(S,data);
%load('scale_factors')
%PB_whole=[PB_whole,PB];
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


N_avg_coh(:,i)=N_Avg;
N_avg_diff_coh(i)=N_Avg_diff;
%Prob_dist(:,i)=P_n(:);
Variance_coh(:,i)=Var';
Variance_Delta_coh(i)=Var_Delta;
end
%%

save('C:\Data\TES_Data\Sep20\Coh','N_avg_coh','N_avg_diff_coh','Variance_coh','Variance_Delta_coh')

%% Addidng extra data from sep 21 for coherent state to complete coherent state calibration

Variance_Delta_coh_sep21=[ 3.7849,4.2455,  4.8135]
N_avg_coh_sep21=[3.8697,4.3512,4.9455]

%% This part plot and fit the Variances

%save(['C:\Data\TES_Data\Variance_vs_n_'],'Variance','N_avg');
%load(['C:\Data\TES_Data\Variance_vs_n_'],'Variance','N_avg');
figure(7)
hold off

for i=1:2
subplot(2,1,i)
plot(Variance(i,:),N_avg(i,:),'*')
hold on
plot(Variance(i,:),N_avg(i,:))
%f=fit(Variance(i,:)',N_avg(i,:)','poly1');
hold on
%plot(f, Variance(i,:), N_avg(i,:))

xlabel('Average Photon Number');
ylabel('Variance')
title(['TES',num2str(2*i),'     slope=',num2str(f.p1)])
end
figure(8) 
plot(Variance_Delta,sum(N_avg),'*')
hold on
plot(Variance_Delta,sum(N_avg))
%%

for i=1:11;
    m=i+44;
load(['C:\Data\TES_Data\Probability_distributions_',num2str(m),'dB'])
g20(:,i)=sum(n.*(n-1).*P_n,2)./sum(n.*P_n,2).^2;
end

figure(9)
plot([n;n])