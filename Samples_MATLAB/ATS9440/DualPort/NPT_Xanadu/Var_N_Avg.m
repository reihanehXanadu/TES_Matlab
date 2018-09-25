
clear N_avg 
clear Variance
for i=1:11;
    m=i+44;
filename=['C:\Data\TES_Data\Ch2_4_',num2str(m),'dB_1MHz_1000Hz'];
fid = fopen([filename,'.daq'], 'rb');
data = fread(fid,'uint16','ieee-le');

Calibrate(S,data);
%Calibrate_max(S,data);
[N_Avg,N2_Avg,N,N2,c] = Analyze_real_time(S,data);

maxpnumber=size(N,2);
P_n=N./repmat(sum(N,2),1,maxpnumber)
n=[0:maxpnumber-1];
save(['C:\Data\TES_Data\Probability_distributions_',num2str(m),'dB'],'P_n','n')
figure(6)
plot(n,P_n,'*')
hold on
%n= n=[0:.1:maxpnumber-1];
plot(n,(exp(-N_Avg).*((N_Avg).^n))./factorial(n))
xlabel('n')
ylabel('P_n')


for j=1:S.numchannels2record-2
   
   Var(:,j)=std(n,P_n(j,:),2)^2;
   
end

N_avg(:,i)=N_Avg';
N2_avg(:,i)=N2_Avg';
%Prob_dist(:,i)=P_n(:);
Variance(:,i)=Var';

end

%% This part plot and fit the Variances

%save(['C:\Data\TES_Data\Variance_vs_n_'],'Variance','N_avg');
%load(['C:\Data\TES_Data\Variance_vs_n_'],'Variance','N_avg');
figure(7)
hold off

for i=1:2
subplot(2,1,i)
plot(Variance(i,:),N_avg(i,:),'*')
f=fit(Variance(i,:)',N_avg(i,:)','poly1');
hold on
plot(f, Variance(i,:), N_avg(i,:))
xlabel('Average Photon Number');
ylabel('Variance')
title(['TES',num2str(2*i),'     slope=',num2str(f.p1)])
end

%%

for i=1:11;
    m=i+44;
load(['C:\Data\TES_Data\Probability_distributions_',num2str(m),'dB'])
g20(:,i)=sum(n.*(n-1).*P_n,2)./sum(n.*P_n,2).^2;
end

figure(9)
plot([n;n])