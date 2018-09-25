[N_Avg,N2_Avg,N,N2,c] = Analyze_real_time(S,data)


maxpnumber=size(N,2);
P_n=N./repmat(sum(N,2),1,maxpnumber)
n=[0:maxpnumber-1];

for i=1:S.numchannels2record-2
   
   Var(i)=std(n,P_n(i,:),2)^2;
   
end
figure(6)
plot(n,P_n,'*')
hold on
%n= n=[0:.1:maxpnumber-1];
plot(n,(exp(-N_Avg).*((N_Avg).^n))./factorial(n))
xlabel('n')
ylabel('P_n')