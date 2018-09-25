load('C:\Data\TES_Data\Sep20\Ring')


figure(8) 
hold off
%plot(sum(N_avg_ring),Variance_Delta_ring,'*')
hold on
%plot(sum(N_avg_ring),Variance_Delta_ring)
 f_ring=fit(sum(N_avg_ring)',Variance_Delta_ring','poly1');
hold on
plot(f_ring,sum(N_avg_ring),Variance_Delta_ring,'*' )
xlabel('<N_S>+<N_I>')
ylabel('Variance(N_S-N_I)')
legend(['slope=',num2str(f_ring.p1)])

%%
load('C:\Data\TES_Data\Sep20\Coh')


figure(8) 
%hold off
%plot(sum(N_avg),Variance_Delta,'*')
hold on

Variance_Delta_coh_sep21=[ 3.7849,4.2455,  4.8135];
N_avg_coh_sep21=[3.8697,4.3512,4.9455];
sumAll=[fliplr(N_avg_coh_sep21 ),sum(N_avg_coh)];
VarAll=[fliplr(Variance_Delta_coh_sep21 ),Variance_Delta_coh];
f_coh=fit(sumAll',VarAll','poly1');
plot(f_coh,sumAll,VarAll,'*' )
%plot(sum(N_avg),Variance_Delta)
xlabel('<N_S>+<N_I>')
ylabel('Variance(N_S-N_I)')
legend(['Ring slope=',num2str(f_ring.p1)],['Coherent state slope=',num2str(f_coh.p1)])

%%

load('Ring2')


figure(8) 

%plot(sum(N_avg_ring2),Variance_Delta_ring2,'*')
hold on
f_ring2=fit(sum(N_avg_ring2)',Variance_Delta_ring2','poly1');
plot(f_ring2,sum(N_avg_ring2),Variance_Delta_ring2,'*' )
%plot(sum(N_avg_ring2),Variance_Delta_ring2)
xlabel('<N_S-N_I>')
ylabel('Variance(N_S-N_I)')
%title('Ring')
legend(['slope=',num2str(f_ring2.p1)])
%legend('Light from Ring','Light from Ring','Coherent State','Coherent State')

%%
load('passband')


figure(8) 

%plot(sum(N_avg_passband),Variance_Delta_passband,'*')
hold on
f_pass=fit(sum(N_avg_passband)',Variance_Delta_passband','poly1');
plot(f_pass,sum(N_avg_passband),Variance_Delta_passband,'*' )
%plot(sum(N_avg_passband),Variance_Delta_passband)
xlabel('<N_S-N_I>')
ylabel('Variance(N_S-N_I)')
%title('Ring')

legend(['Ring1  slope=',num2str(f_ring.p1)],['Coherent state slope=',num2str(f_coh.p1)],['Ring2 slope=',num2str(f_ring2.p1)],['Pass band slope=',num2str(f_pass.p1)])
%legend(' Ring','Light from Ring','Coherent State','Coherent State','Light from Ring2','Light from Ring2','Passband','Passband')