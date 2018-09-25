
%%plotting Sep 15 data, power scaling


 N_avg_passband=[.17,.41,.82;.16,.38,.7];
 N_avg_Ring=[.16,.48,1.35;.20,.6,1.6];
 power=[400,460,1600];
 %powers are multiplied by these factors to calculate the power before the chip. 100 is for 99/1 splitter, 200 is to convert average power to the peak power, and 2 db loss for measurement fiber..001 is to convert  uw to mW. to 100*200*exp(.2)*.001
 power_passband=100*200*exp(.2)*.001*[2.400,5.400,11.300];
 
 for i=1:2
     subplot(2,1,i)
     hold off
  plot(power_passband,(N_avg_Ring(i,:)))
  hold on
  plot(power_passband,(N_avg_Ring(i,:)),'*')
  title(['TES',num2str(2*i)])
  plot(power_passband,(N_avg_passband(i,:)))
  plot(power_passband,(N_avg_passband(i,:)),'*')
  xlabel('peak Power right before the chip (mW)')
  ylabel('Average photon number')
  legend('Ring','Ring','Passband','Passband')
 end