filename1=['C:\Data\Rings_scans\Die8_120um_450nm_Gap_1500_1600.csv'];

data1=csvread(filename1,16,0);

figure(3)
hold off
plot(data1(:,1),data1(:,2),'red')




%%
filename3=['C:\Data\Ming\Ming3_1500_1600_1pm_loss (1).csv'];
filename2=['C:\Data\Ming\Ming2_1500_1600_1pm_opt_pol (1).csv'];
filename1=['C:\Data\Ming\Ming1_1500_1600_High_Extinct (3).csv'];

data1=csvread(filename1,16,0);
data2=csvread(filename2,16,0);
data3=csvread(filename3,16,0);

figure(3)
hold off
plot(data1(:,1),data1(:,2),'red')
hold on
plot(data2(:,1),data2(:,2))
hold on
plot(data3(:,1),data3(:,2))


%% Plotting filter tarnsmission and reflection

filename4=['C:\Data\Filters\drive-download-20180905T144244Z-001\0828_2018_154812.csv'];
filename5=['C:\Data\Filters\drive-download-20180905T144244Z-001\0828_2018_154851.csv'];

data4=csvread(filename4,16,0);
figure(2)
hold off
plot(data4(:,1),data4(:,2),'red')
hold on
plot(data4(:,1),data4(:,3))

data5=csvread(filename5,16,0);

plot(data5(:,1),data5(:,2))
hold on
plot(data5(:,1),data5(:,3))

%% Plotting filter tarnsmission and reflection for 1552

filename4=['C:\Data\Filters\drive-download-20180905T144244Z-001\0828_2018_155212.csv'];
filename5=['C:\Data\Filters\drive-download-20180905T144244Z-001\0828_2018_155252.csv'];
filename6=['C:\Data\Filters\drive-download-20180905T144244Z-001\0828_2018_155293_2.csv'];

data4=csvread(filename4,16,0);
figure(2)
hold off
%plot(data4(:,1),data4(:,2),'red')
%hold on
plot(data4(:,1),data4(:,3))

data5=csvread(filename5,16,0);

%plot(data5(:,1),data5(:,2))
hold on
plot(data5(:,1),data5(:,3))

data6=csvread(filename6,16,0);

%plot(data6(:,1),data6(:,2))
hold on
plot(data6(:,1),data6(:,2))
