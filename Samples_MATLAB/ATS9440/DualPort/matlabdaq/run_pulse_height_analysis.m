clear all
d  = 'C:\daq\Data\';
d2 = [d '130304_WFH_Vaccum_PH_HOMdataStgPsn*'];
p  = dir(d2);

P.hist_bins = 200;
load hist_vector.txt

for chns = 1:2
    bins = linspace(hist_vector(chns,1),hist_vector(chns,2),P.hist_bins);
    bins3d(chns,:) = bins;
end

edges{1} = bins3d(1,:);
edges{2} = bins3d(2,:);
figure(1)
 for i = 1:length(p)
%for i = 3:
    load([d p(i).name])
    S = size(PHtemp2);
    PH1 = zeros(S(1),S(2)/2);
    PH2 = PH1;
    PH1(:,1:end) = PHtemp2(:,1:2:end);
    PH2(:,1:end) = PHtemp2(:,2:2:end);
    PH1 = reshape(PH1,1,S(1)*S(2)/2);
    PH2 = reshape(PH2,1,S(1)*S(2)/2);
    %dPH(i) = var(PH1-PH2);
    PH = [PH1; PH2;]';
    Y3 = hist3(PH,'Edges',edges);
    pcolor(Y3);
    title(sprintf('%s',p(i).name))
    shading flat;
    %pause(0.1)
    drawnow;
    mn1(i) = mean(PH1);
    mn2(i) = mean(PH2);
    v1(i)  = var(PH1);
    v2(i)  = var(PH2);
 end
%%
figure(5)

plot(mn1,'r')
