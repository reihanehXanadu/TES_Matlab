clear all
d  = 'C:\daq\Data\';
d2 = [d '130222LongScan_40degrees_NBAR_One_Mode_HOMdataStgPsn*'];
p  = dir(d2);

% P.hist_bins = 200;
% load hist_vector.txt
% 
% for chns = 1:2
%     bins = linspace(hist_vector(chns,1),hist_vector(chns,2),P.hist_bins);
%     bins3d(chns,:) = bins;
% end
% 
% edges{1} = bins3d(1,:);
% edges{2} = bins3d(2,:);

 for i = 3:length(p)
%for i = 3:4
    
    load([d p(i).name])
    S = size(nbar2);
    nb1 = zeros(S(1),S(2)/2);
    nb2 = nb1;
    nb1(:,1:end) = nbar2(:,1:2:end);
    nb2(:,1:end) = nbar2(:,2:2:end);
    nb1 = reshape(nb1,1,S(1)*S(2)/2);
    nb2 = reshape(nb2,1,S(1)*S(2)/2);
    %dPH(i) = var(PH1-PH2);
    nb = [nb1; nb2;]';
%    Y3 = hist3(nb,'Edges',edges);
 %   bar3(Y3)
%     pcolor(Y3);
%     title(sprintf('%s',p(i).name))
%     shading flat;
    %pause(0.1)
%     drawnow;
    mn1(i) = mean(nb1);
    mn2(i) = mean(nb2);
    v1(i)  = var(nb1);
    v2(i)  = var(nb2);
 end
%%
figure(5)

plot(mn1,'r')

save nbars.mat mn1 mn2

