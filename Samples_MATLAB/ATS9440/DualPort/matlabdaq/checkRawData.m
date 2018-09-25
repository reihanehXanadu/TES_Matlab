clear all

fn = 'c:\daq\data\10-13-11_dual_TES_coincidence_detection_80MHz.daq';


fid = fopen(fn);
data = fread(fid,[2^21 1],'uint16','ieee-le');
fclose(fid);

plot(data(2:4:end),'r')
%plot(data,'r')
% hold on
% plot(data(2:4:end),'b')
% plot(data(3:4:end),'k')
% hold off