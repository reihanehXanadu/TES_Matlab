

 filename=['Chip_ring_DCDC_nW_pump6ATTdB_4loss_20180921-1532'];
fid = fopen([filename,'.daq'], 'rb');
%data = fread(fid,'uint16','ieee-le');
data = fread(fid,'uint16','ieee-le');

 data3D = double(reshape(data, S.numchannels2record, S.RecordLength, S.RecordsPerBuffer*S.buffersPerAcquisition ));

for i=1:2
subplot(2,1,i)
data2D=squeeze(data3D(i,:,1:1000));
plot(data2D)
end