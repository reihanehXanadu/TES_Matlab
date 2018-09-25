


%filename=['Chip_Passband_11300nW_umpATT1dB_s_i_20180915-1950'];
%filename=['Chip_Ring_1600nW_umpATT1dB_s_i_20180915-1952'];
%filename=['Chip_RING_460nW_umpATT4dB_s_i_20180915-1937'];
% filename=['Chip_Passband_5400nW_umpATT4dB_s_i_20180915-1933'];
% filename=['Chip_Ring_400nW_umpATT7dB_s_i_20180915-2000'];
% filename=['Chip_Passband_2400nW_umpATT7dB_s_i_20180915-1959'];
% filename=['Chip_Passband_5400nW_umpATT4dB_s_i_20180915-1933'];
 filename=['Test_pulses1'];
fid = fopen([filename,'.daq'], 'rb');
data = fread(fid,'uint16','ieee-le');

 data3D = double(reshape(data, S.numchannels2record, S.RecordLength, S.RecordsPerBuffer*S.buffersPerAcquisition ));

for i=1:2
subplot(2,1,i)
data2D=squeeze(data3D(i,:,1:1000));
plot(data2D)
end