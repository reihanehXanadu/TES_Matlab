data3D = double(reshape(data, S.numchannels2record, S.RecordLength, S.RecordsPerBuffer ));
data2D=squeeze(data3D(1,:,:));
for i=1:2
subplot(2,1,i)
data2D=squeeze(data3D(i,:,1:1000));
plot(data2D)
end