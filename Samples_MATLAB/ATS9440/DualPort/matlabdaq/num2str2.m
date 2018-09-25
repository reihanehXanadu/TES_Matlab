function str = num2str2(num,numround)
% function str = num2str2(num,numround)
% 
% same as num2str.m, but numround option
% 
% input 'numround' is optional
%		if input a #, rounds to number of decimal places
%		0 rounds to nearest integer


% round input numbers
if exist('numround','var')
	num = round2(num,numround);
else
	num = round(num);
end

str = num2str(num);

