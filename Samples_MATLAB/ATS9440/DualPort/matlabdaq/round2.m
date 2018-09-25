function xout = round2(xin,digits)
% function xout = round2(xin,digits)
% rounds input number to specified digits
% round2(x,0) is the same as round(x)
% 
% ex. round2(3.145,2) = 3.15


xin = xin * power(10,digits);
xin = round(xin);
xout = xin / power(10,digits);