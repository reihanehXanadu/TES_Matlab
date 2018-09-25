
function [g20,varargout] = G20(S,data,varargin)
[N_Avg,N2_Avg,N,N2,c] = Analyze_real_time(S,data)
maxpnumber=size(N,2);
P_n=N./repmat(sum(N,2),1,maxpnumber);
n=[0:maxpnumber-1];
g20=sum(n.*(n-1).*P_n,2)./sum(n.*P_n,2).^2;

end
