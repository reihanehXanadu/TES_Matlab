function [ns,Pl] = generate_squeezed_vacuum_with_loss(SdB,max_photon,eta)
% 
% SdB        = 0.5;
% max_photon = 20;
% eta        = 1;

S  = 10^(SdB/10);
rs = -1/2*log(S);
n_p = sinh(rs)^2

Ps = zeros(max_photon+1,1);
for i = 0:1:max_photon
    ns(i+1) = i;
    if (i/2 - floor(i/2)) == 0
        Ps(i+1) = factorial(i)/factorial(i/2)^2/cosh(rs)*(1/2*tanh(rs))^(i);
    else
        Ps(i+1) = 0;
    end
end

% introduce loss (1-eta)

Pl = zeros(max_photon+1,1);
for m = 0:1:max_photon
    for n = m:1:max_photon
        Pl(m+1) = Pl(m+1)+factorial(n)/(factorial(m)*factorial(n-m))*eta^m*(1-eta).^(n-m)*Ps(n+1);
    end
end


bar(ns,Ps,'r')
hold on
bar(ns,Pl)
hold off

