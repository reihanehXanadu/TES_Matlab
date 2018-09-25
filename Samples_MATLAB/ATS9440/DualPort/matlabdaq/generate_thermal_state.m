function rho = generate_thermal_state(n, max_photon, theta)

r = asinh(sqrt(n));

m = 0:max_photon;
n = (0:max_photon).';


m =      1/cosh(r) .* (-1).^m .*exp(i*m*theta).*tanh(r).^m;
n = conj(1/cosh(r) .* (-1).^n .*exp(i*n*theta).*tanh(r).^n);

rho = repmat(m, max_photon+1, 1) .* repmat(n, 1, max_photon+1);

% Normalizing as a guard against cutoff errors...
rho = rho ./ trace(rho);
