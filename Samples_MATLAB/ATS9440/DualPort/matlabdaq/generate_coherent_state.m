function rho = generate_coherent_state(alpha, max_photon)
% Generates the density matrix for a coherent state alpha (in the Fock
% basis, [photon number]=0,1,...,max_photon).

m = 0:max_photon;
n = (0:max_photon).';

m = alpha .^ m ./ sqrt(factorial(m));
n = conj(alpha) .^ n ./ sqrt(factorial(n));

rho = exp(-1 * abs(alpha) ^ 2) .* repmat(m, max_photon+1, 1) .* repmat(n, 1, max_photon+1);

% Normalizing as a guard against cutoff errors...
rho = rho ./ trace(rho);
