function rho = generate_cat_state(alpha, theta, max_photon)
% Creates an equal superposition of two coherent states with opposite
% phases (namely, e^{i*theta}|-alpha>+|alpha>) and returns its density matrix (in
% the Fock basis, using photon numbers 0,1,...,max_photon).
% 
% input alpha = cat state amplitude
%       theta = cat state phase
%       max_photons = largest number of photons to consider
% output rho = array which is density matrix of the state

cat_vector = generate_cat_vector(alpha,theta,max_photon);
rho = cat_vector*cat_vector';

% check if rho is normalized.  I set a threshold at 10^(-6), but this may
% need to be adjusted depending on the accuracy needed.
tr = trace(rho);
if abs(tr)-1 > 10^(-6)
    warning('densitymatrix:tooSmall','Density matrix needs normalization.')
end

% ensure rho is normalized
rho = rho ./ tr;