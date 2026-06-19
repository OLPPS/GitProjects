function [V_out, model] = pca_analytical_chain(Q, t, Cf, Rf, C_cr, R_cr, C_rc, R_rc)

% Step 1: Generate reference waveform (conv)
V_ref = convolution_chain(Q, t, Cf, Rf, C_cr, R_cr, C_rc, R_rc);

% Build training matrix (small perturbations to mimic variability)
Ntrain = 10;
Vmat = zeros(Ntrain,length(t));

for i=1:Ntrain
    scale = 1 + 0.05*randn;
    Vmat(i,:) = scale * V_ref;
end

% PCA
[coeff,score,~,~,explained,mu] = pca(Vmat);

% keep first few modes
K = 3;
basis = coeff(:,1:K);

% project original
alpha = (V_ref - mu) * basis;

% reconstruct
V_out = mu + alpha * basis';

model.basis = basis;
model.mu = mu;

end