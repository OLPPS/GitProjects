function V = analytical_chain(Q, t, Rvals, Cvals)

dt = t(2)-t(1);

% =============================
% Step 1: Current
% =============================
I = gradient(Q, dt);
I = I(:)';   % ensure row vector

% =============================
% Step 2: Fit exponentials
% =============================

% Use sensible time constants (seconds)
tau = [5e-9 20e-9 100e-9];   % << FIXED

% Detector window (same length as Q)
t_det = t(1:length(I));
I_det = I(1:length(I));

% ---- NORMALISE FOR STABILITY ----
scale = max(abs(I_det));
if scale == 0
    scale = 1;
end

I_scaled = I_det / scale;

% Build matrix
N = length(tau);
A = zeros(length(t_det), N);

for i = 1:N
    A(:,i) = exp(-t_det'/tau(i));
end

% Solve least squares
B_scaled = A \ I_scaled(:);

% Rescale coefficients
B = (B_scaled * scale)';

% =============================
% Step 3: CSA analytical response
% =============================

tau_f = Rvals(1) * Cvals(1);
Cf = Cvals(1);   % IMPORTANT scaling

V = zeros(size(t));

for i = 1:length(tau)

    td = tau(i);

    % Handle near-equal taus safely
    if abs(td - tau_f)/tau_f < 1e-3
        term = (t ./ tau_f) .* exp(-t/tau_f);
    else
        term = (exp(-t/tau_f) - exp(-t/td)) / (1/td - 1/tau_f);
    end

    % Include correct scaling
    V = V + (B(i)/Cf) * term;

end

% =============================
% Step 4: Apply shaping stages
% =============================

for k = 2:length(Rvals)

    tau_s = Rvals(k) * Cvals(k);
    V = apply_RC_analytic(V, t, tau_s);

end




I_fit = zeros(size(I_det));

for i = 1:length(tau)
    I_fit = I_fit + B(i)*exp(-t_det/tau(i));
end

figure;
plot(t_det, I_det, 'k', t_det, I_fit, '--r');
legend('True I','Fit');
title('Exponential Fit Check');


end