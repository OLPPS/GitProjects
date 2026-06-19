function V = convolution_chain(Q, t, Rvals, Cvals)

dt = t(2) - t(1);
nt = length(t);

% =============================
% Step 1: Compute current
% =============================
I = gradient(Q, dt);
I = I(:)';   % force row vector

% =============================
% Step 2: Zero pad to full time
% =============================
I_full = zeros(1, nt);
I_full(1:length(I)) = I;

% =============================
% Step 3: CSA stage
% =============================
tau_f = Rvals(1) * Cvals(1);
Cf    = Cvals(1);

h_csa = (1/Cf) * exp(-t / tau_f);

% =============================
% Step 4: CR stage (high-pass)
% =============================

% Discrete delta approximation
delta = zeros(1, nt);
delta(1) = 1/dt;   % discrete delta

% Use first shaping stage parameters for CR
tau_cr = Rvals(2) * Cvals(2);

h_cr = delta - (1/tau_cr) * exp(-t / tau_cr);

% =============================
% Step 5: RC stages (low-pass)
% =============================
h_total = conv(h_csa, h_cr, 'full');

for k = 2:length(Rvals)

    tau_rc = Rvals(k) * Cvals(k);

    h_rc = (1/tau_rc) * exp(-t / tau_rc);

    h_total = conv(h_total, h_rc, 'full');

end

% Truncate to original length
h_total = h_total(1:nt);

% Ensure row vector
h_total = h_total(:)';

% =============================
% Step 6: Apply convolution (FFT)
% =============================
V = real(ifft(fft(I_full) .* fft(h_total))) * dt;

end
