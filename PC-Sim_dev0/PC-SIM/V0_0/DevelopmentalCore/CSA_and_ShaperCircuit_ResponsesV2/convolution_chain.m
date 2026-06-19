function V = convolution_chain(Q, t, Cf, Rf, C_cr, R_cr, C_rc, R_rc)

dt = t(2)-t(1);
nt = length(t);

% Current
I = gradient(Q, dt);
I = I(:)';

% Pad
I_full = zeros(1,nt);
I_full(1:length(I)) = I;

%% CSA (filter implementation)

tau_f = Rf * Cf;

V_csa = zeros(1, nt);

for n = 2:nt
    V_csa(n) = V_csa(n-1) ...
             + (dt/Cf) * I_full(n) ...
             - (dt/tau_f) * V_csa(n-1);
end

figure(10)
hold off
plot(Q/max(Q))
hold on
plot(I/max(I))
legend('Q','I')
title('Input Q and I')

figure(11)
plot(V_csa)
title('CSA Voltage')

% =============================
% Step 3: CR stage (high-pass)
% =============================

tau_cr = R_cr * C_cr;
alpha = exp(-dt/tau_cr);

V_cr = zeros(1, nt);

for n = 2:nt
    V_cr(n) = alpha * V_cr(n-1)  + V_csa(n) - V_csa(n-1);
end

figure(12)
plot(V_cr)
title('CR output')
alpha

% =============================
% Step 4: RC stages (low-pass)
% =============================

V_rc = V_cr;

for k = 1:length(R_rc)

    tau_rc = R_rc(k) * C_rc(k);
    beta = exp(-dt/tau_rc);

    V_stage = zeros(1, nt);

    for n = 2:nt
        V_stage(n) = beta * V_stage(n-1) ...
                   + (1 - beta) * V_rc(n);
    end

    V_rc = V_stage;

end

figure(13)
plot(V_rc)
title('CSA-CR-RCn output')
beta


% =============================
% Final output
% =============================

V = V_rc;

figure(14)
plot(V)

end
